"""
Fracture-Mechanics-Based Sectional Design (FMSD) workflow.

The proposed procedure (documented in docs/manual, Chapter "Proposed method"):

  Step 1  Material identification: matrix (fctm, GF) + residual fibre law from EN 14651 fR,j
          (MC2010 linear, EC2 Annex L, or inverse analysis with the cracked hinge).
  Step 2  Characteristic and design cohesive laws: sigma_d(w) = kO kG sigma_k(w) / gamma_SF
          (fibres) and fct,d = fctk/gamma_c (matrix), rebars with fyd, fcd for compression.
  Step 3  Cracked-hinge analysis of the critical section: M-theta, crack depth, crack opening at
          the rebar level / at the soffit, rebar stress.
  Step 4  ULS: M_Rd = max M on the design path before failure (crushing, rupture, w > wu);
          alternatively ECOV with mean/characteristic analyses (MC2010 7.11.3.3.2).
  Step 5  Ductility: ratio M_peak,post / M_cr and rotation capacity; snap-back indicator
          (brittle response when the post-cracking peak is lower than the cracking peak).
  Step 6  SLS: crack opening directly from the hinge at M_ser (mean laws, gamma = 1) compared
          with w_lim; stresses in steel.
  Step 7  Cross-check with code formulae (Annex L rigid-plastic, MC2010) and report.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
import numpy as np

from .materials import (CohesiveLaw, Rebar, Steel, compression_annexL, compression_ec2,
                        compression_parabola_rectangle, ecm_from_fcm, fctm_from_fck, gf_mc2010,
                        law_annexL, law_envelope, law_mc2010_frc_linear, law_mc2010_plain, tcm_bond,
                        characteristic_length)
from .hinge import CrackedHinge, Section
from . import codes
from .safety import ecov


@dataclass
class DesignInput:
    b: float = 300.0
    h: float = 500.0
    fck: float = 30.0
    fR1k: float = 0.0
    fR3k: float = 0.0
    bars: list = field(default_factory=list)  # list of (y_from_bottom, As, phi)
    fyk: float = 500.0
    fuk: float = 540.0
    eps_uk: float = 0.05
    kO: float = 1.0
    gamma_c: float = 1.5
    gamma_s: float = 1.15
    gamma_SF: float = 1.5
    law_model: str = "annexL"  # "annexL" or "mc2010"
    M_Ed: float = 0.0  # N mm
    M_ser: float = 0.0  # N mm
    w_lim: float = 0.3  # mm
    V_fR: float = 0.20  # assumed CoV of residual strengths, for mean values (ECOV)
    wu: float = 2.5
    n_layers: int = 300

    @property
    def fcm(self):
        return self.fck + 8.0


def cracking_and_post_peak(M, w, w_cr, Mcr0=None, w_window=0.2):
    """Cracking peak and post-cracking maximum of a moment-rotation curve.

    The cracking peak is the running maximum at the first point where M drops by more than 1 %
    while the soffit crack opening is still small (w <= w_window) and the running maximum does not
    exceed 2.2 f_ct W (upper bound of the modulus of rupture). If no such drop exists the response
    is monotonic through cracking and M_cr is M at w = w_cr.
    """
    M = np.asarray(M)
    w = np.asarray(w)
    runmax = np.maximum.accumulate(M)
    cond = (M < 0.99 * runmax) & (w <= w_window)
    if Mcr0 is not None:
        cond &= runmax <= 2.2 * Mcr0
    idx = np.where(cond)[0]
    if idx.size:
        j = int(idx[0])
        return float(runmax[j]), float(M[j:].max())
    return float(np.interp(w_cr, w, M)), float(M.max())


def _cohesive(inp: DesignInput, level: str, kG: float):
    fctm = fctm_from_fck(inp.fck)
    fctk = 0.7 * fctm
    Ecm = ecm_from_fcm(inp.fcm)
    GF = gf_mc2010(inp.fcm)
    fR1, fR3 = inp.fR1k, inp.fR3k
    if level == "mean":
        # lognormal residual strengths (consistent with the ECOV format): f_Rm = f_Rk exp(1.645 V)
        ft, fac = fctm, inp.kO * kG
        fR1, fR3 = fR1 * np.exp(1.645 * inp.V_fR), fR3 * np.exp(1.645 * inp.V_fR)
    elif level == "char":
        ft, fac = fctk, inp.kO * kG
    else:  # design
        ft, fac = fctk / inp.gamma_c, inp.kO * kG / inp.gamma_SF
    plain = law_mc2010_plain(ft, GF)
    if inp.fR1k <= 0:
        return plain, ft, Ecm, GF
    if inp.law_model == "mc2010":
        fib = law_mc2010_frc_linear(fR1, fR3, wu=inp.wu).scaled(fac)
    else:
        fib = law_annexL(fR1, fR3, ft, Ecm, kO=1.0, kG=1.0, wu=inp.wu,
                         include_matrix=False).scaled(fac)
    return law_envelope(plain, fib, w_max=inp.wu, name=f"{level} composite"), ft, Ecm, GF


def build_hinge(inp: DesignInput, level: str = "design", kG: float = 1.0) -> CrackedHinge:
    law, ft, Ecm, GF = _cohesive(inp, level, kG)
    fctm = fctm_from_fck(inp.fck)
    t0, t1 = tcm_bond(fctm if level == "mean" else 0.7 * fctm)
    sfrc = inp.fR1k > 0
    if level == "design":
        fcd = codes.eta_cc(inp.fck) * inp.fck / inp.gamma_c
        comp = compression_parabola_rectangle(fcd, sfrc=sfrc)
        steel = Steel(inp.fyk / inp.gamma_s, inp.fuk / inp.gamma_s, 200000.0, 0.9 * inp.eps_uk)
    elif level == "char":
        comp = compression_annexL(inp.fck, inp.fR1k, Ecm) if sfrc else compression_ec2(inp.fck, Ecm)
        steel = Steel(inp.fyk, inp.fuk, 200000.0, inp.eps_uk)
    else:
        comp = compression_annexL(inp.fcm, inp.fR1k, Ecm) if sfrc else compression_ec2(inp.fcm, Ecm)
        steel = Steel(1.1 * inp.fyk, 1.1 * inp.fuk, 200000.0, inp.eps_uk)
    rebars = [Rebar(y, As, phi, steel, t0, t1) for (y, As, phi) in inp.bars]
    return CrackedHinge(Section.rectangle(inp.b, inp.h), law, comp, Ecm, rebars=rebars,
                        n_layers=inp.n_layers, w_limit=inp.wu if sfrc and not rebars else None)


def run_design(inp: DesignInput, theta_max: float = 0.06, n: int = 160) -> dict:
    fctm = fctm_from_fck(inp.fck)
    Ecm = ecm_from_fcm(inp.fcm)
    GF = gf_mc2010(inp.fcm)
    out = {"input": asdict(inp)}
    # kappa_G: tension zone area ~ b (h - x); estimated with x = 0.2 h (iterated once below)
    Act = inp.b * 0.8 * inp.h * 1e-6
    kG = codes.kappa_G(Act) if inp.fR1k > 0 else 1.0
    out["kappa_G"] = kG
    res = {}
    for level in ("design", "char", "mean"):
        hg = build_hinge(inp, level, kG)
        r = hg.response(theta_max=theta_max, n=n)
        res[level] = r
    out["responses"] = res
    rd = res["design"]
    M = rd["M"]
    th = rd["theta"]
    Mcr0 = fctm_from_fck(inp.fck) * inp.b * inp.h ** 2 / 6.0
    M_cr, M_post = cracking_and_post_peak(M, rd["w_bottom"], GF / fctm_from_fck(inp.fck), Mcr0)
    M_Rd = float(M.max())
    out.update({"M_Rd_hinge": M_Rd, "M_cr_design": M_cr, "M_post_peak": M_post,
                "ductility_ratio": M_post / M_cr if M_cr > 0 else np.nan,
                "brittle_snapback": bool(M_post < M_cr), "failure_design": rd["failure"],
                "theta_u": float(th[-1])})
    # ECOV (MC2010 7.11.3.3.2)
    Rm, Rk = float(res["mean"]["M"].max()), float(res["char"]["M"].max())
    out["ecov"] = ecov(Rm, Rk)
    # SLS: crack opening at M_ser with mean laws
    rm = res["mean"]
    if inp.M_ser > 0:
        Mm = rm["M"]
        idx = np.where(Mm >= inp.M_ser)[0]
        if idx.size and rm["theta"][idx[0]] > 0:
            j = idx[0]
            if rm["w_s"].ndim == 2 and rm["w_s"].shape[1] > 0:
                w_ser = float(rm["w_s"][j].max())
            else:
                w_ser = float(rm["w_bottom"][j])
            out["w_ser"] = w_ser
        else:
            out["w_ser"] = np.nan
    # code cross-check (Annex L rigid plastic)
    fFtud = inp.kO * kG * 0.33 * inp.fR3k / inp.gamma_SF if inp.fR3k > 0 else 0.0
    layers = [(inp.h - y, As, inp.fyk / inp.gamma_s) for (y, As, phi) in inp.bars]
    out["annexL_uls"] = codes.uls_bending_annexL(inp.b, inp.h, inp.fck, fFtud, layers, inp.gamma_c)
    out["fracture"] = {"lch": characteristic_length(Ecm, GF, fctm), "GF": GF,
                       "h_over_lch": inp.h / characteristic_length(Ecm, GF, fctm)}
    if inp.M_Ed > 0:
        out["ULS_ok"] = M_Rd >= inp.M_Ed
    if inp.M_ser > 0 and "w_ser" in out:
        out["SLS_ok"] = bool(out["w_ser"] <= inp.w_lim)
    return out
