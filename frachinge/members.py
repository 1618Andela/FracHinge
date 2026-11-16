"""Member-level checks based on fracture mechanics (primary) with code formulas kept apart as
reference values.

Fracture-mechanics / mechanical models (primary results of FracHinge)
---------------------------------------------------------------------
* Shear of RC beams with or without stirrups: Compression Chord Capacity Model / Multi-Action Shear
  Model of Mari et al. (2014), Table 1. The cracked-web term uses the fracture energy G_F and the
  size effect of the compression chord follows Bazant's size-effect law.
* Shear of SFRC beams without stirrups: Multi-Action Shear Model of Mari et al. (2020), Eq. (26),
  with the residual tensile stress sigma_tu taken from the cohesive law (fracture-mechanics input).
* Combined bending and axial force: N-M interaction obtained with the cohesive cracked hinge.
* Slender columns: single-crack second-order model (central cohesive hinge + elastic segments, in the
  spirit of Ruiz Carmona et al. 2009); the ultimate load is the limit point of equilibrium.
* Punching of slabs: Critical Shear Crack Theory (Muttoni 2008) with the fibre-bridging term of
  Maya Duque et al. (2012); the bridging stress sigma(w) is read from the cohesive law.

Code formulas (reference only) are in :mod:`frachinge.codes`.
Units: N, mm, MPa.
"""
from __future__ import annotations

import math

import numpy as np

from .hinge import CrackedHinge, Section
from .materials import CohesiveLaw, fctm_from_fck

# ----------------------------------------------------------------------------------------------
# Shear: Multi-Action Shear Model (Mari et al. 2014, 2020)
# ----------------------------------------------------------------------------------------------


def masm_neutral_axis(rho_l: float, Es: float, Ec: float) -> float:
    """Relative flexural neutral-axis depth x/d of a cracked section (no concrete tension, linear
    compression), Mari et al. (2014) / (2020) Eq. (11)."""
    nr = Es / Ec * rho_l
    return nr * (-1.0 + math.sqrt(1.0 + 2.0 / nr))


def masm_zeta_2014(a_mm: float) -> float:
    """Size/slenderness factor of the compression chord, Mari et al. (2014) Eq. (8):
    zeta = 1.2 - 0.2 a >= 0.65 (a = shear span in metres)."""
    return max(1.2 - 0.2 * a_mm / 1000.0, 0.65)


def masm_xi_2020(d: float, a: float) -> float:
    """Combined size and slenderness factor, Mari et al. (2020) Eq. (10):
    xi = 2 / sqrt(1 + d0/200) * (d/a)^0.2 >= 0.45, d0 = d >= 100 mm."""
    d0 = max(d, 100.0)
    return max(2.0 / math.sqrt(1.0 + d0 / 200.0) * (d / a) ** 0.2, 0.45)


def masm_ec_2014(fcm: float) -> float:
    """Ec = 22 (fcm/10)^0.3 GPa, Mari et al. (2014) Eq. (11) [MPa]."""
    return 22000.0 * (fcm / 10.0) ** 0.3


def masm_gf_2014(fcm: float, dmax: float = 20.0) -> float:
    """Gf = 0.028 fcm^0.18 dmax^0.32 [N/mm], Mari et al. (2014) Eq. (12)."""
    return 0.028 * fcm ** 0.18 * dmax ** 0.32


def shear_masm_rc(b: float, d: float, a: float, rho_l: float, fcm: float, rho_w: float = 0.0,
                  fyw: float = 500.0, Es: float = 200000.0, dmax: float = 20.0,
                  GF: float | None = None, fct: float | None = None) -> dict:
    """Shear strength of an RC beam (with or without stirrups) by the Compression Chord Capacity
    Model (Mari et al. 2014, Table 1), mean values. Returns V_u [N] and its components.

    v_w = 167 fct/Ec (1 + 2 Ec Gf / (fct^2 d)); v_s = 0.85 rho_w fyw / fct;
    v_l = 0.25 x/d - 0.05 (only with stirrups); v_c = zeta [(0.88 + 0.70 v_s) x/d + 0.02].
    """
    Ec = masm_ec_2014(fcm)
    fct = fct if fct is not None else fctm_from_fck(min(fcm - 8.0, 60.0))
    Gf = GF if GF is not None else masm_gf_2014(fcm, dmax)
    xd = masm_neutral_axis(rho_l, Es, Ec)
    zeta = masm_zeta_2014(a)
    vw = 167.0 * fct / Ec * (1.0 + 2.0 * Ec * Gf / (fct ** 2 * d))
    vs = 0.85 * rho_w * fyw / fct
    vl = max(0.25 * xd - 0.05, 0.0) if vs > 0 else 0.0
    vc = zeta * ((0.88 + 0.70 * vs) * xd + 0.02)
    base = fct * b * d
    out = {"x_d": xd, "zeta": zeta, "V_c": vc * base, "V_w": vw * base, "V_l": vl * base,
           "V_s": vs * base, "V_u": (vc + vw + vl + vs) * base, "fct": fct, "Gf": Gf, "Ec": Ec}
    out["V_max_strut"] = strut_crushing(b, d, xd * d, fcm - 8.0, gamma_c=1.0)
    out["V_u"] = min(out["V_u"], out["V_max_strut"])
    return out


def shear_masm_sfrc(b: float, d: float, a: float, rho_l: float, fcm: float, sigma_tu: float,
                    Es: float = 200000.0, Ec: float | None = None, fct: float | None = None,
                    rho_w: float = 0.0, fyw: float = 500.0) -> dict:
    """Shear strength of an SFRC beam by the Multi-Action Shear Model (Mari et al. 2020, Eq. 26),
    mean values:  V_u = [xi (0.84 - 0.10 s) x/d + 0.08 + 1.10 s] fctm b d,  s = sigma_tu / fctm.

    sigma_tu is the average residual tensile stress of the SFRC; in FracHinge it is taken from the
    cohesive law (see :func:`sigma_tu_from_law`). If rho_w > 0 the stirrup term of Mari et al.
    (2014) is added by superposition (extension, not validated in the source).
    """
    fct = fct if fct is not None else fctm_from_fck(min(fcm - 8.0, 60.0))
    Ec = Ec if Ec is not None else masm_ec_2014(fcm)
    xd = masm_neutral_axis(rho_l, Es, Ec)
    xi = masm_xi_2020(d, a)
    s = min(sigma_tu / fct, 1.0)
    vu = xi * (0.84 - 0.10 * s) * xd + 0.08 + 1.10 * s
    Vu = vu * fct * b * d
    Vs = 0.85 * rho_w * fyw * b * d if rho_w > 0 else 0.0
    out = {"x_d": xd, "xi": xi, "sigma_tu_over_fct": s, "V_u_concrete_fibres": Vu, "V_s": Vs,
           "V_max_strut": strut_crushing(b, d, xd * d, fcm - 8.0, gamma_c=1.0)}
    out["V_u"] = min(Vu + Vs, out["V_max_strut"])
    out["stirrups_superposed"] = rho_w > 0
    return out


def shear_masm_design(b: float, d: float, a: float, rho_l: float, fck: float, fFtud: float = 0.0,
                      rho_w: float = 0.0, fywd: float = 435.0, gamma_c: float = 1.5,
                      Es: float = 200000.0, dmax: float = 20.0) -> dict:
    """Design shear resistance with the mechanical models (Mari et al. 2020 design proposal):
    the concrete tensile strength is expressed with the design compressive strength,
    fctd = 0.3 (fck/gamma_c)^(2/3) (fck <= 60 MPa), and the fibre residual stress is the design
    value fFtud of the cohesive law. Without fibres the CCCM of Mari et al. (2014) is used with the
    same fctd. The strut-crushing limit uses fcd. NOTE: the partial factors of these mechanical
    models have not been calibrated by a reliability analysis (stated by Mari et al. 2020)."""
    fcd = min(fck, 60.0) / gamma_c
    fctd = 0.3 * fcd ** (2.0 / 3.0)
    fcm = fck + 8.0
    if fFtud > 0:
        r = shear_masm_sfrc(b, d, a, rho_l, fcm, fFtud, Es=Es, fct=fctd, rho_w=rho_w, fyw=fywd)
        model = "MASM-SFRC (Mari et al. 2020, design)"
    else:
        r = shear_masm_rc(b, d, a, rho_l, fcm, rho_w=rho_w, fyw=fywd, Es=Es, dmax=dmax, fct=fctd)
        model = "CCCM-RC (Mari et al. 2014, design)"
    xd = r["x_d"]
    r["V_max_strut"] = strut_crushing(b, d, xd * d, fck, gamma_c=gamma_c)
    r["V_Rd"] = min(r["V_u"], r["V_max_strut"])
    r["model"] = model
    return r


def strut_crushing(b: float, d: float, x: float, fck: float, gamma_c: float = 1.5,
                   alpha_cw: float = 1.0) -> float:
    """Web-strut crushing used by Mari et al. (2014, 2020), Eq. (9): EC2 (2004) plasticity
    formula with the strut angle equal to that of the critical crack,
    cot(theta) = 0.85 d / (d - x) <= 2.5; z = 0.9 d; nu1 = 0.6 (1 - fck/250)."""
    cot = min(0.85 * d / max(d - x, 1e-6), 2.5)
    nu1 = 0.6 * (1.0 - fck / 250.0)
    return alpha_cw * b * 0.9 * d * nu1 * fck / gamma_c * cot / (1.0 + cot ** 2)


def sigma_tu_from_law(law: CohesiveLaw, w_from: float = 0.0, w_to: float = 1.5,
                      fibre_only: CohesiveLaw | None = None) -> float:
    """Average residual tensile stress of a cohesive law over [w_from, w_to] (fracture-mechanics
    input of the MASM). Pass the fibre-only law to exclude the matrix contribution."""
    lw = fibre_only if fibre_only is not None else law
    w = np.linspace(w_from, w_to, 401)
    return float(np.trapezoid(lw.sigma(w), w) / (w_to - w_from))


# ----------------------------------------------------------------------------------------------
# Combined bending and axial force with the cohesive hinge
# ----------------------------------------------------------------------------------------------


def interaction_hinge(section: Section, law: CohesiveLaw, comp, E: float, rebars=(),
                      N_values=None, n_layers: int = 300, theta_max: float = 0.05,
                      n: int = 120) -> dict:
    """N-M interaction of a section with the cohesive cracked hinge: for each axial force N
    (compression negative, N), the maximum moment of the M-theta response before the first failure
    criterion (crushing or bar rupture). Returns arrays N [N] and M [N mm]."""
    h = section.h
    if N_values is None:
        Ac = section.area
        Nc = -comp.fc * Ac
        N_values = np.linspace(0.95 * Nc, 0.0, 12)
    Ms = []
    for N in N_values:
        hg = CrackedHinge(section, law, comp, E, rebars=list(rebars), n_layers=n_layers)
        try:
            r = hg.response(theta_max=theta_max, n=n, N_ext=float(N))
            Ms.append(float(np.nanmax(r["M"])) if len(r["M"]) else np.nan)
        except Exception:
            Ms.append(np.nan)
    return {"N": np.asarray(N_values, float), "M": np.asarray(Ms, float)}


# ----------------------------------------------------------------------------------------------
# Slender columns: single cohesive hinge + elastic segments (second order)
# ----------------------------------------------------------------------------------------------


def column_buckling_hinge(section: Section, law: CohesiveLaw, comp, E: float, L: float,
                          e0: float, rebars=(), n_layers: int = 200, theta_max: float = 0.03,
                          n: int = 100, N_grid=None) -> dict:
    """Ultimate axial load of a pinned column of length L with first-order eccentricity e0 (mm),
    modelled as a central cohesive hinge plus elastic, uncracked end segments.

    Mid-height deflection: delta = theta L/4 + M L^2 / (8 E I) (hinge rotation + elastic bending
    of the segments under the mid-height moment); equilibrium M_hinge(theta, N) = N (e0 + delta).
    The ultimate load is the largest N for which this equation has a solution (limit point).
    Assumptions (uncracked segments, single crack) are those of a single-crack model and must be
    stated with the result."""
    I = section.inertia
    EI = E * I
    if N_grid is None:
        Ac = section.area
        N_grid = np.linspace(0.02, 0.95, 25) * comp.fc * Ac
    feasible = []
    for P in N_grid:
        hg = CrackedHinge(section, law, comp, E, rebars=list(rebars), n_layers=n_layers)
        try:
            r = hg.response(theta_max=theta_max, n=n, N_ext=-float(P), stop_on_failure=True)
        except Exception:
            feasible.append(False)
            continue
        M, th = r["M"], r["theta"]
        f = M * (1.0 - P * L ** 2 / (8.0 * EI)) - P * (e0 + th * L / 4.0)
        feasible.append(bool(np.nanmax(f) >= 0.0))
    feasible = np.asarray(feasible)
    if not feasible.any():
        return {"N_u": 0.0, "N_grid": np.asarray(N_grid), "feasible": feasible}
    k = int(np.where(feasible)[0].max())
    return {"N_u": float(N_grid[k]), "N_grid": np.asarray(N_grid), "feasible": feasible,
            "N_euler": math.pi ** 2 * EI / L ** 2}


# ----------------------------------------------------------------------------------------------
# Punching: Critical Shear Crack Theory with fibre bridging (Maya Duque et al. 2012)
# ----------------------------------------------------------------------------------------------


def punching_csct_sfrc(rc: float, d: float, rs: float, fc: float, fy: float, mR: float,
                       rq: float | None = None, dg: float = 16.0, Es: float = 200000.0,
                       fibre_law: CohesiveLaw | None = None, circular: bool = True) -> dict:
    """Punching strength of an interior slab-column connection by the Critical Shear Crack Theory
    (Muttoni 2008) with the fibre term of Maya Duque, Fernandez Ruiz, Muttoni & Foster (2012,
    Eng. Struct. 40), mean values (gamma = 1):

    failure criterion   V_Rc = 3/4 b0 d sqrt(fc) / (1 + 15 psi d / (dg0 + dg)),  dg0 = 16 mm
    fibre bridging      V_Rf = A_p sigma_f(w = psi d / 6)
    load-rotation       psi = 1.5 rs/d fy/Es (V / V_flex)^(3/2),  V_flex = 2 pi mR rs / (rq - rc)

    The bridging stress sigma_f(w) is read from the fibre cohesive law (FracHinge input) instead of
    the Variable Engagement Model used in the source. b0 is the control perimeter at d/2 from the
    column face, A_p the horizontal projection of a 45-degree failure surface.
    Units: N, mm, MPa; mR in N mm/mm (per unit width)."""
    rq = rq if rq is not None else rs  # load introduced at the contraflexure radius by default
    if circular:
        b0 = 2.0 * math.pi * (rc + 0.5 * d)
        Ap = math.pi * ((rc + d) ** 2 - rc ** 2)
    else:  # square column of side 2 rc
        b0 = 4.0 * 2.0 * rc + math.pi * d
        Ap = (2.0 * rc + 2.0 * d) ** 2 - (2.0 * rc) ** 2 - (4.0 - math.pi) * d ** 2
    Vflex = 2.0 * math.pi * mR * rs / max(rq - rc, 1e-9)

    def VR(psi):
        Vc = 0.75 * b0 * d * math.sqrt(fc) / (1.0 + 15.0 * psi * d / (16.0 + dg))
        Vf = float(fibre_law.sigma(np.array([psi * d / 6.0]))[0]) * Ap if fibre_law is not None else 0.0
        return Vc, Vf

    def VE(psi):  # load-rotation relationship (inverse form)
        return Vflex * (psi * d * Es / (1.5 * rs * fy)) ** (2.0 / 3.0)

    from scipy.optimize import brentq
    g = lambda p: sum(VR(p)) - VE(p)
    psi_hi = 1.5 * rs / d * fy / Es  # rotation at V = V_flex
    if g(psi_hi) > 0:  # punching not reached before flexural capacity
        Vc, Vf = VR(psi_hi)
        return {"psi": psi_hi, "V_R": Vflex, "V_Rc": Vc, "V_Rf": Vf, "V_flex": Vflex, "b0": b0,
                "A_p": Ap, "governing": "flexure (V_flex)"}
    psi = brentq(g, 1e-9, psi_hi)
    Vc, Vf = VR(psi)
    return {"psi": psi, "V_R": Vc + Vf, "V_Rc": Vc, "V_Rf": Vf, "V_flex": Vflex, "b0": b0, "A_p": Ap,
            "w_crack": psi * d / 6.0, "governing": "punching"}
