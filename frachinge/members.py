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
* Shear of beams and punching of slabs: unified approach of Frosch, Yu, Cusatis & Bazant (2017) with
  the energetic size-effect law of Bazant (ACI 446), plus fibre bridging of the cohesive law across the
  cracked part of the inclined crack, multiplied by the same size-effect law.

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
    Default fct = 0.30 min(fcm, 60)^(2/3), the convention used by the authors in their validation.

    v_w = 167 fct/Ec (1 + 2 Ec Gf / (fct^2 d)); v_s = 0.85 rho_w fyw / fct;
    v_l = 0.25 x/d - 0.05 (only with stirrups); v_c = zeta [(0.88 + 0.70 v_s) x/d + 0.02].
    """
    Ec = masm_ec_2014(fcm)
    # Mari et al. (2014), Eq. (9) with fcm in place of fck, limited to 60 MPa (as in their validation)
    fct = fct if fct is not None else 0.30 * min(fcm, 60.0) ** (2.0 / 3.0)
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


def size_factor_bazant(d: float, d0: float = 254.0) -> float:
    """Energetic size-effect law of Bazant (1984), in the form adopted by ACI 318-19 (22.5.5.1.3)
    for one-way shear: lambda_s = sqrt(2 / (1 + d/d0)) <= 1, d0 = 10 in = 254 mm."""
    return min(math.sqrt(2.0 / (1.0 + d / d0)), 1.0)


def shear_masm_sfrc(b: float, d: float, a: float, rho_l: float, fcm: float, sigma_tu: float,
                    Es: float = 200000.0, Ec: float | None = None, fct: float | None = None,
                    rho_w: float = 0.0, fyw: float = 500.0, size_effect: str | None = "bazant",
                    d0: float = 254.0, orientation: str = "aligned") -> dict:
    """Shear strength of an SFRC beam by the Multi-Action Shear Model (Mari et al. 2020, Eq. 26),
    mean values:  V_u = [xi (0.84 - 0.10 s) x/d + 0.08 + 1.10 s] fctm b d,  s = sigma_tu / fctm.

    sigma_tu is the average residual tensile stress of the SFRC; in FracHinge it is taken from the
    cohesive law (see :func:`sigma_tu_from_law`). If rho_w > 0 the stirrup term of Mari et al.
    (2014) is added by superposition (extension, not validated in the source).

    size_effect="bazant" (default, FracHinge 1.2): the fibre-bridging term 1.10 s is multiplied by
    Bazant's energetic size-effect factor lambda_s = sqrt(2/(1 + d/d0)) <= 1 (d0 = 254 mm, as in
    ACI 318-19). The original Eq. (26) has no size effect in this term, which made it unconservative
    for deep beams (Minelli et al. 2014). size_effect=None reproduces the original equation.

    orientation (FracHinge 1.3): sigma_tu comes from EN 14651 prisms, whose crack is normal to the
    preferential fibre direction (beam axis). The critical shear crack is inclined at theta to the axis
    (cot theta = 0.85 d / (d - x) <= 2.5, Mari et al.), so it is crossed by fewer aligned fibres.
    "aligned"  : fibres fully aligned with the beam axis -> bridging stress x sin(theta) (lower bound
                 of the fibre term, default);
    "isotropic": 3D random orientation -> no correction (upper bound).
    Experimental support: orientation numbers about 0.70 in prisms against 0.27-0.45 in cores
    (Cunha 2010), flow-induced alignment in SCFRC (Grunewald 2004; Ozyurt 2006).
    """
    fct = fct if fct is not None else fctm_from_fck(min(fcm - 8.0, 60.0))
    Ec = Ec if Ec is not None else masm_ec_2014(fcm)
    xd = masm_neutral_axis(rho_l, Es, Ec)
    xi = masm_xi_2020(d, a)
    cot = min(0.85 * d / max(d - xd * d, 1e-9), 2.5)
    eta = 1.0 / math.sqrt(1.0 + cot ** 2) if orientation == "aligned" else 1.0
    s = min(eta * sigma_tu / fct, 1.0)
    lam = size_factor_bazant(d, d0) if size_effect == "bazant" else 1.0
    vu = xi * (0.84 - 0.10 * s) * xd + 0.08 + 1.10 * s * lam
    Vu = vu * fct * b * d
    Vs = 0.85 * rho_w * fyw * b * d if rho_w > 0 else 0.0
    out = {"x_d": xd, "xi": xi, "lambda_F": lam, "eta_orientation": eta, "cot_theta_cr": cot, "sigma_tu_over_fct": s, "V_u_concrete_fibres": Vu,
           "V_fibres": 1.10 * s * lam * fct * b * d, "V_s": Vs,
           "V_max_strut": strut_crushing(b, d, xd * d, fcm - 8.0, gamma_c=1.0)}
    out["V_u"] = min(Vu + Vs, out["V_max_strut"])
    out["stirrups_superposed"] = rho_w > 0
    return out


def shear_masm_design(b: float, d: float, a: float, rho_l: float, fck: float, fFtud: float = 0.0,
                      rho_w: float = 0.0, fywd: float = 435.0, gamma_c: float = 1.5,
                      Es: float = 200000.0, dmax: float = 20.0, size_effect: str | None = "bazant",
                      orientation: str = "aligned") -> dict:
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
        r = shear_masm_sfrc(b, d, a, rho_l, fcm, fFtud, Es=Es, fct=fctd, rho_w=rho_w, fyw=fywd,
                            size_effect=size_effect, orientation=orientation)
        model = "MASM-SFRC (Mari et al. 2020, design)" + (" + Bazant size effect on fibres"
                                                         if size_effect == "bazant" else "") + (
                                                         ", aligned fibres" if orientation == "aligned" else "")
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
# Punching: CCCM for SFRC slabs (Mari, Spinella & Rodriguez 2023) -- second mechanical model
# ----------------------------------------------------------------------------------------------


def fctm_mc2010(fcm: float) -> float:
    """fctm of MC2010 as used by Mari et al. (2020, 2023): 0.3 fck^(2/3) (fck <= 50 MPa),
    2.12 ln(1 + 0.1 fcm) otherwise, fck = fcm - 8."""
    fck = fcm - 8.0
    return 0.3 * fck ** (2.0 / 3.0) if fck <= 50.0 else 2.12 * math.log(1.0 + 0.1 * fcm)


def sigma_tu_fibre_factor(Vf: float, Lf_df: float, fct: float, hooked: bool = True,
                          eta_l: float = 0.5) -> float:
    """Average residual tensile stress from the fibre factor (Lim et al.; Spinella et al.), Eq. (20) of
    Mari et al. (2023): sigma_tu / fctm = 2 eta0 eta_l beta_tau Vf Lf/df <= 1, eta0 = 0.405,
    beta_tau = 0.58 x (2.5 hooked, 1.2 straight). Vf as a fraction. Micromechanical reference value;
    FracHinge's primary route takes sigma_tu from the cohesive law (EN 14651)."""
    beta = 0.58 * (2.5 if hooked else 1.2)
    return min(2.0 * 0.405 * eta_l * beta * Vf * Lf_df, 1.0) * fct


def punching_cccm_sfrc(rcol: float, d: float, h: float, fcm: float, rho_l: float, sigma_tu: float,
                       fy: float | None = None, Es: float = 200000.0, size_reduction: bool = True,
                       fct: float | None = None) -> dict:
    """Punching strength of an interior SFRC slab-column connection without shear reinforcement,
    Compression Chord Capacity Model extended to SFRC by Mari, Spinella & Rodriguez (2023), Table 1,
    mean values. rcol: radius of the equivalent circular column (same perimeter: 2 b_c / pi for a
    square column). sigma_tu: residual tensile stress (from the cohesive law in FracHinge).

    Vcu = xi (0.425 + 1.125 x/d) fctm u_crit d >= Vcu,min = [xi (0.375 + 1.125 k_c) + 10/d0] fctm u_crit d,
      xi = 2 / sqrt(1 + d0/200) >= 0.45 (no slenderness term), u_crit = 2 pi (rcol + 0.5 d), k_c = x/d <= 0.2;
    Vcf = sigma_tu A_cr F_red, A_cr = pi (r_crack^2 - r_crit^2), r_crack = r_crit + d / (1 + 2/3 sigma_tu/fctm),
      F_red = 0.005 h <= 1 (reduction for thin slabs proposed by the authors);
    optional yielding limit Vy = 2 pi rho fy d^2 (1 - rho fy / (2 fck)).
    The slenderness term of xi is omitted because the slab span is not part of the published input."""
    fct = fct if fct is not None else fctm_mc2010(fcm)
    Ec = 22000.0 * (fcm / 10.0) ** 0.3
    ae_r = Es / Ec * rho_l
    xd = ae_r * (-1.0 + math.sqrt(1.0 + 2.0 / ae_r))
    d0 = max(d, 100.0)
    xi = max(2.0 / math.sqrt(1.0 + d0 / 200.0), 0.45)
    rcrit = rcol + 0.5 * d
    u = 2.0 * math.pi * rcrit
    Vcu = xi * (0.425 + 1.125 * xd) * fct * u * d
    Vmin = (xi * (0.375 + 1.125 * min(xd, 0.2)) + 10.0 / d0) * fct * u * d
    Vcu = max(Vcu, Vmin)
    s = min(sigma_tu / fct, 1.0)
    rcrack = rcrit + d / (1.0 + 2.0 / 3.0 * s)
    Acr = math.pi * (rcrack ** 2 - rcrit ** 2)
    Fred = min(0.005 * h, 1.0) if size_reduction else 1.0
    Vcf = s * fct * Acr * Fred
    out = {"x_d": xd, "xi": xi, "u_crit": u, "V_cu": Vcu, "V_cf": Vcf, "A_cr": Acr, "F_red": Fred,
           "V_u": Vcu + Vcf}
    if fy is not None:
        fck = fcm - 8.0
        out["V_y"] = 2.0 * math.pi * rho_l * fy * d ** 2 * (1.0 - rho_l * fy / (2.0 * fck))
    return out


# ----------------------------------------------------------------------------------------------
# Shear and punching: unified approach of Bazant et al. (energetic size-effect law)
# ----------------------------------------------------------------------------------------------

SQRT_PSI = 1.0 / 12.0416  # sqrt(f'c [psi]) expressed in MPa: k sqrt(psi) = (k/12.04) sqrt(MPa)


def neutral_axis_elastic(rho_l: float, fc: float, Es: float = 200000.0) -> float:
    """Elastic cracked-section neutral-axis ratio k = c/d = sqrt(2 rho n + (rho n)^2) - rho n with
    n = Es/Ec, Ec = 4700 sqrt(fc) (Frosch et al. 2017, Eq. 3)."""
    rn = rho_l * Es / (4700.0 * math.sqrt(fc))
    return math.sqrt(2.0 * rn + rn * rn) - rn


def size_factor_unified(d: float, stirrups: bool = False) -> float:
    """Size-effect factor of the unified approach (Frosch, Yu, Cusatis & Bazant 2017, Eqs. 6-7), based
    on the energetic size-effect law proposed by ACI Committee 446:
    gamma_d = 1.4 / sqrt(1 + d/d0) <= 1, d0 = 10 in (254 mm) without minimum stirrups and
    100 in (2540 mm) with A_v >= A_v,min; gamma_d = 1 for d < 10 in (or d < 100 in with stirrups)."""
    d0 = 2540.0 if stirrups else 254.0
    if d < d0:
        return 1.0
    return min(1.4 / math.sqrt(1.0 + d / d0), 1.0)


def _fibre_bridging(b: float, d: float, c: float, sigma_tu: float, orientation: str, cot_theta: float = 1.0):
    """Vertical resultant of the cohesive bridging stress on the cracked part (d - c) of an inclined
    crack: V_F = eta_or lambda_s sigma_tu b (d - c) cot(theta) (equilibrium of the crack face,
    theta = 45 deg). lambda_s: Bazant energetic size-effect law (d0 = 254 mm); eta_or = sin(theta)
    for fibres aligned with the member axis (or lying in the slab plane), 1 for isotropic fibres."""
    lam = size_factor_bazant(d, 254.0)
    eta = 1.0 / math.sqrt(1.0 + cot_theta ** 2) if orientation == "aligned" else 1.0
    return eta * lam * sigma_tu * b * max(d - c, 0.0) * cot_theta, lam, eta


def shear_bazant_unified(b: float, d: float, rho_l: float, fc: float, sigma_tu: float = 0.0,
                         rho_w: float = 0.0, fyw: float = 500.0, stirrups_min: bool = False,
                         orientation: str = "aligned", Es: float = 200000.0, phi: float = 1.0) -> dict:
    """Shear strength of a beam (RC, SFRC or both) with the unified approach of Frosch, Yu, Cusatis &
    Bazant (2017): V_c = 5 sqrt(f'c) b_w c gamma_d (psi, in) = 0.415 sqrt(fc) b c gamma_d (MPa, mm),
    c = elastic cracked neutral-axis depth, gamma_d from the energetic size-effect law; stirrups
    V_s = A_v f_y d / s; upper limit V <= 5 V_c (Eq. 9). FracHinge adds the fibre bridging term of
    the cohesive law (see _fibre_bridging). The concrete term is a nominal (lower-bound) expression;
    phi = 0.75 gives the ACI design format."""
    k = neutral_axis_elastic(rho_l, fc, Es)
    c = k * d
    gd = size_factor_unified(d, stirrups_min)
    Vc = 5.0 * SQRT_PSI * math.sqrt(fc) * b * c * gd
    VF, lam, eta = _fibre_bridging(b, d, c, sigma_tu, orientation)
    Vs = rho_w * b * fyw * d
    Vn = min(Vc + VF + Vs, 5.0 * (Vc + VF))
    return {"c": c, "k": k, "gamma_d": gd, "lambda_F": lam, "eta_orientation": eta, "V_c": Vc,
            "V_F": VF, "V_s": Vs, "V_n": Vn, "V_Rd": phi * Vn,
            "model": "Unified approach (Frosch, Yu, Cusatis & Bazant 2017) + cohesive fibre bridging"}


def punching_bazant_unified(b0: float, d: float, rho_l: float, fc: float, sigma_tu: float = 0.0,
                            orientation: str = "aligned", Es: float = 200000.0, phi: float = 1.0) -> dict:
    """Two-way (punching) shear strength with the unified approach: V_c = 10 sqrt(f'c) b0 c (psi, in)
    = 0.83 sqrt(fc) b0 c (MPa, mm) (Frosch et al. 2017; ACI 440.1R-15), b0 = control perimeter at d/2.
    FracHinge multiplies it by the energetic size-effect factor gamma_d (punching size effect verified
    by Bazant & Cao 1987, ACI 446 1992) and adds cohesive fibre bridging on b0 (d - c)."""
    k = neutral_axis_elastic(rho_l, fc, Es)
    c = k * d
    gd = size_factor_unified(d, False)
    Vc = 10.0 * SQRT_PSI * math.sqrt(fc) * b0 * c * gd
    VF, lam, eta = _fibre_bridging(b0, d, c, sigma_tu, orientation)
    Vn = Vc + VF
    return {"c": c, "gamma_d": gd, "lambda_F": lam, "eta_orientation": eta, "V_c": Vc, "V_F": VF,
            "V_n": Vn, "V_Rd": phi * Vn}
