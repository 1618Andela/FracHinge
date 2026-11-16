"""
Code-based (normative) formulations used as reference in FracHinge.

  * fib Model Code 2010 (fib, 2013): classification (5.6-1), (5.6-2), (5.6-3); FRC shear (7.7-5);
    FRC crack width (7.7-22)-(7.7-24); flexural tensile strength (5.1-8).
  * EN 1992-1-1:2023 (FprEN working draft 2022-05-20 used here): Annex L (L.1)-(L.6), (L.19),
    (L.20), (L.26)-(L.28), (L.29); main text (8.11), (8.16), 5.1.6 (eta_cc).
  * ACI 544.4R-88: Henager and Doherty (1976) flexural strength of beams with bars and fibres,
    Eqs. (3-3)-(3-6).

All functions use N, mm, MPa.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

# --------------------------------------------------------------------------- classification

STRENGTH_SERIES_MC2010 = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
STRENGTH_SERIES_EC2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]


def classify_mc2010(fR1k: float, fR3k: float) -> str:
    """MC2010 5.6.3: strength interval + letter from fR3k/fR1k, Eq. (5.6-1)."""
    cls = None
    for v in STRENGTH_SERIES_MC2010:
        if fR1k >= v:
            cls = v
    if cls is None:
        return "below 1.0 (not classifiable)"
    r = fR3k / fR1k
    letter = ("a" if 0.5 <= r < 0.7 else "b" if 0.7 <= r < 0.9 else "c" if 0.9 <= r < 1.1
              else "d" if 1.1 <= r < 1.3 else "e" if r >= 1.3 else "-")
    return f"{cls:g}{letter}"


def classify_ec2(fR1k: float, fR3k: float) -> str:
    """EN 1992-1-1:2023 Table L.2: class N with fR1k >= N; letter from fR3k >= {0.5..1.3} N."""
    N = None
    for v in STRENGTH_SERIES_EC2:
        if fR1k >= v:
            N = v
    if N is None:
        return "below 1.0 (not classifiable)"
    letter = "-"
    for f, l in zip([0.5, 0.7, 0.9, 1.1, 1.3], "abcde"):
        if fR3k >= f * N:
            letter = l
    return f"{N:g} {letter}"


def mc2010_substitution_ok(fR1k: float, fR3k: float, fLk: float) -> bool:
    """MC2010 Eqs. (5.6-2), (5.6-3): fR1k/fLk > 0.4 and fR3k/fR1k > 0.5."""
    return (fR1k / fLk > 0.4) and (fR3k / fR1k > 0.5)


def annexL_ductility_ok(fR1k: float, fctk005: float) -> bool:
    """Annex L Formula (L.1): fR1k / fctk,0.05 >= 0.5."""
    return fR1k / fctk005 >= 0.5


def kappa_G(Act_m2: float, kmax07: bool = False) -> float:
    """Annex L Formula (L.6): kappa_G = 1 + 0.5 Act (Act in m^2) <= 1.5 (1.3 if k_max = 0.7)."""
    return min(1.0 + 0.5 * Act_m2, 1.3 if kmax07 else 1.5)


def gamma_SF(situation: str = "persistent") -> float:
    """Annex L Table L.1 (NDP)."""
    return {"persistent": 1.5, "transient": 1.5, "accidental": 1.2, "sls": 1.0}[situation]


def eta_cc(fck: float) -> float:
    """EN 1992-1-1:2023 5.1.6: eta_cc = (40/fck)^(1/3) <= 1."""
    return min((40.0 / fck) ** (1.0 / 3.0), 1.0)


def fct_flexural_mc2010(fctm: float, h: float) -> float:
    """MC2010 Eq. (5.1-8): fctm = alpha_fl fctm,fl with alpha_fl = 0.06 h^0.7 / (1 + 0.06 h^0.7)."""
    a = 0.06 * h ** 0.7 / (1.0 + 0.06 * h ** 0.7)
    return fctm / a


# --------------------------------------------------------------------------- ULS bending

def uls_bending_annexL(b: float, h: float, fck: float, fFtud: float = 0.0, layers=(),
                       gamma_c: float = 1.5, N_Ed: float = 0.0, sfrc: bool | None = None,
                       eps_Ftu: float = 0.02, eps_ud: float = 0.0225, Es: float = 200000.0,
                       tension_block: str = "rigid-plastic", fFt1d: float | None = None,
                       fFt3d: float | None = None, lcs: float = 125.0, n_strips: int = 400):
    """Design moment resistance of a rectangular (S)FRC section (EN 1992-1-1:2023 8.1 + Annex L.8.1).

    layers: iterable of (d_i [mm from top], As_i [mm^2], fyd_i [MPa]).
    Compression: parabola-rectangle with eps_c2/eps_cu = 0.002/0.0035 (plain) or 0.0025/0.006 (SFRC).
    Tension in SFRC: rigid-plastic block fFtud over the whole tension zone (Figure L.2a) or
    bilinear (Figure L.2b) with fFt1d, fFt3d. Steel: elastic - perfectly plastic.
    Returns dict(M_Rd [N mm], x [mm], eps_top, eps_bot, governing).
    """
    if sfrc is None:
        sfrc = fFtud > 0
    ec2, ecu = (0.0025, 0.006) if sfrc else (0.002, 0.0035)
    fcd = eta_cc(fck) * fck / gamma_c
    ys = (np.arange(n_strips) + 0.5) * h / n_strips  # depth from top
    dA = b * h / n_strips

    def strains(x):
        # pivot selection: top crushing, or bottom tensile limit (fibre strain / steel strain)
        k1 = ecu / x
        lim = []
        if fFtud > 0 or (fFt1d is not None):
            lim.append(eps_Ftu / max(h - x, 1e-9))
        for d, As, fyd in layers:
            lim.append(eps_ud / max(d - x, 1e-9))
        k = min([k1] + lim)
        return k

    def forces(x):
        k = strains(x)
        eps = k * (ys - x)  # tension positive, below the neutral axis
        sig_c = np.where(eps < 0, -fcd * (1.0 - (1.0 - np.clip(-eps / ec2, 0, 1)) ** 2), 0.0)
        if tension_block == "rigid-plastic":
            sig_f = np.where(eps > 0, fFtud, 0.0)
        else:  # bilinear (Figure L.2b), strains converted with lcs
            e1 = 0.5 / lcs
            e3 = eps_Ftu
            sig_f = np.where(eps > 0, np.interp(eps, [0.0, e1, e3], [fFt1d, fFt1d, fFt3d]), 0.0)
        N = np.sum((sig_c + sig_f) * dA)
        M = np.sum((sig_c + sig_f) * dA * (h / 2.0 - ys))
        for d, As, fyd in layers:
            es = k * (d - x)
            ss = np.clip(Es * es, -fyd, fyd)
            N += As * ss
            M += As * ss * (h / 2.0 - d)
        return N, -M, k

    g = lambda x: forces(x)[0] - N_Ed
    x = brentq(g, 1e-3 * h, 0.999 * h, xtol=1e-6)
    N, M, k = forces(x)
    gov = "concrete crushing" if abs(k - ecu / x) < 1e-12 else "tensile strain limit"
    return {"M_Rd": abs(M), "x": x, "eps_top": -k * x, "eps_bot": k * (h - x), "governing": gov, "fcd": fcd}


# --------------------------------------------------------------------------- shear

def ddg(D_lower: float = 16.0, fck: float = 30.0) -> float:
    """EN 1992-1-1:2023 8.2.1(4): ddg = 16 + D_lower <= 40 mm (fck <= 60 MPa)."""
    if fck <= 60:
        return min(16.0 + D_lower, 40.0)
    return min(16.0 + D_lower * (60.0 / fck) ** 2, 40.0)


def tau_Rdc_min(fck, fyd, dv, D_lower=16.0, gamma_V=1.4):
    """EN 1992-1-1:2023 Formula (8.11)."""
    return 11.0 / gamma_V * np.sqrt(fck / fyd * ddg(D_lower, fck) / dv)


def tau_Rdc(rho_l, fck, dv, D_lower=16.0, gamma_V=1.4, fyd=435.0):
    """EN 1992-1-1:2023 Formula (8.16) (members without shear reinforcement)."""
    t = 0.66 / gamma_V * (100.0 * rho_l * fck * ddg(D_lower, fck) / dv) ** (1.0 / 3.0)
    return max(t, tau_Rdc_min(fck, fyd, dv, D_lower, gamma_V))


def tau_RdcF_annexL(rho_l, fck, dv, fFtud, fFtuk, D_lower=16.0, gamma_V=1.4, fyd=435.0):
    """Annex L Formulae (L.19), (L.20): tau_Rd,cF = eta_cF tau_Rd,c + eta_F fFtud."""
    eta_cF = min(max(1.2 - 0.5 * fFtuk, 0.4), 1.0)
    t = 0.66 / gamma_V * (100.0 * rho_l * fck * ddg(D_lower, fck) / dv) ** (1.0 / 3.0)
    tmin = tau_Rdc_min(fck, fyd, dv, D_lower, gamma_V)
    return max(eta_cF * t + fFtud, eta_cF * tmin + fFtud), eta_cF


def VRdF_mc2010(bw, d, rho_l, fck, fFtuk, fctk, sigma_cp=0.0, gamma_c=1.5):
    """MC2010 Eq. (7.7-5): V_Rd,F for members with longitudinal bars and without stirrups [N].

    fFtuk is the characteristic ultimate residual strength for wu = 1.5 mm (Eq. 5.6-6).
    """
    k = min(1.0 + np.sqrt(200.0 / d), 2.0)
    rho = min(rho_l, 0.02)
    V = (0.18 / gamma_c * k * (100.0 * rho * (1.0 + 7.5 * fFtuk / fctk) * fck) ** (1.0 / 3.0)
         + 0.15 * sigma_cp) * bw * d
    vmin = 0.035 * k ** 1.5 * np.sqrt(fck)
    return max(V, (vmin + 0.15 * sigma_cp) * bw * d)


# --------------------------------------------------------------------------- SLS cracking

def crack_spacing_annexL(c, phi, rho_eff, fFt1_ef, fctm, h, x, s_bars, kfl=1.0, kb=0.9):
    """Annex L (L.26)-(L.28): mean crack spacing of members with bars and fibres."""
    if s_bars <= 10.0 * phi:
        alpha_f = min(fFt1_ef / fctm, 1.0)
        return 1.5 * c + kfl * kb / 7.2 * phi / rho_eff * (1.0 - alpha_f), 1.7
    return h - x, 1.3


def crack_width_mc2010_frc(sigma_s, phi, rho_eff, c, fctm, fFtsk, Es=200000.0, Ec=33000.0,
                           beta=0.6, tau_factor=1.8, k=1.0):
    """MC2010 Eqs. (7.7-22)-(7.7-24) (short-term, stabilised or crack formation, no shrinkage)."""
    fFtsm = fFtsk / 0.7
    tau = tau_factor * fctm
    ae = Es / Ec
    ls = k * c + 0.25 * (fctm - fFtsm) / tau * phi / rho_eff
    ls = max(ls, k * c)
    sig_sr = max(fctm - fFtsm, 0.0) * (1.0 + ae * rho_eff) / rho_eff
    eps = max(sigma_s - beta * sig_sr, 0.0) / Es
    return 2.0 * ls * eps, ls


# --------------------------------------------------------------------------- minimum reinforcement

def as_min_bending(b, h, d, fctm, fyk, fck=30.0, fFtu_ef=0.0):
    """Annex L (L.29) / EC2 (12.1): smallest As with M_R,min >= M_cr (characteristic values)."""
    Mcr = fctm * b * h ** 2 / 6.0

    def MR(As):
        return uls_bending_annexL(b, h, fck, fFtud=fFtu_ef, sfrc=fFtu_ef > 0,
                                  layers=[(d, As, fyk)], gamma_c=1.0)["M_Rd"]

    if fFtu_ef > 0 and MR(1e-6) >= Mcr:
        return 0.0, Mcr
    return brentq(lambda A: MR(A) - Mcr, 1e-6, 0.05 * b * h), Mcr


# --------------------------------------------------------------------------- ACI 544.4R

def henager_doherty(b, h, d, As, fy, fc, lf, df, Vf_percent, Fbe=1.0, tau_d=2.3, Es=200000.0,
                    eps_cu=0.003):
    """ACI 544.4R-88 Eqs. (3-3)-(3-6): nominal moment of a beam with bars and steel fibres [N mm].

    sigma_t = 0.00772 (l/d) rho_f Fbe  [MPa] (rho_f in percent);
    e = [eps_s(fibres) + 0.003] c / 0.003 with eps_s(fibres) = sigma_f / Es, where sigma_f is the
    fibre stress at pull-out, here taken as 2 tau_d (l/d) (half embedment length, tau_d = 2.3 MPa).
    """
    sigma_t = 0.00772 * lf / df * Vf_percent * Fbe
    sigma_f = 2.0 * tau_d * lf / df
    eps_f = sigma_f / Es
    beta1 = 0.85 if fc <= 28 else max(0.85 - 0.05 * (fc - 28) / 7.0, 0.65)

    def resid(c):
        a = beta1 * c
        e = (eps_f + eps_cu) * c / eps_cu
        C = 0.85 * fc * a * b
        T = As * fy + sigma_t * b * max(h - e, 0.0)
        return C - T

    c = brentq(resid, 1e-3, h)
    a = beta1 * c
    e = (eps_f + eps_cu) * c / eps_cu
    Mn = As * fy * (d - a / 2.0) + sigma_t * b * max(h - e, 0.0) * (h / 2.0 + e / 2.0 - a / 2.0)
    return {"Mn": Mn, "c": c, "a": a, "e": e, "sigma_t": sigma_t}
