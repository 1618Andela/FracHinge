"""
Constitutive laws used by FracHinge.

Units: N, mm, MPa (N/mm^2). Crack openings w in mm. Fracture energy in N/mm.

References (see docs/bib/references.bib):
  - fib Model Code 2010 (fib, 2013): Eqs. (5.1-9), (5.1-31), (5.1-32), (5.6-5), (5.6-6).
  - EN 1992-1-1:2023 (FprEN 2022 working draft), Formulae (5.6)-(5.10), Annex L (L.2), (L.3),
    (L.7)-(L.14).
  - Marti, Alvarez, Kaufmann, Sigrist (1998) Tension Chord Model, SEI 8(4):287-298.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

# --------------------------------------------------------------------------------------
# Concrete: basic properties (fib MC2010 / EN 1992-1-1)
# --------------------------------------------------------------------------------------

def fctm_from_fck(fck: float) -> float:
    """Mean axial tensile strength, MC2010 Eq. (5.1-3a/b) [MPa]."""
    if fck <= 50.0:
        return 0.3 * fck ** (2.0 / 3.0)
    return 2.12 * np.log(1.0 + 0.1 * (fck + 8.0))


def ecm_from_fcm(fcm: float, alpha_e: float = 1.0) -> float:
    """Secant modulus, EN 1992-1-1:2023 Formula (5.1): Ecm = 22000 (fcm/10)^0.3 [MPa]."""
    return alpha_e * 22000.0 * (fcm / 10.0) ** 0.3


def gf_mc2010(fcm: float) -> float:
    """Fracture energy of plain concrete, MC2010 Eq. (5.1-9): GF = 73 fcm^0.18 [N/m] -> N/mm."""
    return 73.0 * fcm ** 0.18 / 1000.0


def characteristic_length(E: float, GF: float, ft: float) -> float:
    """Hillerborg's characteristic length lch = E GF / ft^2 [mm]."""
    return E * GF / ft ** 2


# --------------------------------------------------------------------------------------
# Cohesive (stress - crack opening) laws
# --------------------------------------------------------------------------------------

@dataclass
class CohesiveLaw:
    """Piecewise-linear softening/hardening law sigma(w).

    Nodes (w_i, s_i) with w_0 = 0. For w beyond the last node the stress is zero
    (complete separation / fibre pull-out) unless ``tail`` is set to "constant".
    """

    w: np.ndarray
    s: np.ndarray
    name: str = "cohesive"
    tail: str = "zero"  # "zero" or "constant"

    def __post_init__(self):
        self.w = np.asarray(self.w, dtype=float)
        self.s = np.asarray(self.s, dtype=float)
        if self.w[0] != 0.0:
            raise ValueError("First node must be at w = 0")
        if np.any(np.diff(self.w) <= 0):
            raise ValueError("Crack-opening nodes must be strictly increasing")
        if np.any(self.s < 0):
            raise ValueError("Cohesive stresses must be non-negative")

    @property
    def ft(self) -> float:
        """Stress at w = 0 (tensile strength of the law)."""
        return float(self.s[0])

    def sigma(self, w):
        w = np.asarray(w, dtype=float)
        out = np.interp(w, self.w, self.s)
        beyond = w > self.w[-1]
        if self.tail == "zero":
            out = np.where(beyond, 0.0, out)
        else:
            out = np.where(beyond, self.s[-1], out)
        return np.where(w < 0, self.s[0], out)

    def energy(self, w_max: float | None = None) -> float:
        """Area under sigma(w) up to w_max [N/mm]."""
        if w_max is None:
            w_max = self.w[-1]
        ww = np.concatenate([self.w[self.w < w_max], [w_max]])
        return float(np.trapezoid(self.sigma(ww), ww))

    def segments(self):
        """Return arrays (w_lo, w_hi, a, b) with sigma = a + b w on each segment."""
        w0, w1 = self.w[:-1], self.w[1:]
        s0, s1 = self.s[:-1], self.s[1:]
        b = (s1 - s0) / (w1 - w0)
        a = s0 - b * w0
        return w0, w1, a, b

    def scaled(self, factor: float, name: str | None = None) -> "CohesiveLaw":
        return CohesiveLaw(self.w.copy(), self.s * factor, name or f"{self.name} x{factor:g}", self.tail)


def law_mc2010_plain(fctm: float, GF: float) -> CohesiveLaw:
    """Bilinear softening of plain concrete, MC2010 Eqs. (5.1-31)-(5.1-32).

    w1 = GF/fctm (sigma = 0.2 fctm), wc = 5 GF/fctm (sigma = 0).
    """
    w1 = GF / fctm
    wc = 5.0 * GF / fctm
    return CohesiveLaw([0.0, w1, wc], [fctm, 0.2 * fctm, 0.0], name="MC2010 plain (bilinear)")


def law_linear_softening(ft: float, GF: float) -> CohesiveLaw:
    """Linear softening law with the same fracture energy (wc = 2 GF/ft)."""
    return CohesiveLaw([0.0, 2.0 * GF / ft], [ft, 0.0], name="linear softening")


def mc2010_residual(fR1: float, fR3: float, wu: float = 2.5, cmod3: float = 2.5):
    """Residual tensile strengths of the MC2010 linear model, Eqs. (5.6-5), (5.6-6)."""
    fFts = 0.45 * fR1
    fFtu = fFts - wu / cmod3 * (fFts - 0.5 * fR3 + 0.2 * fR1)
    return fFts, max(fFtu, 0.0)


def law_mc2010_frc_linear(fR1: float, fR3: float, wu: float = 2.5, w_end: float | None = None) -> CohesiveLaw:
    """MC2010 linear post-cracking law for FRC (Fig. 5.6-7), sigma(0) = fFts, sigma(wu) = fFtu.

    ``w_end`` (> wu) optionally prolongs the same straight line (clipped at zero), which is
    useful to simulate tests beyond wu; by default the law stops at wu.
    """
    fFts, fFtu = mc2010_residual(fR1, fR3, wu)
    if w_end is None or w_end <= wu:
        return CohesiveLaw([0.0, wu], [fFts, fFtu], name="MC2010 FRC linear")
    slope = (fFtu - fFts) / wu
    s_end = max(fFts + slope * w_end, 0.0)
    return CohesiveLaw([0.0, wu, w_end], [fFts, fFtu, s_end], name="MC2010 FRC linear (prolonged)")


def law_mc2010_rigid_plastic(fR3: float, wu: float = 2.5) -> CohesiveLaw:
    """MC2010 rigid-plastic model, Eq. (5.6-4): fFtu = fR3/3 up to wu."""
    return CohesiveLaw([0.0, wu], [fR3 / 3.0, fR3 / 3.0], name="MC2010 rigid-plastic")


def annexL_residual(fR1k: float, fR3k: float, kO: float = 1.0, kG: float = 1.0):
    """EN 1992-1-1:2023 Annex L, Formulae (L.2), (L.3), (L.7), (L.8)."""
    return {
        "fFts_ef": kO * kG * 0.37 * fR1k,
        "fFtu_ef": kO * kG * 0.33 * fR3k,
        "fFt1_ef": kO * kG * 0.37 * fR1k,
        "fFt3_ef": kO * kG * max(0.57 * fR3k - 0.26 * fR1k, 0.0),
    }


def law_annexL(fR1k: float, fR3k: float, fctm: float, Ecm: float, lcs: float = 125.0,
               kO: float = 1.0, kG: float = 1.0, wu: float = 2.5, include_matrix: bool = True) -> CohesiveLaw:
    """Annex L constitutive law (Figure L.1) rewritten in terms of crack opening w = eps*lcs.

    The residual branch is the straight line through (w = 0.5 mm, fFt1,ef) and
    (w = 2.5 mm, fFt3,ef), cut at wu (eps_Ftu = wu/lcs). The elastic peak fctm is reached at
    eps_ctm; the drop to the residual branch ends at eps_F,0 = 2 eps_ctm, i.e. w0 = eps_ctm*lcs.
    """
    r = annexL_residual(fR1k, fR3k, kO, kG)
    f1, f3 = r["fFt1_ef"], r["fFt3_ef"]

    def line(w):
        return max(f1 + (f3 - f1) * (w - 0.5) / 2.0, 0.0)

    w0 = fctm / Ecm * lcs
    if include_matrix:
        nodes_w = [0.0, w0, wu]
        nodes_s = [fctm, line(w0), line(wu)]
    else:
        nodes_w = [0.0, wu]
        nodes_s = [line(0.0), line(wu)]
    return CohesiveLaw(nodes_w, nodes_s, name="EC2 Annex L (w-form)")


def law_envelope(*laws: CohesiveLaw, w_max: float | None = None, n: int = 400, name="envelope") -> CohesiveLaw:
    """Upper envelope sigma(w) = max_i sigma_i(w).

    Implements the MC2010 5.6.5 rule: the plain-concrete branch governs until it intersects the
    residual (fibre) branch.
    """
    if w_max is None:
        w_max = max(l.w[-1] for l in laws)
    nodes = np.unique(np.concatenate([l.w for l in laws] + [np.geomspace(1e-4, w_max, n)]))
    nodes = nodes[nodes <= w_max]
    if nodes[0] != 0.0:
        nodes = np.concatenate([[0.0], nodes])
    s = np.max(np.vstack([l.sigma(nodes) for l in laws]), axis=0)
    # simplify collinear points
    keep = [0]
    for i in range(1, len(nodes) - 1):
        a, b, c = keep[-1], i, i + 1
        pred = s[a] + (s[c] - s[a]) * (nodes[b] - nodes[a]) / (nodes[c] - nodes[a])
        if abs(pred - s[b]) > 1e-6 * max(1.0, s.max()):
            keep.append(i)
    keep.append(len(nodes) - 1)
    return CohesiveLaw(nodes[keep], s[keep], name=name)


# --------------------------------------------------------------------------------------
# Concrete in compression
# --------------------------------------------------------------------------------------

@dataclass
class CompressionLaw:
    """sigma_c(eps) for eps <= 0 (compression negative). Returns negative stresses."""

    kind: str
    fc: float
    Ec: float
    eps_c1: float  # positive number (strain at peak)
    eps_cu: float  # positive number (ultimate strain)
    k_post: float | None = None  # Annex L softening parameter (L.13)
    n: float = 2.0
    eps_c2: float = 0.002

    def sigma(self, eps):
        e = -np.minimum(np.asarray(eps, dtype=float), 0.0)  # positive compressive strain
        if self.kind == "elastic-plastic":
            return -np.minimum(self.Ec * e, self.fc)
        if self.kind == "parabola-rectangle":
            ratio = np.clip(e / self.eps_c2, 0.0, 1.0)
            s = self.fc * (1.0 - (1.0 - ratio) ** self.n)
        else:  # EC2 Formula (5.6) (Sargin-type) with optional Annex L post-peak parameter
            eta = e / self.eps_c1
            k = 1.05 * self.Ec * self.eps_c1 / self.fc
            with np.errstate(divide="ignore", invalid="ignore"):
                s = self.fc * (k * eta - eta ** 2) / (1.0 + (k - 2.0) * eta)
                # the descending branch reaches zero at eta = k; beyond it (and, for k < 2,
                # across the pole at eta = 1/(2-k) >= k) the formula has no meaning -> zero
                s = np.where(eta >= k, 0.0, s)
                if self.k_post is not None:
                    kp = self.k_post
                    s_post = self.fc * (kp * eta - eta ** 2) / (1.0 + (kp - 2.0) * eta)
                    s_post = np.where(eta >= kp, 0.0, s_post)
                    s = np.where(eta > 1.0, s_post, s)
            s = np.clip(np.nan_to_num(s), 0.0, self.fc)
        s = np.where(e > self.eps_cu, 0.0, s)
        return -s


def compression_ec2(fcm: float, Ecm: float | None = None) -> CompressionLaw:
    """EN 1992-1-1:2023 Formulae (5.6)-(5.10) for non-linear analysis (plain concrete)."""
    Ecm = Ecm or ecm_from_fcm(fcm)
    eps_c1 = min(0.7 * fcm ** (1.0 / 3.0), 2.8) / 1000.0
    eps_cu1 = min(2.8 + 14.0 * (1.0 - fcm / 108.0) ** 4, 3.5) / 1000.0
    return CompressionLaw("ec2", fcm, Ecm, eps_c1, eps_cu1)


def compression_annexL(fcm: float, fR1k: float, Ecm: float | None = None) -> CompressionLaw:
    """EN 1992-1-1:2023 Annex L, Formulae (L.12)-(L.14) (Ruiz et al., 2023, Hormigón y Acero)."""
    Ecm = Ecm or ecm_from_fcm(fcm)
    eps_c1 = 0.7 * fcm ** (1.0 / 3.0) * (1.0 + 0.03 * fR1k) / 1000.0
    k = 1.0 + 20.0 / np.sqrt(82.0 - 2.2 * fR1k)
    return CompressionLaw("annexL", fcm, Ecm, eps_c1, k * eps_c1, k_post=k)


def compression_elastic_plastic(fc: float, Ec: float) -> CompressionLaw:
    """Elastic - perfectly plastic compression without crushing limit.

    Used to simulate EN 14651 prisms, where crushing does not govern (Olesen, 2001, adopted a
    linear-elastic compression zone in the hinge).
    """
    return CompressionLaw("elastic-plastic", fc, Ec, fc / Ec, np.inf)


def compression_parabola_rectangle(fcd: float, sfrc: bool = False) -> CompressionLaw:
    """Design parabola-rectangle; Annex L L.8.1(4): eps_c2 = 0.0025, eps_cu = 0.006 for SFRC."""
    if sfrc:
        return CompressionLaw("parabola-rectangle", fcd, 1.0, 0.0025, 0.006, eps_c2=0.0025)
    return CompressionLaw("parabola-rectangle", fcd, 1.0, 0.002, 0.0035, eps_c2=0.002)


# --------------------------------------------------------------------------------------
# Reinforcing steel and crack bridging by bonded bars (Tension Chord Model)
# --------------------------------------------------------------------------------------

@dataclass
class Steel:
    fy: float = 500.0
    fu: float = 540.0
    Es: float = 200000.0
    eps_u: float = 0.05

    @property
    def Esh(self) -> float:
        return (self.fu - self.fy) / max(self.eps_u - self.fy / self.Es, 1e-9)

    def sigma(self, eps):
        eps = np.asarray(eps, dtype=float)
        ey = self.fy / self.Es
        a = np.abs(eps)
        s = np.where(a <= ey, self.Es * a, self.fy + self.Esh * (a - ey))
        s = np.minimum(s, self.fu)
        return np.sign(eps) * s


@dataclass
class Rebar:
    """A layer of bonded bars crossing the crack.

    y   : height of the bar centroid measured from the tension (bottom) face [mm]
    As  : steel area [mm^2]
    phi : bar diameter [mm]
    tau_b0, tau_b1 : bond stresses before/after yielding (TCM: 2 fctm and fctm).
    """

    y: float
    As: float
    phi: float
    steel: Steel = field(default_factory=Steel)
    tau_b0: float = 5.0
    tau_b1: float = 2.5

    def w_of_sigma(self, sig):
        """Crack opening produced by a bar stress sig at the crack (both crack faces)."""
        st = self.steel
        sig = np.asarray(sig, dtype=float)
        se = np.minimum(sig, st.fy)
        w = se ** 2 * self.phi / (4.0 * self.tau_b0 * st.Es)
        dp = np.maximum(sig - st.fy, 0.0)
        w_p = 2.0 * (st.fy / st.Es * dp + dp ** 2 / (2.0 * st.Esh)) * self.phi / (4.0 * self.tau_b1)
        return w + w_p

    def sigma_of_w(self, w):
        st = self.steel
        w = np.maximum(np.asarray(w, dtype=float), 0.0)
        wy = st.fy ** 2 * self.phi / (4.0 * self.tau_b0 * st.Es)
        s_el = np.sqrt(4.0 * self.tau_b0 * st.Es * w / self.phi)
        # post-yield: 2*(ey*dp + dp^2/(2Esh))*phi/(4 tau1) = w - wy -> quadratic in dp
        c = (w - wy) * 4.0 * self.tau_b1 / (2.0 * self.phi)
        A = 1.0 / (2.0 * st.Esh)
        B = st.fy / st.Es
        dp = (-B + np.sqrt(B ** 2 + 4.0 * A * np.maximum(c, 0.0))) / (2.0 * A)
        s_pl = st.fy + dp
        s = np.where(w <= wy, s_el, s_pl)
        return np.minimum(s, st.fu)

    @property
    def w_rupture(self) -> float:
        return float(self.w_of_sigma(self.steel.fu))


def tcm_bond(fctm: float):
    """Bond stresses of the Tension Chord Model (Marti et al., 1998): tau_b0 = 2 fctm, tau_b1 = fctm."""
    return 2.0 * fctm, fctm
