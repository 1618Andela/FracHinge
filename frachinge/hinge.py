"""
Cohesive cracked-hinge model for plain, fibre-reinforced and hybrid (bars + fibres) sections.

The formulation follows the non-linear hinge of Olesen (2001, J. Eng. Mech. 127(3):272-280),
adopted by RILEM TC 162-TDF (2002), in which a band of width ``s`` around the crack is divided
into independent layers. Each layer is either elastic, in compression, or cracked; in a cracked
layer the band elongation u is shared between the elastic deformation of the band and the
crack opening w,  u = s*sigma/E + w,  with sigma = sigma_w(w)  (Hillerborg's fictitious crack).

Extensions implemented here (documented in docs/manual):
  * non-linear concrete in compression (EC2 Formula (5.6) or Annex L (L.12)-(L.14));
  * bonded rebars crossing the crack, whose force depends on the crack opening at the bar level
    through the rigid-plastic bond law of the Tension Chord Model (Marti et al., 1998);
  * arbitrary axial force N (compression negative);
  * variable width b(y) (rectangular, T and I sections as stacked rectangles).

Sign conventions: y is measured from the tension (bottom) face; theta is the total angular
deformation of the hinge (theta = 2*phi in Olesen's notation); tensile stresses positive.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import brentq

from .materials import CohesiveLaw, CompressionLaw, Rebar


@dataclass
class Section:
    """Section made of stacked rectangles: list of (y_bottom, y_top, width)."""

    parts: list
    name: str = "section"

    @staticmethod
    def rectangle(b: float, h: float) -> "Section":
        return Section([(0.0, h, b)], name=f"rect {b:g}x{h:g}")

    @staticmethod
    def tee(bw: float, h: float, bf: float, hf: float) -> "Section":
        return Section([(0.0, h - hf, bw), (h - hf, h, bf)], name=f"T bw={bw:g} h={h:g} bf={bf:g} hf={hf:g}")

    @property
    def h(self) -> float:
        return max(p[1] for p in self.parts)

    def width(self, y):
        y = np.asarray(y, dtype=float)
        b = np.zeros_like(y)
        for y0, y1, w in self.parts:
            b = np.where((y >= y0) & (y <= y1), w, b)
        return b

    @property
    def area(self) -> float:
        return sum((p[1] - p[0]) * p[2] for p in self.parts)

    @property
    def centroid(self) -> float:
        return sum((p[1] - p[0]) * p[2] * 0.5 * (p[0] + p[1]) for p in self.parts) / self.area

    @property
    def inertia(self) -> float:
        yc = self.centroid
        I = 0.0
        for y0, y1, w in self.parts:
            hh = y1 - y0
            I += w * hh ** 3 / 12.0 + w * hh * (0.5 * (y0 + y1) - yc) ** 2
        return I


@dataclass
class HingeState:
    theta: float
    y0: float
    M: float
    N: float
    w_bottom: float
    crack_tip: float
    eps_top: float
    sigma_s: np.ndarray
    w_s: np.ndarray
    y: np.ndarray
    sigma: np.ndarray
    w: np.ndarray
    u: np.ndarray


@dataclass
class CrackedHinge:
    section: Section
    law: CohesiveLaw
    comp: CompressionLaw
    E: float
    s: float | None = None
    rebars: list = field(default_factory=list)
    n_layers: int = 400
    w_limit: float | None = None  # optional crack-opening limit (e.g. wu)
    c_fdur: float = 0.0  # Annex L L.6(3)-(4): residual strength disregarded for y < c_fdur

    def __post_init__(self):
        h = self.section.h
        if self.s is None:
            self.s = 0.5 * h  # Olesen (2001), RILEM TC 162-TDF (2002)
        dy = h / self.n_layers
        self.y = (np.arange(self.n_layers) + 0.5) * dy
        self.dA = self.section.width(self.y) * dy
        w0, w1, a, b = self.law.segments()
        self._seg = (w0, w1, a, b)
        self.yref = 0.5 * h

    # ------------------------------------------------------------------ layer response
    def _layer_stress(self, u):
        """Stresses and crack openings of all layers for band elongations u (vector)."""
        s, E, law = self.s, self.E, self.law
        sig = np.empty_like(u)
        w = np.zeros_like(u)
        comp = u <= 0.0
        sig[comp] = self.comp.sigma(u[comp] / s)
        u_cr = s * law.ft / E
        el = (~comp) & (u <= u_cr)
        sig[el] = E * u[el] / s
        cr = (~comp) & (u > u_cr)
        if np.any(cr):
            uc = u[cr][:, None]
            w0, w1, a, b = self._seg
            den = 1.0 + s * b / E
            with np.errstate(divide="ignore", invalid="ignore"):
                ww = (uc - s * a / E) / den
            valid = (den > 0) & (ww >= w0 - 1e-12) & (ww <= w1 + 1e-12)
            ww = np.where(valid, ww, np.inf)
            wmin = ww.min(axis=1)
            # beyond the last node: sigma = 0 (or constant tail)
            if law.tail == "zero":
                w_tail = uc[:, 0]
                tail_ok = w_tail >= law.w[-1]
            else:
                w_tail = uc[:, 0] - s * law.s[-1] / E
                tail_ok = w_tail >= law.w[-1]
            wsel = np.where(np.isfinite(wmin), wmin, np.where(tail_ok, w_tail, uc[:, 0]))
            w[cr] = wsel
            sig[cr] = law.sigma(wsel)
            if self.c_fdur > 0:
                kill = cr & (self.y < self.c_fdur)
                sig[kill] = 0.0
                w[kill] = u[kill]
        return sig, w

    def _resultants(self, theta, y0):
        u = theta * (y0 - self.y)
        sig, w = self._layer_stress(u)
        F = sig * self.dA
        N = F.sum()
        M = np.sum(F * (self.yref - self.y))
        ss, ws = [], []
        for rb in self.rebars:
            u_s = theta * (y0 - rb.y)
            if u_s <= 0:
                sig_s = float(rb.steel.sigma(u_s / self.s))
                w_s = 0.0
            else:
                # crack opening at bar level: interpolate the concrete layer solution
                w_s = float(np.interp(rb.y, self.y, w))
                if w_s <= 0.0:
                    sig_s = float(rb.steel.sigma(u_s / self.s))
                else:
                    sig_s = float(rb.sigma_of_w(w_s))
            N += rb.As * sig_s
            M += rb.As * sig_s * (self.yref - rb.y)
            ss.append(sig_s)
            ws.append(w_s)
        return N, M, u, sig, w, np.array(ss), np.array(ws)

    def solve(self, theta: float, N_ext: float = 0.0, y0_guess: float | None = None) -> HingeState:
        h = self.section.h
        f = lambda y0: self._resultants(theta, y0)[0] - N_ext
        ref = y0_guess if y0_guess is not None else 0.5 * h
        y0 = None
        # 1) local bracketing around the previous solution (continuation)
        fr = f(ref)
        if fr == 0.0:
            y0 = ref
        else:
            step = 0.01 * h
            for _ in range(12):
                lo, hi = ref - step, ref + step
                flo, fhi = f(lo), f(hi)
                if flo * fr < 0:
                    y0 = brentq(f, lo, ref, xtol=1e-9 * h)
                    break
                if fhi * fr < 0:
                    y0 = brentq(f, ref, hi, xtol=1e-9 * h)
                    break
                step *= 1.8
                if step > 2.0 * h:
                    break
        # 2) global scan (fallback)
        if y0 is None:
            grid = np.linspace(-0.5 * h, 1.5 * h, 161)
            vals = np.array([f(g) for g in grid])
            roots = []
            for i in range(len(grid) - 1):
                if vals[i] == 0:
                    roots.append(grid[i])
                elif vals[i] * vals[i + 1] < 0:
                    roots.append(brentq(f, grid[i], grid[i + 1], xtol=1e-9 * h))
            if not roots:
                raise RuntimeError(f"No equilibrium found for theta={theta:g}")
            y0 = min(roots, key=lambda r: abs(r - ref))
        N, M, u, sig, w, ss, ws = self._resultants(theta, y0)
        cracked = w > 0
        tip = float(self.y[cracked].max()) if np.any(cracked) else 0.0
        return HingeState(theta, y0, M, N, float(w[0]), tip, theta * (y0 - h) / self.s,
                          ss, ws, self.y, sig, w, u)

    def response(self, theta_max: float | None = None, n: int = 120, N_ext: float = 0.0,
                 stop_on_failure: bool = True):
        """Moment-rotation response under monotonically increasing rotation.

        Returns a dict of arrays and a ``failure`` string (crushing, rupture, w_limit, or None).
        """
        h = self.section.h
        th_cr = 2.0 * self.law.ft / self.E * self.s / (0.5 * h) if theta_max is None else None
        if theta_max is None:
            theta_max = 0.04 * self.s / h * 50
        theta0 = max(1e-3 * self.law.ft / self.E * self.s / h, 1e-9)
        # rotation at first cracking of an elastic section (ft at the soffit): theta = 2 ft s/(E h)
        th_c = 2.0 * self.law.ft * self.s / (self.E * h)
        dense = np.linspace(0.6 * th_c, 4.0 * th_c, max(n // 2, 40))
        thetas = np.unique(np.concatenate([np.geomspace(theta0, theta_max, n), dense[dense < theta_max]]))
        thetas = np.concatenate([[0.0], thetas])
        out = {k: [] for k in ["theta", "M", "y0", "w_bottom", "crack_tip", "eps_top", "sigma_s", "w_s"]}
        failure = None
        y0 = 0.5 * h
        for th in thetas:
            if th == 0.0:
                for k in out:
                    out[k].append(0.0 if k not in ("sigma_s", "w_s") else np.zeros(len(self.rebars)))
                out["y0"][-1] = self.section.centroid
                continue
            st = self.solve(th, N_ext, y0)
            y0 = st.y0
            out["theta"].append(th)
            out["M"].append(st.M)
            out["y0"].append(st.y0)
            out["w_bottom"].append(st.w_bottom)
            out["crack_tip"].append(st.crack_tip)
            out["eps_top"].append(st.eps_top)
            out["sigma_s"].append(st.sigma_s)
            out["w_s"].append(st.w_s)
            if stop_on_failure:
                if -st.eps_top > self.comp.eps_cu:
                    failure = "concrete crushing"
                elif any(ss >= rb.steel.fu - 1e-6 for ss, rb in zip(st.sigma_s, self.rebars)):
                    failure = "rebar rupture"
                elif self.w_limit is not None and st.w_bottom > self.w_limit:
                    failure = "crack opening limit"
                if failure:
                    break
        res = {k: np.array(v) for k, v in out.items()}
        res["failure"] = failure
        return res


# ------------------------------------------------------------------------------------------
# Beams built around a midspan cracked hinge
# ------------------------------------------------------------------------------------------

@dataclass
class SimplyBeam:
    """Simply supported beam with the cracked hinge at midspan.

    load = "3pb" (central point load) or "4pb" (two loads at L/3).
    a0 : notch depth (the hinge section is then the ligament h - a0); I uses the full depth.
    """

    hinge: CrackedHinge
    L: float
    load: str = "3pb"
    a0: float = 0.0
    h_full: float | None = None
    b_full: float | None = None

    def response(self, **kw):
        r = self.hinge.response(**kw)
        E = self.hinge.E
        if self.h_full is not None:
            I = (self.b_full or self.hinge.section.width(0.0)) * self.h_full ** 3 / 12.0
        else:
            I = self.hinge.section.inertia
        M = r["M"]
        if self.load == "3pb":
            P = 4.0 * M / self.L
            d_el = P * self.L ** 3 / (48.0 * E * I)
        else:
            P = 6.0 * M / self.L
            d_el = 23.0 * P * self.L ** 3 / (1296.0 * E * I)
        th_el = M * self.hinge.s / (E * self.hinge.section.inertia)
        d = d_el + np.maximum(r["theta"] - th_el, 0.0) * self.L / 4.0
        # crack mouth opening displacement (opening at the beam soffit)
        cmod = []
        for th, y0 in zip(r["theta"], r["y0"]):
            u_tip = th * (y0 - self.hinge.y[0])
            cmod.append(max(u_tip, 0.0) + th * self.a0)
        r.update({"P": P, "delta": d, "CMOD": np.array(cmod)})
        return r


def en14651_prism(law: CohesiveLaw, comp: CompressionLaw, E: float, b: float = 150.0,
                  h: float = 150.0, a0: float = 25.0, L: float = 500.0, n_layers: int = 300,
                  cmod_max: float = 4.0, n: int = 160):
    """Simulates the EN 14651 notched-prism test and returns F-CMOD and fR,j (j = 1..4)."""
    hsp = h - a0
    hinge = CrackedHinge(Section.rectangle(b, hsp), law, comp, E, s=0.5 * hsp, n_layers=n_layers)
    beam = SimplyBeam(hinge, L, "3pb", a0=a0, h_full=h, b_full=b)
    theta_max = 2.0 * cmod_max / hsp
    r = beam.response(theta_max=theta_max, n=n, stop_on_failure=False)
    F, cmod = r["P"], r["CMOD"]
    fR = {}
    for j, c in zip((1, 2, 3, 4), (0.5, 1.5, 2.5, 3.5)):
        Fj = float(np.interp(c, cmod, F))
        fR[j] = 3.0 * Fj * L / (2.0 * b * hsp ** 2)
    mask = cmod <= 0.05
    FL = float(F[mask].max()) if np.any(mask) else float(F.max())
    fL = 3.0 * FL * L / (2.0 * b * hsp ** 2)
    return {"F": F, "CMOD": cmod, "delta": r["delta"], "fR": fR, "FL": FL, "fL": fL, "raw": r}
