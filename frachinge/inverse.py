"""
Inverse identification of cohesive laws from the EN 14651 notched-prism test.

Two procedures are provided:
  * ``identify_fct``: tensile strength of the matrix that reproduces the limit of
    proportionality F_L (with the MC2010 plain-concrete softening, GF from Eq. 5.1-9 or given).
  * ``identify_law``: residual (fibre) law sigma_F(w) defined by nodes at w_j; the nodal stresses
    are fitted (least squares) so that the simulated test reproduces the given fR,j values.
    With nodes at the crack openings corresponding to CMOD_j it is a point-by-point inverse
    analysis in the spirit of Amin, Foster & Muttoni (2015) and the RILEM TC 162-TDF approach.

The composite law is the upper envelope of the plain-concrete law and the fibre law
(MC2010, 5.6.5).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, least_squares

from .materials import CohesiveLaw, law_mc2010_plain, law_envelope
from .hinge import en14651_prism

CMODS = np.array([0.5, 1.5, 2.5, 3.5])


def composite_law(fct: float, GF: float, w_nodes, s_nodes, w_end: float = 4.5) -> CohesiveLaw:
    w_nodes = np.asarray(w_nodes, float)
    s_nodes = np.maximum(np.asarray(s_nodes, float), 0.0)
    # fibre law: linear extrapolation from the first two nodes to w = 0, and from the last
    # two nodes up to w_end (clipped at zero)
    s0 = max(s_nodes[0] - (s_nodes[1] - s_nodes[0]) / (w_nodes[1] - w_nodes[0]) * w_nodes[0], 0.0)
    sl = s_nodes[-1] + (s_nodes[-1] - s_nodes[-2]) / (w_nodes[-1] - w_nodes[-2]) * (w_end - w_nodes[-1])
    wf = np.concatenate([[0.0], w_nodes, [w_end]])
    sf = np.concatenate([[min(s0, fct)], s_nodes, [max(sl, 0.0)]])
    fibre = CohesiveLaw(wf, sf, name="fibre (identified)")
    plain = law_mc2010_plain(fct, GF)
    return law_envelope(plain, fibre, w_max=w_end, name="identified composite law"), fibre


def identify_fct(FL: float, E: float, GF: float, comp, b=150.0, h=150.0, a0=25.0, L=500.0):
    """Matrix tensile strength reproducing the measured F_L (plain-concrete softening)."""
    def g(ft):
        law = law_mc2010_plain(ft, GF)
        return en14651_prism(law, comp, E, b, h, a0, L, cmod_max=0.1, n=60)["FL"] - FL
    return brentq(g, 0.5, 15.0, xtol=1e-4)


def identify_law(fR: dict, fct: float, GF: float, E: float, comp, w_nodes=None,
                 use=(1, 2, 3, 4), b=150.0, h=150.0, a0=25.0, L=500.0, n=90):
    """Least-squares identification of the nodal stresses of the fibre law.

    fR : dict {1: fR1, 2: fR2, 3: fR3, 4: fR4} (MPa). ``use`` selects which fR,j are fitted.
    w_nodes : crack openings of the nodes (default: one node per fitted fR,j at w ~ 0.85 CMOD_j).
    """
    use = tuple(use)
    if w_nodes is None:
        w_nodes = 0.85 * CMODS[[j - 1 for j in use]]
    w_nodes = np.asarray(w_nodes, float)
    target = np.array([fR[j] for j in use])

    def sim(p):
        law, _ = composite_law(fct, GF, w_nodes, p)
        r = en14651_prism(law, comp, E, b, h, a0, L, cmod_max=3.6, n=n)
        return r

    def res(p):
        r = sim(p)
        return np.array([r["fR"][j] for j in use]) - target

    p0 = 0.37 * np.interp(w_nodes, 0.85 * CMODS, [fR.get(1, target[0]), fR.get(2, target[0]),
                                                   fR.get(3, target[-1]), fR.get(4, target[-1])])
    sol = least_squares(res, p0, bounds=(0.0, 3.0 * max(target)), xtol=1e-8, ftol=1e-8)
    law, fibre = composite_law(fct, GF, w_nodes, sol.x)
    r = en14651_prism(law, comp, E, b, h, a0, L, cmod_max=3.6, n=n)
    return {"law": law, "fibre": fibre, "nodes_w": w_nodes, "nodes_s": sol.x, "sim": r,
            "residual": sol.fun, "success": sol.success}
