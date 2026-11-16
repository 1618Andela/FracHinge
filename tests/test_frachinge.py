"""Unit and verification tests for FracHinge (run with: pytest -q)."""
import numpy as np
import pytest

from frachinge.materials import (CohesiveLaw, Rebar, Steel, compression_annexL, compression_ec2,
                                 compression_elastic_plastic, ecm_from_fcm, fctm_from_fck,
                                 gf_mc2010, law_annexL, law_envelope, law_mc2010_frc_linear,
                                 law_mc2010_plain, law_mc2010_rigid_plastic, mc2010_residual,
                                 annexL_residual)
from frachinge.hinge import CrackedHinge, Section, en14651_prism
from frachinge import codes, fracture
from frachinge.safety import ecov
from frachinge.design import DesignInput, run_design

FCK = 30.0
FCM = 38.0
E = ecm_from_fcm(FCM)
FT = fctm_from_fck(FCK)
GF = gf_mc2010(FCM)


def test_basic_properties():
    assert FT == pytest.approx(0.3 * 30 ** (2 / 3))
    assert GF == pytest.approx(73 * 38 ** 0.18 / 1000)
    assert E == pytest.approx(22000 * 3.8 ** 0.3)


def test_mc2010_plain_law_energy():
    law = law_mc2010_plain(FT, GF)
    # area of the MC2010 bilinear law = GF (by construction of w1 and wc)
    assert law.energy() == pytest.approx(GF, rel=1e-9)
    assert law.sigma(GF / FT) == pytest.approx(0.2 * FT)


def test_mc2010_linear_frc():
    fFts, fFtu = mc2010_residual(4.0, 3.6)
    assert fFts == pytest.approx(1.8)
    assert fFtu == pytest.approx(0.5 * 3.6 - 0.2 * 4.0)  # wu = CMOD3 = 2.5 mm


def test_annexL_values():
    r = annexL_residual(4.0, 3.6)
    assert r["fFtu_ef"] == pytest.approx(0.33 * 3.6)
    assert r["fFt3_ef"] == pytest.approx(0.57 * 3.6 - 0.26 * 4.0)
    law = law_annexL(4.0, 3.6, FT, E)
    assert law.sigma(0.5) == pytest.approx(0.37 * 4.0, rel=1e-6)
    assert law.sigma(2.5) == pytest.approx(0.57 * 3.6 - 0.26 * 4.0, rel=1e-6)


def test_annexL_compression():
    c = compression_annexL(FCM, 4.0)
    assert c.eps_c1 == pytest.approx(0.7 * FCM ** (1 / 3) * (1 + 0.12) / 1000)
    assert c.k_post == pytest.approx(1 + 20 / np.sqrt(82 - 8.8))
    # continuity at the peak
    assert -c.sigma(-c.eps_c1) == pytest.approx(FCM, rel=1e-9)


def test_envelope_is_max():
    a = law_mc2010_plain(FT, GF)
    b = law_mc2010_frc_linear(4.0, 3.6)
    env = law_envelope(a, b)
    w = np.linspace(0, 2.5, 300)
    assert np.allclose(env.sigma(w), np.maximum(a.sigma(w), b.sigma(w)), atol=2e-3)


def test_tcm_rebar_inverse():
    rb = Rebar(40, 400, 16, Steel(500, 540, 200000, 0.05), 2 * FT, FT)
    for s in [50, 200, 499, 510, 535]:
        assert rb.sigma_of_w(rb.w_of_sigma(s)) == pytest.approx(s, rel=1e-8)


def test_elastic_limit_moment():
    """Before cracking, M = E I theta / s (elastic band of width s)."""
    law = law_mc2010_plain(FT, GF)
    sec = Section.rectangle(100, 200)
    hg = CrackedHinge(sec, law, compression_elastic_plastic(FCM, E), E, n_layers=800)
    th = 0.3 * 2 * FT * hg.s / (E * 200)
    st = hg.solve(th)
    assert st.M == pytest.approx(E * sec.inertia * th / hg.s, rel=2e-3)


def test_rigid_plastic_limit():
    """A perfectly plastic cohesive law gives M -> f b h^2 / 2 at large rotation."""
    f = 1.2
    law = CohesiveLaw([0.0, 50.0], [f, f], tail="constant")
    hg = CrackedHinge(Section.rectangle(100, 200), law, compression_elastic_plastic(1e4, E), E,
                      n_layers=800)
    st = hg.solve(0.2)
    assert st.M == pytest.approx(f * 100 * 200 ** 2 / 2, rel=0.01)


def test_size_effect_vs_mc2010():
    law = law_mc2010_plain(FT, GF)
    for h, tol in [(100, 0.06), (400, 0.08), (800, 0.03)]:
        hg = CrackedHinge(Section.rectangle(100, h), law, compression_ec2(FCM, E), E, n_layers=400)
        r = hg.response(theta_max=40 * FT / E, n=120)
        ratio = r["M"].max() / (100 * h ** 2 / 6) / FT
        assert ratio == pytest.approx(codes.fct_flexural_mc2010(FT, h) / FT, rel=tol)


def test_en14651_reproduces_rigid_plastic():
    """With a constant residual stress f, fR3 approaches 3 f (MC2010 Eq. 5.6-4 rationale)."""
    f = 1.0
    law = law_envelope(law_mc2010_plain(FT, GF), CohesiveLaw([0, 10.0], [f, f]))
    r = en14651_prism(law, compression_elastic_plastic(FCM, E), E, cmod_max=3.0, n=80)
    assert r["fR"][3] == pytest.approx(3.0 * f, rel=0.08)


def test_classification():
    assert codes.classify_mc2010(4.0, 3.6) == "4c"
    assert codes.classify_ec2(4.0, 3.6) == "4 c"
    assert codes.classify_ec2(3.6, 2.6) == "3.5 b"


def test_kappa_G_and_shear():
    assert codes.kappa_G(0.1) == pytest.approx(1.05)
    assert codes.kappa_G(5.0) == pytest.approx(1.5)
    t, eta = codes.tau_RdcF_annexL(0.01, 30, 400, 0.8, 1.2)
    assert eta == pytest.approx(0.6)
    assert t > codes.tau_Rdc(0.01, 30, 400) * 0.6


def test_uls_annexL_rc_matches_hand_calc():
    # RC rectangle, plain concrete: compare with rectangular-block hand calculation (approx.)
    b, h, d, As, fyd = 300, 500, 460, 603, 435
    r = codes.uls_bending_annexL(b, h, 30, 0.0, [(d, As, fyd)], 1.5)
    fcd = codes.eta_cc(30) * 30 / 1.5
    x = As * fyd / (0.8 * fcd * b)
    M_hand = As * fyd * (d - 0.4 * x)
    assert r["M_Rd"] == pytest.approx(M_hand, rel=0.02)


def test_ecov():
    e = ecov(100.0, 100.0 * np.exp(-1.65 * 0.1))
    assert e["V_R"] == pytest.approx(0.1)
    assert e["gamma_R"] == pytest.approx(np.exp(0.8 * 3.8 * 0.1))


def test_senb_factor():
    # Srawley (1976), S/W = 4, load form K = P S /(B W^1.5) f(alpha), f = 1.5 sqrt(alpha) F(alpha);
    # f(0.5) = 2.66 (standard tabulated value)
    a = 0.5
    assert 1.5 * np.sqrt(a) * fracture.senb_geometry_factor(a) == pytest.approx(2.66, abs=0.01)


def test_design_workflow_runs():
    out = run_design(DesignInput(b=1000, h=200, fR1k=3.0, fR3k=2.7, bars=[(40, 393, 10)],
                                 M_Ed=20e6, M_ser=35e6), n=80)
    assert out["M_Rd_hinge"] > 0
    assert out["ULS_ok"]
    assert 0 < out["ecov"]["gamma_R"] < 3


def test_cfdur_reduces_capacity():
    law = law_envelope(law_mc2010_plain(FT, GF), law_mc2010_frc_linear(4.0, 3.6, w_end=3.0))
    comp = compression_elastic_plastic(FCM, E)
    a = CrackedHinge(Section.rectangle(1000, 200), law, comp, E, n_layers=400).solve(0.02).M
    b = CrackedHinge(Section.rectangle(1000, 200), law, comp, E, n_layers=400, c_fdur=20.0).solve(0.02).M
    assert b < a


def test_ec2_compression_no_pole_when_k_below_2():
    # low Ec (Carpinteri et al. 2015: Ec = 22199 MPa, fc = 50.95 MPa) gives k < 2 in EC2 (5.6)
    from frachinge.materials import compression_ec2
    c = compression_ec2(50.95, 22199.0)
    eps = -np.linspace(0.0, c.eps_cu, 2001)
    s = -c.sigma(eps)
    assert np.all(np.isfinite(s)) and s.max() <= c.fc + 1e-9 and s.min() >= 0.0


def test_carpinteri_beam_2phi8_10kg():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
    import data_carpinteri2015 as C
    from frachinge.materials import law_linear_softening, compression_ec2, Rebar, Steel, tcm_bond, fctm_from_fck
    m = C.MAT[10]
    t0, t1 = tcm_bond(fctm_from_fck(m["fc"] - 8.0))
    rb = [Rebar(C.C, 2 * np.pi * 16.0, 8.0, Steel(0.999 * 630.0, 630.0, C.ES, 0.10), t0, t1)]
    comp = compression_ec2(m["fc"], C.EC)
    hg = CrackedHinge(Section.rectangle(C.B, C.H), law_linear_softening(m["ft3"], m["GF3"]), comp, C.EC,
                      rebars=rb, n_layers=300)
    r = hg.response(theta_max=0.06, n=150, stop_on_failure=False)
    ok = -r["eps_top"] <= comp.eps_cu
    P = 4.0 * (r["M"][ok].max() - 24e-6 * C.B * C.H * C.L ** 2 / 8.0) / C.L / 1000.0
    assert abs(P / C.PEXP[10]["2Ø8"][0] - 1.0) < 0.10


def test_cuenca_database_mc2010_mean_ratio():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
    from data_cuenca2018 import DB
    from frachinge import codes
    from frachinge.materials import fctm_from_fck
    r = []
    for (_, _, _, _, bw, d, _, rho, fc, fR1, fR3, fm, Vu) in DB:
        if fm in ("S", "SY"):
            fFts = 0.45 * fR1
            fFtu = max(fFts - 0.6 * (fFts - 0.5 * fR3 + 0.2 * fR1), 0.0)
            r.append(Vu * 1e3 / codes.VRdF_mc2010(bw, d, rho / 100, fc, fFtu, fctm_from_fck(max(fc - 8, 12)), gamma_c=1.0))
    assert len(r) == 43 and 0.95 < np.mean(r) < 1.20


def test_tiberti_database_integrity():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
    import data_tiberti2017 as T
    assert len(T.DB) == 81 and sum(r["N"] for r in T.DB) == 528
    for r in T.DB:
        assert r["fk"][0] <= r["fL"] and all(r["fk"][j + 1] <= r["fR"][j] for j in range(4))
