"""Tests of the member checks (1.1): mechanical shear, N-M, columns, punching, code references."""
import math

import numpy as np
import pytest

from frachinge import codes, members
from frachinge.hinge import Section
from frachinge.materials import (Rebar, Steel, compression_annexL, ecm_from_fcm, fctm_from_fck, gf_mc2010,
                                 law_envelope, law_linear_softening, law_mc2010_frc_linear)


def test_neutral_axis_closed_form():
    # x/d = n rho (-1 + sqrt(1 + 2/(n rho)))  -> check equilibrium of the cracked elastic section
    xd = members.masm_neutral_axis(0.01, 200000.0, 30000.0)
    nr = 200000.0 / 30000.0 * 0.01
    assert xd ** 2 / 2 == pytest.approx(nr * (1 - xd), rel=1e-12)


def test_masm_eq26_hand_value():
    # Mari et al. (2020) Eq. (26), mean values; hand calculation
    b, d, a, rho, fcm, st = 250.0, 440.0, 1320.0, 0.0112, 32.1, 1.0
    r = members.shear_masm_sfrc(b, d, a, rho, fcm, st, size_effect=None, orientation="isotropic")
    fct = 0.3 * (fcm - 8.0) ** (2 / 3)
    xi = max(2 / math.sqrt(1 + d / 200) * (d / a) ** 0.2, 0.45)
    s = st / fct
    Vu = (xi * (0.84 - 0.10 * s) * r["x_d"] + 0.08 + 1.10 * s) * fct * b * d
    assert r["V_u_concrete_fibres"] == pytest.approx(Vu, rel=1e-12)


def test_masm_fibres_increase_shear_and_stirrups_add():
    v0 = members.shear_masm_design(300, 450, 1350, 0.012, 30, 0.0)["V_Rd"]
    vf = members.shear_masm_design(300, 450, 1350, 0.012, 30, 0.8)["V_Rd"]
    vs = members.shear_masm_design(300, 450, 1350, 0.012, 30, 0.0, rho_w=0.002)["V_Rd"]
    assert vf > v0 and vs > v0


def test_cccm_size_effect():
    # web term decreases with depth through Gf/d (fracture-mechanics size effect)
    r1 = members.shear_masm_rc(200, 200, 600, 0.01, 35)
    r2 = members.shear_masm_rc(200, 1200, 3600, 0.01, 35)
    assert r2["V_u"] / (200 * 1200) < r1["V_u"] / (200 * 200)


def _column_setup():
    fcm = 38.0
    E = ecm_from_fcm(fcm)
    law = law_envelope(law_linear_softening(fctm_from_fck(30), gf_mc2010(fcm)), law_mc2010_frc_linear(4.0, 3.5))
    sec = Section.rectangle(300, 300)
    comp = compression_annexL(fcm, 4.0)
    rb = [Rebar(40, 452, 12, Steel(500, 540)), Rebar(260, 452, 12, Steel(500, 540))]
    return sec, law, comp, E, rb, fcm


def test_interaction_hinge_shape():
    sec, law, comp, E, rb, fcm = _column_setup()
    Ns = np.linspace(-0.9 * fcm * sec.area, 0.0, 6)
    r = members.interaction_hinge(sec, law, comp, E, rebars=rb, N_values=Ns, n_layers=100, n=40)
    assert np.all(np.isfinite(r["M"]))
    # moderate compression increases the moment capacity (balanced-point bulge)
    assert np.nanmax(r["M"]) > r["M"][-1]


def test_column_slenderness_reduces_capacity():
    sec, law, comp, E, rb, fcm = _column_setup()
    grid = np.linspace(0.05, 1.0, 12) * fcm * sec.area
    n1 = members.column_buckling_hinge(sec, law, comp, E, 2000, 20, rebars=rb, n_layers=80, n=40, N_grid=grid)
    n2 = members.column_buckling_hinge(sec, law, comp, E, 9000, 20, rebars=rb, n_layers=80, n=40, N_grid=grid)
    assert 0 < n2["N_u"] < n1["N_u"]
    assert n2["N_u"] < n2["N_euler"]


def test_punching_fibres_increase_capacity():
    p0 = members.punching_csct_sfrc(150, 200, 1100, 38, 550, 150000)
    p1 = members.punching_csct_sfrc(150, 200, 1100, 38, 550, 150000, fibre_law=law_mc2010_frc_linear(4, 3.5))
    assert p1["V_R"] > p0["V_R"] and p1["psi"] > p0["psi"]
    # failure criterion of Muttoni (2008) satisfied at the solution
    Vc = 0.75 * p0["b0"] * 200 * math.sqrt(38) / (1 + 15 * p0["psi"] * 200 / 32)
    assert p0["V_Rc"] == pytest.approx(Vc, rel=1e-9)


def test_code_references():
    assert codes.VRd_max_strut(300, 400, 30, 1.5, 1.0) == pytest.approx(300 * 400 * 0.6 * 0.88 * 20 / 2)
    assert codes.annexL_stirrup_replacement_ok(0.5, 30) and not codes.annexL_stirrup_replacement_ok(0.4, 30)
    nc = codes.nominal_curvature_moment(1000e3, 0.0, 5000, 260, 300, 9e4, 900, 30, 500)
    assert nc["e_i"] == pytest.approx(12.5) and 0 < nc["Kr"] <= 1
    k = 1 + math.sqrt(200 / 200)
    v = 0.18 / 1.5 * k * (100 * 0.008 * 30) ** (1 / 3)
    assert codes.punching_ec2_2004(1000, 200, 0.008, 30) == pytest.approx(v * 1000 * 200)


def test_bazant_size_factor_on_fibre_term():
    assert members.size_factor_bazant(254.0) == pytest.approx(1.0)
    assert members.size_factor_bazant(100.0) == 1.0
    assert members.size_factor_bazant(1270.0) == pytest.approx(math.sqrt(2 / 6))
    r0 = members.shear_masm_sfrc(250, 1440, 4320, 0.0101, 32.1, 1.65, size_effect=None)
    r1 = members.shear_masm_sfrc(250, 1440, 4320, 0.0101, 32.1, 1.65)
    assert r1["V_fibres"] == pytest.approx(r0["V_fibres"] * members.size_factor_bazant(1440))
    assert r1["V_u"] < r0["V_u"]


def test_orientation_bounds():
    ra = members.shear_masm_sfrc(300, 450, 1350, 0.012, 38, 1.2)
    ri = members.shear_masm_sfrc(300, 450, 1350, 0.012, 38, 1.2, orientation="isotropic")
    assert ri["eta_orientation"] == 1.0
    assert ra["eta_orientation"] == pytest.approx(1 / math.sqrt(1 + ra["cot_theta_cr"] ** 2))
    assert ra["V_u"] < ri["V_u"]


def test_cccm_reproduces_sarkhosh_database_statistics():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
    from data_sarkhosh2014 import DB
    r = np.array([Vu * 1e3 / members.shear_masm_rc(bw, d, a, As / (bw * d), fc)["V_u"]
                  for (_, bw, a, d, As, fc, Vu, *_) in DB])
    assert len(r) == 393
    assert 1.0 < r.mean() < 1.1 and r.std(ddof=1) / r.mean() < 0.16
