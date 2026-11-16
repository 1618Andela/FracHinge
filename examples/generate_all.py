"""
Generates every figure (ES/EN), LaTeX table and numerical macro used in the documents.

    python examples/generate_all.py            # writes docs/figures, docs/tables

All numbers quoted in the manual, article, slides and course come from this script
(docs/tables/numbers.tex), so text and code cannot diverge.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import style  # noqa: E402
from data_gao2025 import GAO2025  # noqa: E402
from data_bosco1990 import GRADES as BOSCO_GRADES, BEAMS as BOSCO_BEAMS  # noqa: E402
import data_carpinteri2015 as CARP  # noqa: E402
import data_minelli2014 as MINELLI  # noqa: E402
from data_cuenca2018 import DB as CUENCA_DB  # noqa: E402
import data_tiberti2017 as TIB  # noqa: E402
import data_delarosa2025 as DLR  # noqa: E402
from frachinge.materials import (  # noqa: E402
    CohesiveLaw, Rebar, Steel, characteristic_length, compression_annexL, compression_ec2,
    compression_parabola_rectangle, compression_elastic_plastic, ecm_from_fcm, fctm_from_fck, gf_mc2010, law_annexL,
    law_envelope, law_linear_softening, law_mc2010_frc_linear, law_mc2010_plain, law_mc2010_rigid_plastic,
    mc2010_residual, annexL_residual, tcm_bond)
from frachinge.hinge import CrackedHinge, Section, SimplyBeam, en14651_prism  # noqa: E402
from frachinge.inverse import identify_fct, identify_law  # noqa: E402
from frachinge import codes, fracture, members  # noqa: E402
from frachinge.safety import ecov  # noqa: E402

FIG = os.path.join(ROOT, "docs", "figures")
TAB = os.path.join(ROOT, "docs", "tables")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)
style.apply()
P, LS = style.PALETTE, style.LINESTYLES

NUM: dict[str, str] = {}


def put(key, value, fmt="{:.2f}"):
    NUM[key] = fmt.format(value) if not isinstance(value, str) else value
    return value


def T(lang, en, es):
    return en if lang == "en" else es


def save(fig, name, lang):
    fig.savefig(os.path.join(FIG, f"{name}_{lang}.pdf"))
    plt.close(fig)


def fmt_dec(x, nd=2, lang="en"):
    s = f"{x:.{nd}f}"
    return s


# ======================================================================================
# Reference material used in the worked examples (Example 1)
# ======================================================================================
FCK = 30.0
FCM = FCK + 8.0
ECM = ecm_from_fcm(FCM)
FCTM = fctm_from_fck(FCK)
FCTK = 0.7 * FCTM
GF = gf_mc2010(FCM)
LCH = characteristic_length(ECM, GF, FCTM)
FR1K, FR3K = 4.0, 3.6


def example1():
    put("fck", FCK, "{:.0f}")
    put("fcm", FCM, "{:.0f}")
    put("Ecm", ECM, "{:.0f}")
    put("fctm", FCTM)
    put("fctk", FCTK)
    put("GF", GF * 1000, "{:.1f}")  # N/m
    put("lch", LCH, "{:.0f}")
    put("fRonek", FR1K, "{:.1f}")
    put("fRthreek", FR3K, "{:.1f}")
    put("classMC", codes.classify_mc2010(FR1K, FR3K))
    put("classEC", codes.classify_ec2(FR1K, FR3K))
    fFts, fFtu = mc2010_residual(FR1K, FR3K)
    put("MCfFts", fFts)
    put("MCfFtu", fFtu)
    put("MCfFtuRP", FR3K / 3)
    r = annexL_residual(FR1K, FR3K)
    put("ALfFts", r["fFts_ef"])
    put("ALfFtu", r["fFtu_ef"])
    put("ALfFtone", r["fFt1_ef"])
    put("ALfFtthree", r["fFt3_ef"])
    put("ALfFtud", r["fFtu_ef"] / 1.5)
    put("ALfFtsd", r["fFts_ef"] / 1.5)
    put("ductLone", FR1K / (0.7 * FCTM))  # (L.1) uses fctk,0.05 = 0.7 fctm
    comp = compression_annexL(FCM, FR1K, ECM)
    put("ALepscone", comp.eps_c1 * 1000)
    put("ALk", comp.k_post)
    put("ALepscuone", comp.eps_cu * 1000)
    pl = compression_ec2(FCM, ECM)
    put("ECepscone", pl.eps_c1 * 1000)
    put("ECepscuone", pl.eps_cu * 1000)

    # table of classes (Annex L Table L.2 analytical formulae check)
    rows = []
    for N in codes.STRENGTH_SERIES_EC2:
        rows.append([f"{N:.1f}"] + [f"{f * N:.1f}" for f in (0.5, 0.7, 0.9, 1.1, 1.3)])
    with open(os.path.join(TAB, "tab_classes.tex"), "w") as f:
        f.write("\\begin{tabular}{c ccccc}\n\\toprule\n$N$ [MPa] & a & b & c & d & e\\\\\n\\midrule\n")
        for r_ in rows:
            f.write(" & ".join(r_) + "\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")


def fig_laws():
    plain = law_mc2010_plain(FCTM, GF)
    mcl = law_mc2010_frc_linear(FR1K, FR3K, w_end=3.5)
    rp = law_mc2010_rigid_plastic(FR3K)
    al = law_annexL(FR1K, FR3K, FCTM, ECM, lcs=125.0, wu=2.5)
    env = law_envelope(plain, mcl)
    w = np.linspace(0, 3.5, 700)
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
        a = ax[0]
        a.plot(w, plain.sigma(w), color=P[0], ls=LS[0], label=T(lang, "plain concrete (MC2010)", "hormigón en masa (MC2010)"))
        a.plot(w, mcl.sigma(w), color=P[1], ls=LS[1], label=T(lang, "MC2010 linear (FRC)", "MC2010 lineal (HRF)"))
        a.plot(w, rp.sigma(w), color=P[2], ls=LS[2], label=T(lang, "MC2010 rigid-plastic", "MC2010 rígido-plástico"))
        a.plot(w, al.sigma(w), color=P[3], ls=LS[3], label=T(lang, "EC2 Annex L (w-form)", "EC2 Anejo L (en w)"))
        a.plot(w, env.sigma(w), color=P[6], ls=LS[4], lw=1.2, label=T(lang, "envelope (MC2010 5.6.5)", "envolvente (MC2010 5.6.5)"))
        a.set_xlabel("$w$ [mm]")
        a.set_ylabel(r"$\sigma$ [MPa]")
        a.set_ylim(0, FCTM * 1.1)
        a.legend(loc="upper right")
        a.set_title(T(lang, "(a) cohesive laws, class 4c", "(a) leyes cohesivas, clase 4c"))
        b = ax[1]
        wz = np.linspace(0, 0.3, 400)
        b.plot(wz, plain.sigma(wz), color=P[0], ls=LS[0])
        b.plot(wz, al.sigma(wz), color=P[3], ls=LS[3])
        b.plot(wz, env.sigma(wz), color=P[6], ls=LS[4], lw=1.2)
        b.set_xlabel("$w$ [mm]")
        b.set_title(T(lang, "(b) zoom: matrix-to-fibre transfer", "(b) detalle: transferencia matriz-fibra"))
        b.set_ylim(0, FCTM * 1.1)
        fig.tight_layout()
        save(fig, "fig_laws", lang)
    put("plainwone", GF / FCTM, "{:.3f}")
    put("plainwc", 5 * GF / FCTM, "{:.3f}")


def fig_compression():
    pl = compression_ec2(FCM, ECM)
    e = np.linspace(0, 0.012, 600)
    put("ECkSargin", 1.05 * ECM * pl.eps_c1 / FCM)
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(4.6, 3.2))
        ax.plot(e * 1000, -pl.sigma(-e), color=P[0], ls=LS[0], label=T(lang, "plain (EC2 5.6)", "en masa (EC2 5.6)"))
        for i, fr in enumerate([2.0, 4.0, 8.0]):
            c = compression_annexL(FCM, fr, ECM)
            ax.plot(e * 1000, -c.sigma(-e), color=P[i + 1], ls=LS[i + 1], label=f"SFRC $f_{{R,1k}}$={fr:g} MPa")
        ax.set_xlabel(r"$\varepsilon_c$ [‰]")
        ax.set_ylabel(r"$\sigma_c$ [MPa]")
        ax.legend()
        fig.tight_layout()
        save(fig, "fig_compression", lang)


# ======================================================================================
# Example 2 - EN 14651 forward simulation with code laws ("consistency check")
# ======================================================================================

def example2():
    comp = compression_elastic_plastic(FCM, ECM)
    plain = law_mc2010_plain(FCTM, GF)
    lawA = law_envelope(plain, law_mc2010_frc_linear(FR1K, FR3K, w_end=4.5))
    lawB = law_annexL(FR1K, FR3K, FCTM, ECM, lcs=125.0, wu=4.5)
    lawC = law_envelope(plain, law_mc2010_rigid_plastic(FR3K, wu=4.5))
    rs = {}
    for k, law in zip("ABC", (lawA, lawB, lawC)):
        rs[k] = en14651_prism(law, comp, ECM, cmod_max=4.0, n=140)
    for k in "ABC":
        for j in (1, 2, 3, 4):
            put(f"extwo{k}fR{'abcd'[j - 1]}", rs[k]["fR"][j])
        put(f"extwo{k}fL", rs[k]["fL"])
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        names = [T(lang, "MC2010 linear + plain", "MC2010 lineal + matriz"),
                 T(lang, "EC2 Annex L", "EC2 Anejo L"),
                 T(lang, "MC2010 rigid-plastic + plain", "MC2010 rígido-plástico + matriz")]
        for i, k in enumerate("ABC"):
            ax.plot(rs[k]["CMOD"], rs[k]["F"] / 1000, color=P[i], ls=LS[i], label=names[i])
        hsp = 125.0
        Fj = [f * 2 * 150 * hsp ** 2 / (3 * 500) / 1000 for f in (FR1K, FR3K)]
        ax.plot([0.5, 2.5], Fj, ls="none", marker="s", color="#0b0b0b",
                label=T(lang, "input $F_{R,1}$, $F_{R,3}$", "datos $F_{R,1}$, $F_{R,3}$"))
        ax.set_xlabel("CMOD [mm]")
        ax.set_ylabel("$F$ [kN]")
        ax.set_xlim(0, 4)
        ax.legend()
        fig.tight_layout()
        save(fig, "fig_en14651_forward", lang)


# ======================================================================================
# Validation 1 - Gao et al. (2025) 14 SFRC mixes
# ======================================================================================

def validation_gao():
    rows = []
    ident = {}
    for (id_, fcu, f1, f2, f3, f4, F1, EG) in GAO2025:
        E = EG * 1000.0
        fcm = 0.8 * fcu
        GFm = gf_mc2010(fcm)
        comp = compression_elastic_plastic(fcm, E)
        fct = identify_fct(F1 * 1000.0, E, GFm, comp)
        fR = {1: f1, 2: f2, 3: f3, 4: f4}
        A = en14651_prism(law_envelope(law_mc2010_plain(fct, GFm), law_mc2010_frc_linear(f1, f3, w_end=4.5)),
                          comp, E, cmod_max=3.6, n=90)["fR"]
        B = en14651_prism(law_annexL(f1, f3, fct, E, lcs=125.0, wu=4.5), comp, E, cmod_max=3.6, n=90)["fR"]
        C = identify_law(fR, fct, GFm, E, comp, use=(1, 3))
        D = identify_law(fR, fct, GFm, E, comp, use=(1, 2, 3, 4))
        ident[id_] = D
        rows.append({"id": id_, "fct": fct, "E": E, "exp": fR, "A": A, "B": B, "C": C["sim"]["fR"],
                     "D": D["sim"]["fR"], "Dnodes": (D["nodes_w"].tolist(), D["nodes_s"].tolist()),
                     "Cnodes": (C["nodes_w"].tolist(), C["nodes_s"].tolist())})
        print("  validated", id_, f"fct={fct:.2f}")
    # statistics
    stats = {}
    for m in "ABCD":
        for j in (1, 2, 3, 4):
            r = np.array([row[m][j] / row["exp"][j] for row in rows])
            stats[(m, j)] = (r.mean(), r.std(ddof=1) / r.mean())
            put(f"val{m}{'abcd'[j - 1]}mean", r.mean())
            put(f"val{m}{'abcd'[j - 1]}cov", 100 * r.std(ddof=1) / r.mean(), "{:.1f}")
    # pooled prediction errors for blind points of model C (fR2 and fR4)
    rC = np.array([[row["C"][2] / row["exp"][2], row["C"][4] / row["exp"][4]] for row in rows]).ravel()
    put("valCblindmean", rC.mean())
    put("valCblindcov", 100 * rC.std(ddof=1) / rC.mean(), "{:.1f}")
    put("valCblindmaxerr", 100 * np.max(np.abs(rC - 1)), "{:.0f}")
    nv = np.array([[0.5 * (row["exp"][1] + row["exp"][3]) / row["exp"][2],
                    (row["exp"][3] + 0.5 * (row["exp"][3] - row["exp"][1])) / row["exp"][4]] for row in rows]).ravel()
    put("valNaivemean", nv.mean())
    put("valNaivecov", 100 * nv.std(ddof=1) / nv.mean(), "{:.1f}")
    put("valfctmin", min(r["fct"] for r in rows))
    put("valfctmax", max(r["fct"] for r in rows))
    # table
    with open(os.path.join(TAB, "tab_validation_gao.tex"), "w") as f:
        f.write("\\begin{tabular}{l c cccc cccc cc}\n\\toprule\n")
        f.write(" & & \\multicolumn{4}{c}{$f_{R,j}^{\\mathrm{exp}}$ [MPa]} & \\multicolumn{4}{c}{$f_{R,j}^{\\mathrm{A}}/f_{R,j}^{\\mathrm{exp}}$}"
                " & \\multicolumn{2}{c}{C (blind)}\\\\\n")
        f.write("\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}\\cmidrule(lr){11-12}\n")
        f.write("Mix & $f_{ct}$ & 1 & 2 & 3 & 4 & 1 & 2 & 3 & 4 & 2 & 4\\\\\n\\midrule\n")
        for r in rows:
            e = r["exp"]
            f.write(f"{r['id']} & {r['fct']:.2f} & " + " & ".join(f"{e[j]:.2f}" for j in (1, 2, 3, 4)) + " & "
                    + " & ".join(f"{r['A'][j] / e[j]:.2f}" for j in (1, 2, 3, 4)) + " & "
                    + f"{r['C'][2] / e[2]:.2f} & {r['C'][4] / e[4]:.2f}\\\\\n")
        f.write("\\midrule\n")
        f.write("Mean & & & & & & " + " & ".join(f"{stats[('A', j)][0]:.2f}" for j in (1, 2, 3, 4))
                + f" & {stats[('C', 2)][0]:.2f} & {stats[('C', 4)][0]:.2f}\\\\\n")
        f.write("CoV [\\%] & & & & & & " + " & ".join(f"{100 * stats[('A', j)][1]:.1f}" for j in (1, 2, 3, 4))
                + f" & {100 * stats[('C', 2)][1]:.1f} & {100 * stats[('C', 4)][1]:.1f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    with open(os.path.join(TAB, "tab_validation_summary.tex"), "w") as f:
        f.write("\\begin{tabular}{l cccc}\n\\toprule\n & $f_{R,1}$ & $f_{R,2}$ & $f_{R,3}$ & $f_{R,4}$\\\\\n\\midrule\n")
        labels = {"A": "A: MC2010 linear + plain", "B": "B: EC2 Annex L", "C": "C: hinge inverse (1,3)",
                  "D": "D: hinge inverse (1--4)"}
        for m in "ABCD":
            f.write(labels[m] + " & " + " & ".join(f"{stats[(m, j)][0]:.2f} ({100 * stats[(m, j)][1]:.1f})"
                                                  for j in (1, 2, 3, 4)) + "\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    # figure: predicted vs experimental
    for lang in ("en", "es"):
        fig, axs = plt.subplots(1, 4, figsize=(8.0, 2.4), sharex=True, sharey=True)
        mlab = {"A": T(lang, "A MC2010", "A MC2010"), "B": T(lang, "B Annex L", "B Anejo L"),
                "C": T(lang, "C inverse (1,3)", "C inverso (1,3)")}
        for j, ax in zip((1, 2, 3, 4), axs):
            ex = np.array([r["exp"][j] for r in rows])
            for i, m in enumerate("ABC"):
                pr = np.array([r[m][j] for r in rows])
                ax.plot(ex, pr, ls="none", marker=style.MARKERS[i], color=P[i], label=mlab[m], mfc="none" if i == 2 else P[i])
            ax.plot([0, 8], [0, 8], color="#52514e", lw=0.8)
            ax.set_title(f"$f_{{R,{j}}}$")
            ax.set_xlabel(T(lang, "experimental [MPa]", "experimental [MPa]"))
        axs[0].set_ylabel(T(lang, "predicted [MPa]", "predicho [MPa]"))
        axs[0].set_xlim(0, 8)
        axs[0].set_ylim(0, 8)
        axs[-1].legend(loc="lower right", fontsize=7)
        fig.tight_layout()
        save(fig, "fig_validation_gao", lang)
    # identified laws for a subset
    sel = ["3D65/60-10", "3D65/60-40", "4D65/60-40", "5D65/60-40"]
    w = np.linspace(0, 4.0, 500)
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.0, 3.3))
        for i, k in enumerate(sel):
            ax.plot(w, ident[k]["law"].sigma(w), color=P[i], ls=LS[i], label=k)
        ax.set_xlabel("$w$ [mm]")
        ax.set_ylabel(r"$\sigma$ [MPa]")
        ax.legend(title=T(lang, "identified (fR1-fR4)", "identificada (fR1-fR4)"))
        fig.tight_layout()
        save(fig, "fig_identified_laws", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_gao.json"), "w"), indent=1, default=float)
    return rows


# ======================================================================================
# Validation 2 - size effect on the flexural strength of plain concrete
# ======================================================================================

def validation_size_effect():
    law = law_mc2010_plain(FCTM, GF)
    comp = compression_ec2(FCM, ECM)
    hs = np.array([50, 100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400])
    mor = []
    for h in hs:
        hg = CrackedHinge(Section.rectangle(100, h), law, comp, ECM, n_layers=500)
        r = hg.response(theta_max=40 * FCTM / ECM * hg.s / h * 2, n=140)
        mor.append(r["M"].max() / (100 * h ** 2 / 6) / FCTM)
    mor = np.array(mor)
    mc = np.array([codes.fct_flexural_mc2010(FCTM, h) / FCTM for h in hs])
    err = np.abs(mor / mc - 1)
    put("seMaxErr", 100 * err[(hs >= 100) & (hs <= 1200)].max(), "{:.0f}")
    put("seMaxErrAll", 100 * err.max(), "{:.0f}")
    put("seHmin", hs.min(), "{:d}")
    put("seHmax", hs.max(), "{:d}")
    for h, m, c in zip(hs, mor, mc):
        put(f"seH{int(h)}hinge", m)
        put(f"seH{int(h)}mc", c)
    with open(os.path.join(TAB, "tab_size_effect.tex"), "w") as f:
        f.write("\\begin{tabular}{r r cc c}\n\\toprule\n$h$ [mm] & $h/l_{ch}$ & hinge & MC2010 (5.1-8) & ratio\\\\\n\\midrule\n")
        for h, m, c in zip(hs, mor, mc):
            f.write(f"{h:.0f} & {h / LCH:.2f} & {m:.3f} & {c:.3f} & {m / c:.3f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    hh = np.geomspace(40, 3000, 200)
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.0, 3.3))
        ax.semilogx(hs / LCH, mor, marker="o", color=P[0], ls=LS[0], label=T(lang, "cracked hinge (this work)", "rótula fisurada (este trabajo)"))
        ax.semilogx(hh / LCH, [codes.fct_flexural_mc2010(FCTM, h) / FCTM for h in hh], color=P[1], ls=LS[1],
                    label=T(lang, "MC2010 Eq. (5.1-8)", "MC2010 Ec. (5.1-8)"))
        ax.axhline(1.0, color="#52514e", lw=0.8)
        ax.set_xlabel(r"$h/l_{ch}$")
        ax.set_ylabel(r"$f_{ct,fl}/f_{ct}$")
        ax.legend()
        fig.tight_layout()
        save(fig, "fig_size_effect", lang)


# ======================================================================================
# Example 3 - hinge mechanics: stress profiles (didactic)
# ======================================================================================

def fig_profiles():
    plain = law_mc2010_plain(FCTM, GF)
    env = law_envelope(plain, law_mc2010_frc_linear(FR1K, FR3K, w_end=4.0))
    comp = compression_annexL(FCM, FR1K, ECM)
    hg = CrackedHinge(Section.rectangle(150, 125), env, comp, ECM, n_layers=300)
    thetas = [0.5 * 2 * FCTM / ECM * hg.s / 62.5, 0.001, 0.004, 0.012, 0.04]
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.2), sharey=True)
        y0 = 62.5
        for i, th in enumerate(thetas):
            st = hg.solve(th, 0.0, y0)
            y0 = st.y0
            ax[0].plot(st.sigma, st.y, color=P[i], ls=LS[i], label=fr"$\theta$={th:.4f}")
            ax[1].plot(st.w, st.y, color=P[i], ls=LS[i])
        ax[0].axvline(0, color="#52514e", lw=0.8)
        ax[0].set_xlabel(r"$\sigma$ [MPa]")
        ax[0].set_ylabel("$y$ [mm]")
        ax[0].set_xlim(-30, 6)
        ax[0].legend(fontsize=7)
        ax[1].set_xlabel("$w(y)$ [mm]")
        ax[0].set_title(T(lang, "(a) normal stresses in the hinge", "(a) tensiones normales en la rótula"))
        ax[1].set_title(T(lang, "(b) crack opening profile", "(b) perfil de apertura de fisura"))
        fig.tight_layout()
        save(fig, "fig_profiles", lang)


# ======================================================================================
# Example 4 - hybrid R/SFRC beam: brittle-ductile transition and minimum reinforcement
# ======================================================================================

def hinge_rc(b, h, As, phi, fR1=0.0, fR3=0.0, level="mean", cover=40.0, n_layers=400, fck=FCK):
    fcm = fck + 8
    E = ecm_from_fcm(fcm)
    ft = fctm_from_fck(fck)
    gf = gf_mc2010(fcm)
    plain = law_mc2010_plain(ft, gf)
    if fR1 > 0:
        law = law_envelope(plain, law_annexL(fR1, fR3, ft, E, lcs=125.0, wu=2.5, include_matrix=False))
        comp = compression_annexL(fcm, fR1, E)
    else:
        law, comp = plain, compression_ec2(fcm, E)
    t0, t1 = tcm_bond(ft)
    rb = [Rebar(cover, As, phi, Steel(550.0, 594.0, 200000.0, 0.075), t0, t1)] if As > 0 else []
    return CrackedHinge(Section.rectangle(b, h), law, comp, E, rebars=rb, n_layers=n_layers)


def first_peak_and_post(M, w=None, w_cr=None, w_window=0.2, Mcr0=None):
    """Cracking peak M_cr and post-cracking maximum M_post.

    The cracking peak is the first local maximum of M(theta) reached while the soffit crack
    opening is still small (w <= w_window, default 0.2 mm). If no such local maximum exists the
    response is monotonic through cracking (no snap-back) and M_cr is taken as M at w = w_cr.
    """
    M = np.asarray(M)
    if w is None:
        w = np.zeros_like(M)
    if w_cr is None:
        w_cr = GF / FCTM
    runmax = np.maximum.accumulate(M)
    cond = (M < 0.99 * runmax) & (w <= w_window)  # drop of more than 1 % after a maximum
    if Mcr0 is not None:
        # a cracking peak cannot exceed ~2.2 times the elastic cracking moment f_ct W
        cond &= runmax <= 2.2 * Mcr0
    idx = np.where(cond)[0]
    if idx.size:
        j = int(idx[0])
        return float(runmax[j]), float(M[j:].max())
    Mc = float(np.interp(w_cr, w, M))
    return Mc, float(M.max())


FR1_MIN, FR3_MIN = 2.0, 1.4  # class 2b used in the minimum-reinforcement study


def example4():
    b, h = 300.0, 500.0
    ratios = [0.0, 0.0005, 0.001, 0.0015, 0.0025, 0.005]
    out = {}
    for case, (f1, f3) in {"plain": (0, 0), "sfrc": (FR1_MIN, FR3_MIN)}.items():
        curves = []
        for rho in ratios:
            As = rho * b * (h - 40)
            hg = hinge_rc(b, h, As, 12, f1, f3)
            r = hg.response(theta_max=0.03, n=200)
            curves.append((rho, r))
        out[case] = curves
    for lang in ("en", "es"):
        fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
        for a, case, ttl in zip(axs, ("plain", "sfrc"), (T(lang, "(a) RC (no fibres)", "(a) HA (sin fibras)"),
                                                         T(lang, "(b) R/SFRC, class 2b", "(b) HA+fibras, clase 2b"))):
            for i, (rho, r) in enumerate(out[case]):
                a.semilogx(r["theta"][1:] * 1000, r["M"][1:] / 1e6, color=P[i], ls=LS[i], label=fr"$\rho$={100 * rho:.2f}%")
            a.set_xlabel(r"$\theta$ [mrad] (log)")
            a.set_title(ttl)
            a.set_xlim(0.03, 30)
        axs[0].set_ylabel("$M$ [kNm]")
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        save(fig, "fig_brittle_ductile", lang)
    # minimum reinforcement vs depth: fracture-based criterion vs EC2 / Annex L
    depths = [200, 300, 500, 800, 1200]
    res = []
    for hh in depths:
        row = {"h": hh}
        for case, (f1, f3) in {"plain": (0, 0), "sfrc": (FR1_MIN, FR3_MIN)}.items():
            def margin(rho):
                As = rho * b * (hh - 40)
                hg = hinge_rc(b, hh, As, 12, f1, f3, n_layers=300)
                r = hg.response(theta_max=0.03, n=160)
                Mc, Mp = first_peak_and_post(r["M"], r["w_bottom"], GF / FCTM, Mcr0=FCTM * b * hh ** 2 / 6)
                return Mp / Mc - 1.0
            lo, hi = 0.0, 0.006
            if margin(lo) >= 0:
                row[case] = 0.0
            else:
                for _ in range(16):
                    mid = 0.5 * (lo + hi)
                    if margin(mid) >= 0:
                        hi = mid
                    else:
                        lo = mid
                row[case] = hi
        As_ec, _ = codes.as_min_bending(b, hh, hh - 40, FCTM, 550.0, FCK)
        row["ec2"] = As_ec / (b * (hh - 40))
        As_l, _ = codes.as_min_bending(b, hh, hh - 40, FCTM, 550.0, FCK, fFtu_ef=0.33 * FR3_MIN)
        row["annexL"] = As_l / (b * (hh - 40))
        res.append(row)
        print("  rho_min h=", hh, row)
    with open(os.path.join(TAB, "tab_rhomin.tex"), "w") as f:
        f.write("\\begin{tabular}{r cc cc}\n\\toprule\n & \\multicolumn{2}{c}{RC} & \\multicolumn{2}{c}{R/SFRC 2b}\\\\\n"
                "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n$h$ [mm] & hinge & EC2 (12.1) & hinge & Annex L (L.29)\\\\\n\\midrule\n")
        for r_ in res:
            f.write(f"{r_['h']} & {100 * r_['plain']:.3f} & {100 * r_['ec2']:.3f} & {100 * r_['sfrc']:.3f} & {100 * r_['annexL']:.3f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for r_ in res:
        put(f"rhominH{r_['h']}plain", 100 * r_["plain"], "{:.3f}")
        put(f"rhominH{r_['h']}ec", 100 * r_["ec2"], "{:.3f}")
        put(f"rhominH{r_['h']}sfrc", 100 * r_["sfrc"], "{:.3f}")
        put(f"rhominH{r_['h']}annexL", 100 * r_["annexL"], "{:.3f}")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.0, 3.3))
        H = [r_["h"] for r_ in res]
        ax.plot(H, [100 * r_["plain"] for r_ in res], marker="o", color=P[0], ls=LS[0], label=T(lang, "hinge, RC", "rótula, HA"))
        ax.plot(H, [100 * r_["ec2"] for r_ in res], marker="s", color=P[1], ls=LS[1], label="EC2 (12.1)")
        ax.plot(H, [100 * r_["sfrc"] for r_ in res], marker="^", color=P[2], ls=LS[2], label=T(lang, "hinge, R/SFRC 2b", "rótula, HA+fibras 2b"))
        ax.plot(H, [100 * r_["annexL"] for r_ in res], marker="D", color=P[3], ls=LS[3], label=T(lang, "Annex L (L.29)", "Anejo L (L.29)"))
        ax.set_xlabel("$h$ [mm]")
        ax.set_ylabel(r"$\rho_{min}$ [%]")
        ax.legend()
        fig.tight_layout()
        save(fig, "fig_rhomin", lang)
    return res


# ======================================================================================
# Example 5 - SFRC slab strip without bars (ULS, ECOV, ductility)
# ======================================================================================

def example5():
    b, h = 1000.0, 200.0
    Act = b * h * 1e-6  # whole tensile zone ~ b h (conservative upper bound used for kG)
    kG = codes.kappa_G(b * 0.9 * h * 1e-6)
    put("exfivekG", kG)
    # Annex L rigid-plastic design moment (no bars): equilibrium with parabola-rectangle
    fFtud = kG * 0.33 * FR3K / 1.5
    uls = codes.uls_bending_annexL(b, h, FCK, fFtud, [], 1.5)
    put("exfiveMRdAL", uls["M_Rd"] / 1e6)
    put("exfivexAL", uls["x"], "{:.1f}")
    put("exfivefFtud", fFtud)
    # simple MC2010 formula (7.7-16): m_Rd = fFtud h^2 / 2 (per unit width)
    put("exfiveMRdMC", fFtud * h ** 2 / 2 * b / 1e6)
    # hinge: design (partial factors, Annex L), characteristic and mean (for ECOV)
    # Mean residual strengths (lognormal, as ECOV): fR,m = fR,k exp(1.645 V_fR). V_fR = 0.20 is an ASSUMED
    # coefficient of variation (user input; to be replaced by production-control data).
    VfR = 0.20
    put("exfiveVfR", VfR, "{:.2f}")
    fRm = (FR1K * np.exp(1.645 * VfR), FR3K * np.exp(1.645 * VfR))
    put("exfivefRonem", fRm[0])
    put("exfivefRthreem", fRm[1])
    res = {}
    for lev in ("design", "char", "mean"):
        if lev == "design":
            ft, fac, f1, f3 = FCTK / 1.5, kG / 1.5, FR1K, FR3K
            comp = compression_parabola_rectangle(codes.eta_cc(FCK) * FCK / 1.5, sfrc=True)
        elif lev == "char":
            ft, fac, f1, f3 = FCTK, kG, FR1K, FR3K
            comp = compression_annexL(FCK, FR1K, ECM)
        else:
            ft, fac, f1, f3 = FCTM, kG, fRm[0], fRm[1]
            comp = compression_annexL(FCM, fRm[0], ECM)
        fib = law_annexL(f1, f3, ft, ECM, lcs=125, wu=2.5, include_matrix=False).scaled(fac)
        law = law_envelope(law_mc2010_plain(ft, GF), fib, w_max=2.5)
        hg = CrackedHinge(Section.rectangle(b, h), law, comp, ECM, n_layers=300, w_limit=2.5)
        r = hg.response(theta_max=0.08, n=200)
        res[lev] = r
    for lev in res:
        Mc, Mp = first_peak_and_post(res[lev]["M"], res[lev]["w_bottom"], GF / FCTM)
        put(f"exfive{lev}Mmax", res[lev]["M"].max() / 1e6)
        put(f"exfive{lev}Mcr", Mc / 1e6)
        put(f"exfive{lev}Mpost", Mp / 1e6)
        put(f"exfive{lev}Mwu", res[lev]["M"][-1] / 1e6)
    E = ecov(res["mean"]["M"].max(), res["char"]["M"].max())
    put("exfiveVR", E["V_R"], "{:.3f}")
    put("exfivegR", E["gamma_R"], "{:.3f}")
    put("exfiveRd", E["R_d"] / 1e6)
    # hinge design moment at w_u (rigid-plastic-like) and residual capacity ratio
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.2, 3.3))
        lab = {"design": T(lang, "design laws", "leyes de cálculo"), "char": T(lang, "characteristic", "características"),
               "mean": T(lang, "mean", "medias")}
        for i, lev in enumerate(("mean", "char", "design")):
            ax.plot(res[lev]["w_bottom"], res[lev]["M"] / 1e6, color=P[i], ls=LS[i], label=lab[lev])
        ax.axhline(uls["M_Rd"] / 1e6, color=P[3], ls=LS[3], label=T(lang, "Annex L rigid-plastic $M_{Rd}$", "Anejo L rígido-plástico $M_{Rd}$"))
        ax.axhline(E["R_d"] / 1e6, color=P[6], ls=LS[4], label=T(lang, "ECOV $R_d$ (peak)", "ECOV $R_d$ (pico)"))
        ax.set_xlabel(T(lang, "crack opening at the soffit $w$ [mm]", "apertura en la cara traccionada $w$ [mm]"))
        ax.set_ylabel(T(lang, "$M$ [kNm/m]", "$M$ [kNm/m]"))
        ax.set_xlim(0, 2.6)
        ax.set_ylim(0, None)
        ax.legend(fontsize=6.5, ncol=2, loc="lower center")
        fig.tight_layout()
        save(fig, "fig_slab", lang)


# ======================================================================================
# Example 6 - hybrid beam: ULS (hinge vs Annex L vs ACI) and SLS crack width
# ======================================================================================

def example6():
    b, h, c = 300.0, 500.0, 40.0
    phi = 16.0
    As = 3 * np.pi * phi ** 2 / 4
    d = h - c
    put("exsixAs", As, "{:.0f}")
    put("exsixd", d, "{:.0f}")
    resM = {}
    for case, (f1, f3) in {"RC": (0, 0), "RSFRC": (FR1K, FR3K)}.items():
        hg = hinge_rc(b, h, As, phi, f1, f3, cover=c)
        r = hg.response(theta_max=0.04, n=220)
        resM[case] = r
        put(f"exsix{case}Mmax", r["M"].max() / 1e6)
        put(f"exsix{case}fail", r["failure"] or "none")
        # design-level hinge
        fcd = codes.eta_cc(FCK) * FCK / 1.5
        if f1 > 0:
            fib = law_annexL(f1, f3, FCTK / 1.5, ECM, lcs=125, wu=2.5, include_matrix=False).scaled(1.0 / 1.5)
            law = law_envelope(law_mc2010_plain(FCTK / 1.5, GF), fib, w_max=2.5)
        else:
            law = law_mc2010_plain(FCTK / 1.5, GF)
        comp = compression_parabola_rectangle(fcd, sfrc=f1 > 0)
        t0, t1 = tcm_bond(FCTK)
        rb = [Rebar(c, As, phi, Steel(500 / 1.15, 540 / 1.15, 200000, 0.9 * 0.05), t0, t1)]
        hd = CrackedHinge(Section.rectangle(b, h), law, comp, ECM, rebars=rb, n_layers=400)
        rd = hd.response(theta_max=0.04, n=220)
        put(f"exsix{case}MRdHinge", rd["M"].max() / 1e6)
        fFtud = 0.33 * f3 / 1.5
        uls = codes.uls_bending_annexL(b, h, FCK, fFtud, [(d, As, 500 / 1.15)], 1.5)
        put(f"exsix{case}MRdAL", uls["M_Rd"] / 1e6)
        put(f"exsix{case}xAL", uls["x"], "{:.1f}")
    # ACI 544.4R (Henager & Doherty) for the hybrid beam: hooked fibres 50/1.0? -> use
    # 40 kg/m3 (Vf = 0.51 %), lf/df = 65, Fbe = 1.0 (illustrative fibre data)
    hd = codes.henager_doherty(b, h, d, As, 500.0, FCK, 60.0, 0.9, 0.51, 1.0)
    put("exsixACIMn", hd["Mn"] / 1e6)
    put("exsixACIsigt", hd["sigma_t"], "{:.3f}")
    hd0 = codes.henager_doherty(b, h, d, As, 500.0, FCK, 60.0, 0.9, 0.0, 1.0)
    put("exsixACIMnzero", hd0["Mn"] / 1e6)
    # SLS: crack width at the rebar level vs M
    Ms = np.linspace(40e6, 140e6, 21)
    rho_eff = As / (b * min(2.5 * (h - d), h / 2))
    wh = {}
    wmc = {}
    for case in ("RC", "RSFRC"):
        r = resM[case]
        M = r["M"]
        ws = r["w_s"][:, 0]
        ss = r["sigma_s"][:, 0]
        i_pk = np.argmax(np.diff(M) < 0) if np.any(np.diff(M) < 0) else 0
        # ascending branch after cracking (monotone part used for interpolation)
        mask = np.arange(len(M)) > i_pk
        Mm, wm, sm = M[mask], ws[mask], ss[mask]
        order = np.maximum.accumulate(Mm)
        wh[case] = np.interp(Ms, order, wm)
        sig = np.interp(Ms, order, sm)
        fFtsk = 0.45 * FR1K if case == "RSFRC" else 0.0
        wmc[case] = np.array([codes.crack_width_mc2010_frc(s_, phi, rho_eff, c - phi / 2, FCTM, fFtsk, Ec=ECM)[0]
                              for s_ in sig])
    Mser = 90e6
    for case in ("RC", "RSFRC"):
        put(f"exsix{case}wser", float(np.interp(Mser, Ms, wh[case])), "{:.3f}")
        put(f"exsix{case}wserMC", float(np.interp(Mser, Ms, wmc[case])), "{:.3f}")
    put("exsixMser", Mser / 1e6, "{:.0f}")
    put("exsixrhoeff", 100 * rho_eff, "{:.2f}")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2))
        for i, case in enumerate(("RC", "RSFRC")):
            nm = T(lang, "RC", "HA") if case == "RC" else T(lang, "R/SFRC", "HA+fibras")
            ax[0].plot(resM[case]["theta"] * 1000, resM[case]["M"] / 1e6, color=P[i], ls=LS[i], label=nm)
            ax[1].plot(Ms / 1e6, wh[case], color=P[i], ls=LS[i], marker="o", ms=3, label=nm + T(lang, " hinge", " rótula"))
            ax[1].plot(Ms / 1e6, wmc[case], color=P[i], ls=LS[i + 2], marker="s", mfc="none", ms=3, label=nm + " MC2010 (7.7-22)")
        ax[0].set_xlabel(r"$\theta$ [mrad]")
        ax[0].set_ylabel("$M$ [kNm]")
        ax[0].set_xlim(0, 20)
        ax[0].legend()
        ax[0].set_title(T(lang, "(a) moment-rotation (mean)", "(a) momento-giro (medios)"))
        ax[1].set_xlabel("$M$ [kNm]")
        ax[1].set_ylabel(T(lang, "crack width at bar level [mm]", "abertura de fisura en la armadura [mm]"))
        ax[1].legend(fontsize=6.5)
        ax[1].set_title(T(lang, "(b) SLS crack width", "(b) abertura de fisura en ELS"))
        fig.tight_layout()
        save(fig, "fig_hybrid_beam", lang)


# ======================================================================================
# Example 7 - shear without stirrups (code formulae)
# ======================================================================================

def example7():
    d = 460.0
    fR3 = np.linspace(0, 6, 61)
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.0, 3.3))
        for i, rho in enumerate([0.005, 0.01, 0.02]):
            tau = [codes.tau_RdcF_annexL(rho, FCK, 0.9 * d, 0.33 * f / 1.5, 0.33 * f)[0] for f in fR3]
            ax.plot(fR3, tau, color=P[i], ls=LS[i], label=fr"Annex L, $\rho_l$={100 * rho:.1f}%" if lang == "en" else fr"Anejo L, $\rho_l$={100 * rho:.1f}%")
            vmc = [codes.VRdF_mc2010(1.0, d, rho, FCK, (0.5 * f - 0.2 * (f / 0.9)) if f > 0 else 0.0, FCTK) / d for f in fR3]
            ax.plot(fR3, vmc, color=P[i], ls=LS[i + 3], lw=1.0, label=fr"MC2010 (7.7-5), $\rho_l$={100 * rho:.1f}%")
        ax.set_xlabel("$f_{R,3k}$ [MPa]")
        ax.set_ylabel(T(lang, "design shear stress [MPa]", "tensión tangencial de cálculo [MPa]"))
        ax.legend(fontsize=6.5)
        fig.tight_layout()
        save(fig, "fig_shear", lang)
    t, eta = codes.tau_RdcF_annexL(0.01, FCK, 0.9 * d, 0.33 * FR3K / 1.5, 0.33 * FR3K)
    put("exsevtau", t)
    put("exseveta", eta)
    put("exsevtauRC", codes.tau_Rdc(0.01, FCK, 0.9 * d))
    put("exsevtaumin", codes.tau_Rdc_min(FCK, 435.0, 0.9 * d))
    put("exsevddg", codes.ddg(16.0, FCK), "{:.0f}")
    # stirrup replacement checks: L.32 fFtu,ef >= 0.08 sqrt(fck); MC2010 (7.7-14) fFtuk >= 0.08 sqrt(fck)
    put("exsevLthirtytwo", 0.08 * np.sqrt(FCK))
    put("exsevfFtuef", 0.33 * FR3K)


# ======================================================================================
# LEFM figure (course) - SENB geometry factor and nominal strength
# ======================================================================================

def fig_lefm():
    a = np.linspace(0.05, 0.8, 200)
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(4.6, 3.1))
        ax.plot(a, [fracture.senb_geometry_factor(x) for x in a], color=P[0])
        ax.set_xlabel(r"$\alpha=a/W$")
        ax.set_ylabel(r"$F(\alpha)$")
        ax.set_title(T(lang, "SENB, S/W = 4 (Srawley, 1976)", "Probeta SENB, S/W = 4 (Srawley, 1976)"))
        fig.tight_layout()
        save(fig, "fig_lefm_F", lang)
    KIC = fracture.kic_from_gf(GF, ECM)
    put("KIC", KIC, "{:.1f}")
    put("KICSI", KIC / np.sqrt(1000.0), "{:.2f}")


def validation_bosco():
    """Member-level validation: lightly reinforced beams of Bosco, Carpinteri, Debernardi (1990)."""
    rows = []
    for (bid, g, b, h, L, n, phi, rho, NP, Pc, Py, Pu) in BOSCO_BEAMS:
        G = BOSCO_GRADES[g]
        E, fcm, GFb = G["Ec"], G["fcm"], G["GF"]
        ft = fctm_from_fck(fcm - 8.0)
        As = n * np.pi * phi ** 2 / 4.0
        fy = NP * G["KIC"] / (rho / 100.0 * np.sqrt(h))  # back-calculated with Eq. (15)
        t0, t1 = tcm_bond(ft)
        rb = [Rebar(0.1 * h, As, phi, Steel(fy, 1.08 * fy, 200000.0, 0.05), t0, t1)]
        hg = CrackedHinge(Section.rectangle(b, h), law_mc2010_plain(ft, GFb), compression_ec2(fcm, E), E,
                          rebars=rb, n_layers=300)
        r = hg.response(theta_max=0.03, n=160)
        Msw = 24e-6 * b * h * L ** 2 / 8.0  # self-weight 24 kN/m3
        M, w = r["M"], r["w_bottom"]
        w1 = GFb / ft
        Mcr = M[w <= 2.0 * w1].max()
        iy = np.where(r["sigma_s"][:, 0] >= fy - 1e-6)[0]
        My = M[iy[0]] if iy.size else M.max()
        Pcr = 4.0 * (Mcr - Msw) / L / 1000.0
        Pym = 4.0 * (My - Msw) / L / 1000.0
        rows.append({"id": bid, "grade": g, "h": h, "NP": NP, "fy": fy, "Pc": Pc, "Py": Py,
                     "Pcm": Pcr, "Pym": Pym, "bexp": Py < Pc, "bmod": Pym < Pcr})
    rc = np.array([r["Pcm"] / r["Pc"] for r in rows])
    ry = np.array([r["Pym"] / r["Py"] for r in rows])
    agree = sum(r["bexp"] == r["bmod"] for r in rows)
    put("boN", len(rows), "{:d}")
    put("boPcmean", rc.mean())
    put("boPccov", 100 * rc.std(ddof=1) / rc.mean(), "{:.1f}")
    put("boPymean", ry.mean())
    put("boPycov", 100 * ry.std(ddof=1) / ry.mean(), "{:.1f}")
    put("boAgree", agree, "{:d}")
    put("boNexpbrittle", sum(r["bexp"] for r in rows), "{:d}")
    put("boNmodbrittle", sum(r["bmod"] for r in rows), "{:d}")
    put("boConservative", sum((not r["bexp"]) and r["bmod"] for r in rows), "{:d}")
    put("boUnsafe", sum(r["bexp"] and (not r["bmod"]) for r in rows), "{:d}")
    with open(os.path.join(TAB, "tab_validation_bosco.tex"), "w") as f:
        f.write("\\begin{tabular}{l c r c r cc cc cc}\n\\toprule\n")
        f.write("Beam & Gr. & $h$ & $N_P$ & $f_y$ & $P_{cr}^{exp}$ & $P_{cr}^{mod}$ & $P_{y}^{exp}$ & $P_{y}^{mod}$ & exp & mod\\\\\n")
        f.write(" & & [mm] & & [MPa] & [kN] & [kN] & [kN] & [kN] & \\multicolumn{2}{c}{brittle?}\\\\\n\\midrule\n")
        for r in rows:
            f.write(f"{r['id']} & {r['grade']} & {r['h']} & {r['NP']:.3f} & {r['fy']:.0f} & {r['Pc']:.2f} & {r['Pcm']:.2f} & "
                    f"{r['Py']:.2f} & {r['Pym']:.2f} & {'Y' if r['bexp'] else 'N'} & {'Y' if r['bmod'] else 'N'}\\\\\n")
        f.write("\\midrule\n")
        f.write(f"\\multicolumn{{5}}{{l}}{{mean ratio mod/exp (CoV \\%)}} & \\multicolumn{{2}}{{c}}{{{rc.mean():.2f} ({100 * rc.std(ddof=1) / rc.mean():.1f})}} & "
                f"\\multicolumn{{2}}{{c}}{{{ry.mean():.2f} ({100 * ry.std(ddof=1) / ry.mean():.1f})}} & \\multicolumn{{2}}{{c}}{{{agree}/{len(rows)}}}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2))
        for i, (key, kexp, ttl) in enumerate((("Pcm", "Pc", T(lang, "(a) cracking load", "(a) carga de fisuración")),
                                              ("Pym", "Py", T(lang, "(b) yielding load", "(b) carga de plastificación")))):
            for j, gname in enumerate(("I", "II", "III")):
                xs = [r[kexp] for r in rows if r["grade"] == gname]
                ys = [r[key] for r in rows if r["grade"] == gname]
                ax[i].plot(xs, ys, ls="none", marker=style.MARKERS[j], color=P[j], label=T(lang, "grade ", "hormigón ") + gname)
            m = max(max(r[kexp] for r in rows), max(r[key] for r in rows)) * 1.1
            ax[i].plot([0, m], [0, m], color="#52514e", lw=0.8)
            ax[i].set_xlim(0, m)
            ax[i].set_ylim(0, m)
            ax[i].set_xlabel(T(lang, "experimental [kN]", "experimental [kN]"))
            ax[i].set_title(ttl)
        ax[0].set_ylabel(T(lang, "cracked hinge [kN]", "rótula cohesiva [kN]"))
        ax[0].legend()
        fig.tight_layout()
        save(fig, "fig_validation_bosco", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_bosco.json"), "w"), indent=1, default=float)


def _stats(x):
    x = np.asarray(x, dtype=float)
    return x.mean(), 100.0 * x.std(ddof=1) / x.mean()


def validation_carpinteri():
    """R/SFRC beams of Ventura et al. / Carpinteri, Cadamuro & Ventura (2015): 4 fibre contents x 5 bar ratios."""
    rows = []
    Msw = 24e-6 * CARP.B * CARP.H * CARP.L ** 2 / 8.0
    for var in ("A", "B"):
        for vf, m in CARP.MAT.items():
            ft, GFv = (m["ft8"], m["GF8"]) if var == "A" else (m["ft3"], m["GF3"])
            for lab, (n, phi, fy, fu) in CARP.BARS.items():
                rb = []
                if n:
                    t0, t1 = tcm_bond(fctm_from_fck(m["fc"] - 8.0))
                    # Table 8 of the paper: elastic-perfectly plastic steel with fy = fu
                    rb = [Rebar(CARP.C, n * np.pi * phi ** 2 / 4.0, phi, Steel(0.999 * fu, fu, CARP.ES, 0.10), t0, t1)]
                comp = compression_ec2(m["fc"], CARP.EC)
                hg = CrackedHinge(Section.rectangle(CARP.B, CARP.H), law_linear_softening(ft, GFv), comp, CARP.EC,
                                  rebars=rb, n_layers=300)
                r = hg.response(theta_max=0.06, n=150, stop_on_failure=False)
                ok = np.cumprod(-r["eps_top"] <= comp.eps_cu).astype(bool)
                M = r["M"][ok]
                k = int(np.argmax(M))
                if n and r["sigma_s"][ok][k, 0] >= 0.99 * fu:
                    mode = "Y"   # bars at f_u (plateau)
                elif -r["eps_top"][ok][k] >= 0.95 * comp.eps_cu:
                    mode = "C"   # concrete crushing limits the peak
                else:
                    mode = "P"   # cracking peak governs (no bars, or bars below f_u)
                Pm = 4.0 * (M[k] - Msw) / CARP.L / 1000.0
                pe, pmin, pmax = CARP.PEXP[vf][lab]
                rows.append(dict(var=var, vf=vf, bars=lab, rho=100 * n * np.pi * phi ** 2 / 4 / (CARP.B * CARP.H),
                                 Pexp=pe, Pmin=pmin, Pmax=pmax, Pmod=Pm, mode=mode, ratio=Pm / pe,
                                 inrange=pmin <= Pm <= pmax))
    for var in ("A", "B"):
        rr = [r["ratio"] for r in rows if r["var"] == var]
        mu, cov = _stats(rr)
        put(f"ca{var}mean", mu)
        put(f"ca{var}cov", cov, "{:.1f}")
        put(f"ca{var}N", len(rr), "{:d}")
        put(f"ca{var}inrange", sum(r["inrange"] for r in rows if r["var"] == var), "{:d}")
        rb = [r["ratio"] for r in rows if r["var"] == var and r["bars"] != "0Ø0"]
        mu, cov = _stats(rb)
        put(f"ca{var}meanRC", mu)
        put(f"ca{var}covRC", cov, "{:.1f}")
        rp = [r["ratio"] for r in rows if r["var"] == var and r["bars"] == "0Ø0"]
        put(f"ca{var}meanPlain", np.mean(rp))
        rx = [r["ratio"] for r in rows if r["var"] == var and r["bars"] not in ("0Ø0", "2Ø20")]
        mu, cov = _stats(rx)
        put(f"ca{var}meanMid", mu)
        put(f"ca{var}covMid", cov, "{:.1f}")
        r20 = [r["ratio"] for r in rows if r["var"] == var and r["bars"] == "2Ø20"]
        put(f"ca{var}minOver", min(r20))
        put(f"ca{var}maxOver", max(r20))
    with open(os.path.join(TAB, "tab_validation_carpinteri.tex"), "w") as f:
        f.write("\\begin{tabular}{r l r ccc c cc c}\n\\toprule\n")
        f.write("$V_f$ & bars & $\\rho$ & $P^{exp}$ & $P^{exp}_{min}$ & $P^{exp}_{max}$ & $P^{mod}_A$ & $P^{mod}_B$ & $P_B/P^{exp}$ & mode\\\\\n")
        f.write("{[kg/m$^3$]} & & [\\%] & [kN] & [kN] & [kN] & [kN] & [kN] & & \\\\\n\\midrule\n")
        for a, b in zip([r for r in rows if r["var"] == "A"], [r for r in rows if r["var"] == "B"]):
            f.write(f"{a['vf']} & {a['bars'].replace('Ø', '$\\phi$')} & {a['rho']:.2f} & {a['Pexp']:.1f} & {a['Pmin']:.1f} & {a['Pmax']:.1f} & "
                    f"{a['Pmod']:.1f} & {b['Pmod']:.1f} & {b['ratio']:.2f} & {b['mode']}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
        for i, var in enumerate(("A", "B")):
            for j, vf in enumerate((0, 10, 20, 40)):
                sel = [r for r in rows if r["var"] == var and r["vf"] == vf]
                x = np.array([r["Pexp"] for r in sel])
                y = np.array([r["Pmod"] for r in sel])
                xerr = np.array([[r["Pexp"] - r["Pmin"] for r in sel], [r["Pmax"] - r["Pexp"] for r in sel]])
                ax[i].errorbar(x, y, xerr=xerr, ls="none", marker=style.MARKERS[j], color=P[j], ms=5, lw=0.8,
                               capsize=2, label=f"{vf} kg/m$^3$")
            ax[i].plot([8, 200], [8, 200], color="#52514e", lw=0.8)
            for fct in (0.8, 1.2):
                ax[i].plot([8, 200], [8 * fct, 200 * fct], color="#a3a29d", lw=0.6, ls="--")
            ax[i].set_xscale("log")
            ax[i].set_yscale("log")
            ax[i].set_xlim(8, 200)
            ax[i].set_ylim(8, 200)
            ax[i].set_xlabel(T(lang, "experimental peak load [kN]", "carga máxima experimental [kN]"))
            ax[i].set_title(T(lang, "(a) inputs A: $f_t$ of Table 8", "(a) datos A: $f_t$ de la Tabla 8") if var == "A"
                            else T(lang, "(b) inputs B: $f_t$ of Table 3", "(b) datos B: $f_t$ de la Tabla 3"))
        ax[0].set_ylabel(T(lang, "cracked hinge [kN]", "rótula cohesiva [kN]"))
        ax[0].legend(fontsize=7)
        fig.tight_layout()
        save(fig, "fig_validation_carpinteri", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_carpinteri.json"), "w"), indent=1, default=float)


def validation_minelli():
    """Flexural capacity of the 9 full-scale beams of Minelli et al. (2014) by the cracked hinge vs Mu = Vu a."""
    rows = []
    for (H, mix), (fm, Vu, Mu, Mfl) in MINELLI.TESTS.items():
        m = MINELLI.MAT[mix]
        d, L, a, phi, As, fy, fu = MINELLI.GEOM[H]
        plain = law_mc2010_plain(m["fctm"], gf_mc2010(m["fcm"]))
        law = plain if m["fR"] is None else law_envelope(plain, law_mc2010_frc_linear(m["fR"][0], m["fR"][2]))
        t0, t1 = tcm_bond(m["fctm"])
        rb = [Rebar(H - d, As, phi, Steel(fy, fu, 200000.0, 0.08), t0, t1)]
        hg = CrackedHinge(Section.rectangle(250.0, float(H)), law, compression_ec2(m["fcm"], m["Ecm"]), m["Ecm"],
                          rebars=rb, n_layers=300)
        r = hg.response(theta_max=0.03, n=120)
        Mm = r["M"].max() / 1e6
        rows.append(dict(H=H, mix=mix, fm=fm, Vu=Vu, Mu=Mu, Mfl=Mfl, Mfmsd=Mm, r=Mu / Mm, rfl=Mfl / Mm,
                         Vflex=Mm / (a / 1000.0)))
    put("miNflex", sum(r["fm"] == "SY" for r in rows), "{:d}")
    put("miRflexmin", min(r["r"] for r in rows if r["fm"] == "SY"))
    put("miRflexmax", max(r["r"] for r in rows if r["fm"] == "SY"))
    put("miRshearmax", max(r["r"] for r in rows if r["fm"] == "S"))
    put("miRshearmin", min(r["r"] for r in rows if r["fm"] == "S"))
    mu, cov = _stats([r["rfl"] for r in rows])
    put("miMflmean", mu)
    put("miMflcov", cov, "{:.1f}")
    with open(os.path.join(TAB, "tab_validation_minelli.tex"), "w") as f:
        f.write("\\begin{tabular}{l l c r r r r c}\n\\toprule\n")
        f.write("Beam & mix & mode & $V_u$ & $M_u=V_u a$ & $M_{u,fl}$ (MC2010) & $M_{max}$ (FMSD) & $M_u/M_{max}$\\\\\n")
        f.write(" & & & [kN] & [kNm] & [kNm] & [kNm] & \\\\\n\\midrule\n")
        for r in rows:
            f.write(f"H{r['H']} & {r['mix']} & {r['fm']} & {r['Vu']} & {r['Mu']} & {r['Mfl']} & {r['Mfmsd']:.0f} & {r['r']:.2f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    json.dump(rows, open(os.path.join(TAB, "validation_minelli.json"), "w"), indent=1, default=float)
    return rows


def validation_cuenca_shear():
    """Mean-level check of MC2010 (7.7-5) and Annex L (L.19) against the Cuenca et al. (2018) database."""
    rows = []
    for (ref, bid, fib, vf, bw, d, ad, rho, fc, fR1, fR3, fm, Vu) in CUENCA_DB:
        if fm not in ("S", "SY"):
            continue
        rho = rho / 100.0
        fctm = fctm_from_fck(max(fc - 8.0, 12.0))
        fFts = 0.45 * fR1
        fFtu15 = max(fFts - 1.5 / 2.5 * (fFts - 0.5 * fR3 + 0.2 * fR1), 0.0)  # MC2010 (5.6-6), wu = 1.5 mm
        Vmc = codes.VRdF_mc2010(bw, d, rho, fc, fFtu15, fctm, gamma_c=1.0) / 1e3
        tau, _ = codes.tau_RdcF_annexL(rho, fc, 0.9 * d, 0.33 * fR3, 0.33 * fR3, D_lower=16.0, gamma_V=1.0)
        Val = tau * bw * 0.9 * d / 1e3
        rows.append(dict(ref=ref, id=bid, bw=bw, d=d, rho=100 * rho, fc=fc, fR1=fR1, fR3=fR3, Vu=Vu,
                         Vmc=Vmc, Val=Val, rmc=Vu / Vmc, ral=Vu / Val))
    for key, nm in (("rmc", "MC"), ("ral", "AL")):
        mu, cov = _stats([r[key] for r in rows])
        put(f"cu{nm}mean", mu)
        put(f"cu{nm}cov", cov, "{:.1f}")
        put(f"cu{nm}unsafe", sum(r[key] < 1.0 for r in rows), "{:d}")
        put(f"cu{nm}min", min(r[key] for r in rows))
        big = [r[key] for r in rows if r["d"] >= 900]
        put(f"cu{nm}big", np.mean(big))
    put("cuN", len(rows), "{:d}")
    put("cuNall", len(CUENCA_DB), "{:d}")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2))
        for i, (key, ttl) in enumerate((("rmc", T(lang, "(a) MC2010 Eq. (7.7-5)", "(a) MC2010 Ec. (7.7-5)")),
                                        ("ral", T(lang, "(b) EN 1992-1-1 Annex L (L.19)", "(b) EN 1992-1-1 Anejo L (L.19)")))):
            x = [r["d"] for r in rows]
            y = [r[key] for r in rows]
            c = [P[0] if r["fR3"] < 3 else (P[1] if r["fR3"] < 5 else P[2]) for r in rows]
            ax[i].scatter(x, y, c=c, s=22, edgecolors="white", linewidths=0.6, zorder=3)
            ax[i].axhline(1.0, color="#52514e", lw=0.8)
            ax[i].set_xscale("log")
            ax[i].set_xlabel("$d$ [mm]")
            ax[i].set_title(ttl)
        ax[0].set_ylabel("$V_{u,exp}/V_{pred}$ (" + T(lang, "mean values", "valores medios") + ", $\\gamma=1$)")
        for j, lb in enumerate(("$f_{R3}<3$", "$3\\leq f_{R3}<5$", "$f_{R3}\\geq5$ MPa")):
            ax[1].scatter([], [], c=P[j], s=22, label=lb)
        ax[1].legend(fontsize=7)
        fig.tight_layout()
        save(fig, "fig_validation_cuenca", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_cuenca.json"), "w"), indent=1, default=float)


def validation_tiberti():
    """81 EN 14651 series of Tiberti et al. (2017): (i) blind f_L from the measured uniaxial f_ct;
    (ii) consistency of the MC2010 and Annex L laws; (iii) blind fR2, fR4 from inverse analysis on
    fR1, fR3; (iv) scatter of the residual strengths (input of the ECOV format)."""
    rows = []
    for r in TIB.DB:
        E = r["Ec"] * 1000.0
        fcm = 0.8 * r["fcube"]
        GFm = gf_mc2010(fcm)
        comp = compression_elastic_plastic(fcm, E)
        f1, f2, f3, f4 = r["fR"]
        fR = {1: f1, 2: f2, 3: f3, 4: f4}
        row = {"id": r["id"], "vertical": r["id"] in TIB.VERTICAL, "fct_measured": r["fct_measured"],
               "fctm": r["fctm"], "fLexp": r["fL"], "exp": fR, "cv": r["cv"], "Vf": r["Vf"]}
        if r["fct_measured"]:
            pr = en14651_prism(law_annexL(f1, f3, r["fctm"], E, lcs=125.0, wu=4.5), comp, E, cmod_max=0.1, n=60)
            row["fLpred"] = pr["fL"]
            pb = en14651_prism(law_mc2010_plain(r["fctm"], GFm), comp, E, cmod_max=0.1, n=60)
            row["fLblind"] = pb["fL"]
        fct = identify_fct(r["fL"] * 2.0 * 150.0 * 125.0 ** 2 / (3.0 * 500.0), E, GFm, comp)
        row["fct_id"] = fct
        row["A"] = en14651_prism(law_envelope(law_mc2010_plain(fct, GFm), law_mc2010_frc_linear(f1, f3, w_end=4.5)),
                                 comp, E, cmod_max=3.6, n=80)["fR"]
        row["B"] = en14651_prism(law_annexL(f1, f3, fct, E, lcs=125.0, wu=4.5), comp, E, cmod_max=3.6, n=80)["fR"]
        row["C"] = identify_law(fR, fct, GFm, E, comp, use=(1, 3), n=80)["sim"]["fR"]
        rows.append(row)
    hor = [r for r in rows if not r["vertical"]]
    put("tiN", len(rows), "{:d}")
    put("tiNh", len(hor), "{:d}")
    put("tiNbeams", sum(r["N"] for r in TIB.DB), "{:d}")
    # (i) blind fL
    rl = np.array([r["fLpred"] / r["fLexp"] for r in rows if "fLpred" in r])
    put("tiNfL", len(rl), "{:d}")
    mu, cov = _stats(rl)
    put("tifLmean", mu)
    put("tifLcov", cov, "{:.1f}")
    rb_ = np.array([r["fLblind"] / r["fLexp"] for r in rows if "fLblind" in r])
    mu, cov = _stats(rb_)
    put("tifLblindmean", mu)
    put("tifLblindcov", cov, "{:.1f}")
    # (ii)-(iii)
    for m in "ABC":
        for j in (1, 2, 3, 4):
            x = np.array([r[m][j] / r["exp"][j] for r in hor])
            mu, cov = _stats(x)
            put(f"ti{m}{'abcd'[j - 1]}mean", mu)
            put(f"ti{m}{'abcd'[j - 1]}cov", cov, "{:.1f}")
    rC = np.array([[r["C"][2] / r["exp"][2], r["C"][4] / r["exp"][4]] for r in hor]).ravel()
    mu, cov = _stats(rC)
    put("tiCblindmean", mu)
    put("tiCblindcov", cov, "{:.1f}")
    put("tiCblindwithin", 100.0 * np.mean(np.abs(rC - 1.0) <= 0.10), "{:.0f}")
    nv = np.array([[0.5 * (r["exp"][1] + r["exp"][3]) / r["exp"][2],
                    (r["exp"][3] + 0.5 * (r["exp"][3] - r["exp"][1])) / r["exp"][4]] for r in hor]).ravel()
    mu, cov = _stats(nv)
    put("tiNaivemean", mu)
    put("tiNaivecov", cov, "{:.1f}")
    put("tiNaivewithin", 100.0 * np.mean(np.abs(nv - 1.0) <= 0.10), "{:.0f}")
    # (iv) scatter of residual strengths
    cv1 = np.array([r["cv"][1] for r in hor])
    cv3 = np.array([r["cv"][3] for r in hor])
    put("ticvmeanseries", np.median(cv3) / np.sqrt(np.median([TIB.DB[i]["N"] for i in range(len(TIB.DB))])))
    put("ticvonemed", np.median(cv1))
    put("ticvthreemed", np.median(cv3))
    put("ticvthreepninety", np.percentile(cv3, 90))
    put("ticvthreele", 100.0 * np.mean(cv3 <= 0.20), "{:.0f}")
    cvv = np.array([r["cv"][3] for r in rows if r["vertical"]])
    put("ticvvertmin", cvv.min())
    put("ticvvertmax", cvv.max())
    with open(os.path.join(TAB, "tab_validation_tiberti.tex"), "w") as f:
        f.write("\\begin{tabular}{l cccc}\n\\toprule\n & $f_{R,1}$ & $f_{R,2}$ & $f_{R,3}$ & $f_{R,4}$\\\\\n\\midrule\n")
        labels = {"A": "A: MC2010 linear + plain", "B": "B: EC2 Annex L", "C": "C: hinge inverse (1,3)"}
        for m in "ABC":
            cells = []
            for j in (1, 2, 3, 4):
                cells.append(NUM[f"ti{m}{'abcd'[j - 1]}mean"] + " (" + NUM[f"ti{m}{'abcd'[j - 1]}cov"] + ")")
            f.write(labels[m] + " & " + " & ".join(cells) + "\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 3, figsize=(8.0, 2.8))
        x = [r["fLexp"] for r in rows if "fLpred" in r]
        y = [r["fLpred"] for r in rows if "fLpred" in r]
        ax[0].plot(x, y, ls="none", marker="o", color=P[0], ms=4, mfc="none")
        ax[0].plot([2, 9], [2, 9], color="#52514e", lw=0.8)
        ax[0].set_xlim(2, 9)
        ax[0].set_ylim(2, 9)
        ax[0].set_xlabel(T(lang, "$f_L$ experimental [MPa]", "$f_L$ experimental [MPa]"))
        ax[0].set_ylabel(T(lang, "$f_L$ hinge (measured $f_{ct}$) [MPa]", "$f_L$ rótula ($f_{ct}$ medida) [MPa]"))
        ax[0].set_title(T(lang, "(a) blind $f_L$", "(a) $f_L$ ciego"))
        for j, mk in ((2, "o"), (4, "s")):
            ax[1].plot([r["exp"][j] for r in hor], [r["C"][j] for r in hor], ls="none", marker=mk,
                       color=P[1 if j == 2 else 2], ms=4, label=f"$f_{{R,{j}}}$")
        ax[1].plot([0, 12], [0, 12], color="#52514e", lw=0.8)
        ax[1].set_xlim(0, 12)
        ax[1].set_ylim(0, 12)
        ax[1].legend(fontsize=7)
        ax[1].set_xlabel(T(lang, "experimental [MPa]", "experimental [MPa]"))
        ax[1].set_title(T(lang, "(b) blind $f_{R,2}$, $f_{R,4}$", "(b) $f_{R,2}$, $f_{R,4}$ ciegos"))
        bins = np.arange(0.0, 0.50, 0.025)
        ax[2].hist(cv3, bins=bins, color=P[0], edgecolor="white", label=T(lang, "horizontal", "horizontal"))
        ax[2].hist(cvv, bins=bins, color=P[1], edgecolor="white", label=T(lang, "vertical", "vertical"))
        ax[2].axvline(0.20, color="#52514e", lw=0.8, ls="--")
        ax[2].set_xlabel("CoV $f_{R,3}$")
        ax[2].set_ylabel(T(lang, "series", "series"))
        ax[2].legend(fontsize=7)
        ax[2].set_title(T(lang, "(c) scatter (ECOV input)", "(c) dispersión (dato del ECOV)"))
        fig.tight_layout()
        save(fig, "fig_validation_tiberti", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_tiberti.json"), "w"), indent=1, default=float)


def derived_numbers():
    """Percentages quoted in the text, derived from previously computed macros (traceability)."""
    f = lambda k: float(NUM[k])
    put("pcAone", 100 * (f("extwoAfRa") / FR1K - 1), "{:.0f}")
    put("pcAthree", 100 * (f("extwoAfRc") / FR3K - 1), "{:.0f}")
    put("pcBone", 100 * (f("extwoBfRa") / FR1K - 1), "{:.0f}")
    put("pcBthree", 100 * (f("extwoBfRc") / FR3K - 1), "{:.0f}")
    put("pcCthree", 100 * (f("extwoCfRc") / FR3K - 1), "{:.0f}")
    put("pcvalBc", 100 * (f("valBcmean") - 1), "{:.0f}")
    lo = 100 * (min(f("extwoAfRa") / FR1K, f("valAamean"), f("tiAamean")) - 1)
    hi = 100 * (max(f("extwoAfRa") / FR1K, f("valAamean"), f("tiAamean")) - 1)
    put("pcMCoverlo", lo, "{:.0f}")
    put("blindLo", min(f("valCblindmean"), f("tiCblindmean")))
    put("blindHi", max(f("valCblindmean"), f("tiCblindmean")))
    put("pcmDiffLo", 100 * (min(f("pcmRfrcLo"), f("pcmRplaOne")) - 1), "{:.0f}")
    put("pcmDiffHi", 100 * (max(f("pcmRfrcHi"), f("pcmRplaHalf")) - 1), "{:.0f}")
    put("pcMCoverhi", hi, "{:.0f}")


def compare_pcm():
    """Model-to-model comparison of the flexural strength of a section with linear softening:
    cracked hinge (s = h/2) vs planar crack model, Eq. (20) of Carmona et al. (2022), and the
    FEM-based level-III expression of Planas, Guinea & Elices (Bazant & Planas 1998, Eq. 9.3.12)."""
    E, f, h, b = 30000.0, 3.0, 200.0, 100.0
    betas = np.logspace(-3, np.log10(3.0), 16)  # layer stability with s = h/2 requires beta_H < 4
    hin = []
    for beta in betas:
        lch = h / beta
        wu = 2.0 * f * lch / E  # linear softening: l_ch,FRC = E wu / (2 fFts)
        law = CohesiveLaw([0.0, wu], [f, 0.0])
        hg = CrackedHinge(Section.rectangle(b, h), law, compression_elastic_plastic(500.0, E), E, n_layers=400)
        r = hg.response(theta_max=2.4 * wu / h, n=220, stop_on_failure=False)
        hin.append(6.0 * r["M"].max() / (b * h * h) / f)
    hin = np.array(hin)
    pcm = 1.0 + 1.0 / (0.5 + 4.3 * np.sqrt(betas))
    pla = 1.0 + 1.0 / (0.85 + 2.3 * np.sqrt(betas))
    def at(x, arr):
        return float(np.interp(np.log(x), np.log(betas), arr))
    put("pcmRfrcLo", at(0.001, hin) / at(0.001, pcm))
    put("pcmRfrcHi", at(0.01, hin) / at(0.01, pcm))
    put("pcmRplaHalf", at(0.5, hin) / at(0.5, pla))
    put("pcmRplaOne", at(1.0, hin) / at(1.0, pla))
    put("pcmRplaThree", at(3.0, hin) / at(3.0, pla))
    put("pcmRpcmPlaOne", at(1.0, pcm) / at(1.0, pla))
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(5.2, 3.3))
        ax.semilogx(betas, hin, color=P[0], marker="o", ms=3.5, label=T(lang, "cracked hinge, $s=h/2$ (this work)", "rótula fisurada, $s=h/2$ (este trabajo)"))
        ax.semilogx(betas, pcm, color=P[1], ls="--", label=T(lang, "planar crack model, Eq. (20) of [Carmona 2022]", "modelo de fisura plana, Ec. (20) de [Carmona 2022]"))
        m = betas >= 0.1
        ax.semilogx(betas[m], pla[m], color=P[2], ls="-.", label=T(lang, "FEM fit, Planas et al. (level III, $\\beta_H\\geq0.1$)", "ajuste MEF, Planas et al. (nivel III, $\\beta_H\\geq0{,}1$)"))
        ax.axhline(3.0, color="#a3a29d", lw=0.6, ls=":")
        ax.axhline(1.0, color="#a3a29d", lw=0.6, ls=":")
        ax.axvspan(0.001, 0.01, color=P[3], alpha=0.12)
        ax.text(0.0013, 1.12, T(lang, "usual FRC", "HRF usual"), fontsize=7)
        ax.set_xlabel(r"$\beta_H=h/l_{ch,FRC}$")
        ax.set_ylabel(r"$f_R^*=6M_{max}/(bh^2f_{Fts})$")
        ax.set_ylim(0.9, 3.1)
        ax.legend(fontsize=6.5)
        fig.tight_layout()
        save(fig, "fig_compare_pcm", lang)


def validation_delarosa2025():
    """Comparison on the 9 R/SFRC combinations used by De La Rosa, Ruiz & Carmona (2025):
    FMSD (inputs B and inputs C = MC2010 fR-based law with the fR of their Table 1), their
    planar-crack model (Table 2, eta = 4) and the experimental moments with L = 1200 mm."""
    carp = {(r["vf"], r["bars"]): r for r in json.load(open(os.path.join(TAB, "validation_carpinteri.json")))
            if r["var"] == "B"}
    Msw = 24e-6 * CARP.B * CARP.H * CARP.L ** 2 / 8.0
    rows = []
    for (vf, bars), (Mm, Pe, Me14) in DLR.TAB2.items():
        phi_f, fc, E, fR1, fR3, _, _ = DLR.MAT[vf]
        n, phi, fy, fu = CARP.BARS[bars]
        fct = fctm_from_fck(fc - 8.0)
        law = law_envelope(law_mc2010_plain(fct, gf_mc2010(fc)), law_mc2010_frc_linear(fR1, fR3))
        t0, t1 = tcm_bond(fct)
        rb = [Rebar(CARP.C, n * np.pi * phi ** 2 / 4.0, phi, Steel(0.999 * fu, fu, CARP.ES, 0.10), t0, t1)]
        comp = compression_ec2(fc, CARP.EC)
        r = CrackedHinge(Section.rectangle(CARP.B, CARP.H), law, comp, CARP.EC, rebars=rb,
                         n_layers=300).response(theta_max=0.06, n=150, stop_on_failure=False)
        ok = np.cumprod(-r["eps_top"] <= comp.eps_cu).astype(bool)
        MC = r["M"][ok].max() / 1e6
        MB = carp[(vf, bars)]["Pmod"] * CARP.L / 4.0 / 1000.0 + Msw / 1e6
        Mexp = Pe * CARP.L / 4.0 / 1000.0  # their P_m,exp with the 1200 mm span
        rows.append(dict(vf=vf, bars=bars, MB=MB, MC=MC, Mpcm=Mm, Pe=Pe, Mexp12=Mexp, Mexp14=Me14,
                         Pdig=CARP.PEXP[vf][bars][0]))
    def st(key, ref):
        x = np.array([r[key] / r[ref] for r in rows])
        return _stats(x)
    for key, nm in (("MB", "B"), ("MC", "C"), ("Mpcm", "P")):
        mu, cov = st(key, "Mexp12")
        put(f"dlr{nm}mean", mu)
        put(f"dlr{nm}cov", cov, "{:.1f}")
    mu, cov = st("Mpcm", "Mexp14")
    put("dlrPmeanFourteen", mu)
    x = np.array([r["Pdig"] / r["Pe"] for r in rows])
    mu, cov = _stats(x)
    put("dlrDigmean", mu)
    put("dlrDigmax", 100 * np.max(np.abs(x - 1)), "{:.0f}")
    with open(os.path.join(TAB, "tab_validation_delarosa.tex"), "w") as f:
        f.write("\\begin{tabular}{r l r r r r r}\n\\toprule\n")
        f.write("$V_f$ & bars & $M^{exp}$ & FMSD B & FMSD C & PCM-RC & $M^{exp}$\\\\\n")
        f.write("{[kg/m$^3$]} & & ($L=1.2$ m) & & & (2025) & ($l=1.4$ m)\\\\\n\\midrule\n")
        for r in rows:
            f.write(f"{r['vf']} & {r['bars'].replace('Ø', '$\\phi$')} & {r['Mexp12']:.1f} & {r['MB']:.1f} & {r['MC']:.1f} & "
                    f"{r['Mpcm']:.1f} & {r['Mexp14']:.1f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(4.6, 3.4))
        x = [r["Mexp12"] for r in rows]
        for key, lab_en, lab_es, j, mk in (("MB", "FMSD, inputs B", "FMSD, datos B", 0, "o"),
                                          ("MC", "FMSD, inputs C ($f_{R}$, MC2010)", "FMSD, datos C ($f_{R}$, MC2010)", 1, "s"),
                                          ("Mpcm", "planar crack R/SFRC (2025)", "fisura plana HA+HRFA (2025)", 2, "^")):
            ax.plot(x, [r[key] for r in rows], ls="none", marker=mk, color=P[j], ms=5, label=T(lang, lab_en, lab_es))
        ax.plot([4, 32], [4, 32], color="#52514e", lw=0.8)
        ax.set_xlim(4, 32)
        ax.set_ylim(4, 32)
        ax.set_xlabel(T(lang, "experimental $M_{max}$ ($L=1200$ mm) [kNm]", "$M_{max}$ experimental ($L=1200$ mm) [kNm]"))
        ax.set_ylabel(T(lang, "model $M_{max}$ [kNm]", "$M_{max}$ del modelo [kNm]"))
        ax.legend(fontsize=7)
        fig.tight_layout()
        save(fig, "fig_validation_delarosa", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_delarosa.json"), "w"), indent=1, default=float)


def sensitivity_s():
    """Sensitivity of key outputs to the hinge width s (s/h = 0.25, 0.5, 1.0)."""
    rows = []
    plain = law_mc2010_plain(FCTM, GF)
    for ratio in (0.25, 0.5, 1.0):
        out = {"r": ratio}
        for h in (200, 800):
            hg = CrackedHinge(Section.rectangle(100, h), plain, compression_ec2(FCM, ECM), ECM, s=ratio * h, n_layers=400)
            rr = hg.response(theta_max=40 * FCTM / ECM * 2 * ratio, n=140)
            out[f"mor{h}"] = rr["M"].max() / (100 * h ** 2 / 6) / FCTM
        hsp = 125.0
        law = law_annexL(FR1K, FR3K, FCTM, ECM, lcs=125.0, wu=4.5)
        hg = CrackedHinge(Section.rectangle(150, hsp), law, compression_elastic_plastic(FCM, ECM), ECM, s=ratio * hsp, n_layers=300)
        beam = SimplyBeam(hg, 500.0, "3pb", a0=25.0, h_full=150.0, b_full=150.0)
        rr = beam.response(theta_max=2 * 4.0 / hsp, n=140, stop_on_failure=False)
        for j, c in ((1, 0.5), (3, 2.5)):
            out[f"fR{j}"] = 3 * float(np.interp(c, rr["CMOD"], rr["P"])) * 500 / (2 * 150 * hsp ** 2)
        hgb = hinge_rc(300, 500, 603.0, 16, FR1K, FR3K)
        hgb.s = ratio * 500.0
        rb = hgb.response(theta_max=0.04, n=180)
        out["Mhyb"] = rb["M"].max() / 1e6
        rows.append(out)
    with open(os.path.join(TAB, "tab_sensitivity_s.tex"), "w") as f:
        f.write("\\begin{tabular}{c cc cc c}\n\\toprule\n$s/h$ & $f_{ct,fl}/f_{ct}$ ($h$=200) & $f_{ct,fl}/f_{ct}$ ($h$=800) & $f_{R,1}$ & $f_{R,3}$ & $M_{max}$ hybrid [kNm]\\\\\n\\midrule\n")
        for r in rows:
            f.write(f"{r['r']:.2f} & {r['mor200']:.3f} & {r['mor800']:.3f} & {r['fR1']:.2f} & {r['fR3']:.2f} & {r['Mhyb']:.1f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    put("sensMhybVar", 100 * (max(r["Mhyb"] for r in rows) / min(r["Mhyb"] for r in rows) - 1), "{:.1f}")
    put("sensfRthreeVar", 100 * (max(r["fR3"] for r in rows) / min(r["fR3"] for r in rows) - 1), "{:.1f}")
    put("sensfRoneVar", 100 * (max(r["fR1"] for r in rows) / min(r["fR1"] for r in rows) - 1), "{:.1f}")
    put("sensMORVar", 100 * (max(r["mor200"] for r in rows) / min(r["mor200"] for r in rows) - 1), "{:.1f}")


# ------------------------------------------------------------------------------------------
# 1.1: member checks (shear, N-M, columns, punching) -- mechanical models vs code reference
# ------------------------------------------------------------------------------------------

def validation_masm_cuenca():
    """MASM (Mari et al. 2020, Eq. 26) on the Cuenca et al. (2018) beams failing in shear, with the
    residual stress sigma_tu taken from the cohesive laws used by FracHinge (no fitting)."""
    rows = []
    for (ref, bid, fib, vf, bw, d, ad, rho, fc, fR1, fR3, fm, Vu) in CUENCA_DB:
        if fm not in ("S", "SY"):
            continue
        out = {}
        for key, st in (("rp", fR3 / 3.0), ("al", 0.33 * fR3)):
            r = members.shear_masm_sfrc(bw, d, ad * d, rho / 100.0, fc, st)
            out[key] = Vu * 1e3 / r["V_u"]
        rows.append(dict(ref=ref, d=d, fR3=fR3, **out))
    for key, nm in (("rp", "RP"), ("al", "AL")):
        mu, cov = _stats([r[key] for r in rows])
        put(f"masm{nm}mean", mu)
        put(f"masm{nm}cov", cov, "{:.1f}")
        put(f"masm{nm}unsafe", sum(r[key] < 1.0 for r in rows), "{:d}")
        big = [r[key] for r in rows if r["d"] >= 900]
        small = [r[key] for r in rows if r["d"] < 900]
        put(f"masm{nm}big", np.mean(big))
        put(f"masm{nm}small", np.mean(small))
    put("masmN", len(rows), "{:d}")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(figsize=(4.6, 3.2))
        x = [r["d"] for r in rows]
        y = [r["al"] for r in rows]
        c = [P[0] if r["fR3"] < 3 else (P[1] if r["fR3"] < 5 else P[2]) for r in rows]
        ax.scatter(x, y, c=c, s=22, edgecolors="white", linewidths=0.6, zorder=3)
        ax.axhline(1.0, color="#52514e", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("$d$ [mm]")
        ax.set_ylabel("$V_{u,exp}/V_{MASM}$ (" + T(lang, "mean values", "valores medios") + ")")
        for j, lb in enumerate(("$f_{R3}<3$", "$3\\leq f_{R3}<5$", "$f_{R3}\\geq5$ MPa")):
            ax.scatter([], [], c=P[j], s=22, label=lb)
        ax.legend(fontsize=7)
        fig.tight_layout()
        save(fig, "fig_validation_masm", lang)
    json.dump(rows, open(os.path.join(TAB, "validation_masm.json"), "w"), indent=1, default=float)


def example_members():
    """Worked examples of the member checks: N-M interaction, slender column, shear with stirrups
    and fibres, punching. Mechanical / fracture-mechanics result first, code reference second."""
    fck, fcm, fR1k, fR3k = 30.0, 38.0, 4.0, 3.5
    E = ecm_from_fcm(fcm)
    law = law_envelope(law_linear_softening(fctm_from_fck(fck), gf_mc2010(fcm)),
                       law_mc2010_frc_linear(fR1k, fR3k))
    sec = Section.rectangle(300.0, 300.0)
    comp = compression_annexL(fcm, fR1k)
    st = Steel(500.0, 540.0)
    rb = [Rebar(40.0, 452.0, 12.0, st), Rebar(260.0, 452.0, 12.0, st)]
    Ns = np.linspace(-0.95 * fcm * sec.area, 0.0, 14)
    ia = members.interaction_hinge(sec, law, comp, E, rebars=rb, N_values=Ns, n_layers=150, n=60)
    layers = [(300.0 - 40.0, 452.0, 500.0), (40.0, 452.0, 500.0)]
    ra = codes.interaction_annexL(300.0, 300.0, fcm, 0.33 * fR3k, layers, gamma_c=1.0,
                                  N_values=Ns)
    k0 = int(np.argmin(np.abs(Ns)))
    put("exMNzeroHinge", ia["M"][k0] / 1e6, "{:.0f}")
    put("exMNmaxHinge", np.nanmax(ia["M"]) / 1e6, "{:.0f}")
    put("exMNmaxAL", np.nanmax(ra["M"]) / 1e6, "{:.0f}")
    cols = []
    for L in (3000.0, 6000.0, 9000.0):
        cb = members.column_buckling_hinge(sec, law, comp, E, L, 20.0, rebars=rb, n_layers=120, n=50,
                                           N_grid=np.linspace(0.05, 1.0, 20) * fcm * sec.area)
        cols.append((L, cb["N_u"], cb["N_euler"]))
        put(f"exColNu{int(L / 1000)}", cb["N_u"] / 1e3, "{:.0f}")
        put(f"exColNe{int(L / 1000)}", cb["N_euler"] / 1e3, "{:.0f}")
    sh_rc = members.shear_masm_design(300.0, 450.0, 1350.0, 0.012, fck, 0.0)
    sh_f = members.shear_masm_design(300.0, 450.0, 1350.0, 0.012, fck, 0.33 * fR3k / 1.5)
    sh_fs = members.shear_masm_design(300.0, 450.0, 1350.0, 0.012, fck, 0.33 * fR3k / 1.5, rho_w=0.001)
    put("exShRC", sh_rc["V_Rd"] / 1e3, "{:.0f}")
    put("exShF", sh_f["V_Rd"] / 1e3, "{:.0f}")
    put("exShFS", sh_fs["V_Rd"] / 1e3, "{:.0f}")
    fctk = 0.7 * fctm_from_fck(fck)
    ref = codes.VRd_frc_with_stirrups(300.0, 450.0, 0.012, fck, 0.33 * fR3k, fctk, 0.001 * 300.0,
                                      435.0, cot_theta=2.5)
    put("exShRefF", ref["V_Rd,F"] / 1e3, "{:.0f}")
    put("exShRefFS", ref["V_Rd"] / 1e3, "{:.0f}")
    lawf = law_mc2010_frc_linear(fR1k, fR3k)
    p0 = members.punching_csct_sfrc(rc=150.0, d=200.0, rs=1100.0, fc=fcm, fy=550.0, mR=150000.0)
    p1 = members.punching_csct_sfrc(rc=150.0, d=200.0, rs=1100.0, fc=fcm, fy=550.0, mR=150000.0,
                                    fibre_law=lawf)
    put("exPunRC", p0["V_R"] / 1e3, "{:.0f}")
    put("exPunF", p1["V_R"] / 1e3, "{:.0f}")
    put("exPunFf", p1["V_Rf"] / 1e3, "{:.0f}")
    put("exPunW", p1.get("w_crack", np.nan), "{:.2f}")
    for lang in ("en", "es"):
        fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2))
        ax[0].plot(ia["M"] / 1e6, -ia["N"] / 1e3, color=P[0], lw=2,
                   label=T(lang, "cohesive hinge (mean)", "rótula cohesiva (media)"))
        ax[0].plot(ra["M"] / 1e6, -ra["N"] / 1e3, color=P[1], lw=2, ls="--",
                   label=T(lang, "Annex L block (ref.)", "bloque Anejo L (ref.)"))
        ax[0].set_xlabel("$M$ [kN m]")
        ax[0].set_ylabel(T(lang, "$N$ (compression) [kN]", "$N$ (compresión) [kN]"))
        ax[0].legend(fontsize=7)
        ax[0].set_title(T(lang, "(a) N-M interaction 300x300", "(a) Interacción N-M 300x300"))
        ax[1].bar([f"{int(L / 1000)} m" for L, _, _ in cols], [c[1] / 1e3 for c in cols], color=P[0], width=0.5)
        ax[1].set_ylabel("$N_u$ [kN]")
        ax[1].set_title(T(lang, "(b) slender column, $e_0$ = 20 mm", "(b) pilar esbelto, $e_0$ = 20 mm"))
        fig.tight_layout()
        save(fig, "fig_members", lang)


def write_numbers():
    import re
    for fname, dec in (("numbers.tex", "."), ("numbers_es.tex", "{,}")):
        with open(os.path.join(TAB, fname), "w") as f:
            f.write("% Auto-generated by examples/generate_all.py -- do not edit\n")
            f.write("\\makeatletter\n")
            for k, v in NUM.items():
                v = str(v).replace("%", "\\%")
                if dec != ".":
                    v = re.sub(r"(?<=\d)\.(?=\d)", dec, v)
                f.write(f"\\expandafter\\def\\csname fh@{k}\\endcsname{{{v}}}\n")
            f.write("\\makeatother\n")
            f.write("\\providecommand{\\val}[1]{\\csname fh@#1\\endcsname}\n")
    json.dump(NUM, open(os.path.join(TAB, "numbers.json"), "w"), indent=1)


if __name__ == "__main__":
    t0 = time.time()
    steps = [example1, fig_laws, fig_compression, example2, fig_profiles, validation_size_effect,
             example5, example6, example7, fig_lefm, example4, validation_gao, validation_bosco, validation_carpinteri, validation_delarosa2025, validation_minelli,
             validation_cuenca_shear, validation_tiberti, sensitivity_s, compare_pcm, derived_numbers,
             validation_masm_cuenca, example_members]
    if len(sys.argv) > 1:  # run selected steps only, keeping the other numbers already written
        import re as _re
        fn = os.path.join(TAB, "numbers.tex")
        if os.path.exists(fn):
            for k, v in _re.findall(r"csname fh@(.+?)\\endcsname\{(.*)\}$", open(fn).read(), _re.M):
                NUM[k] = v.replace("\\%", "%")
        steps = [st for st in steps if st.__name__ in sys.argv[1:]]
    for s in steps:
        print("->", s.__name__, f"{time.time() - t0:.0f}s")
        s()
        write_numbers()
    print("done", f"{time.time() - t0:.0f}s")
