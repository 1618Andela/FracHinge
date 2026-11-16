"""
FracHinge — graphical application (Streamlit).

Run locally:   streamlit run app/streamlit_app.py
Deploy:        Streamlit Community Cloud / any server with `pip install -r requirements.txt`.

Bilingual (ES/EN) interface for the fracture-mechanics-based design of RC, SFRC and hybrid
sections: material laws, cracked-hinge response, ULS/SLS checks against EN 1992-1-1:2023 Annex L
and fib MC2010, EN 14651 inverse analysis and size effect.
"""
import io
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from frachinge import codes, members  # noqa: E402
from frachinge.design import DesignInput, build_hinge, run_design  # noqa: E402
from dataclasses import replace as dc_replace  # noqa: E402
from frachinge.materials import (compression_elastic_plastic, ecm_from_fcm, fctm_from_fck,  # noqa: E402
                                 gf_mc2010, law_annexL, law_envelope, law_mc2010_frc_linear,
                                 law_mc2010_plain, characteristic_length, compression_ec2)
from frachinge.hinge import CrackedHinge, Section, en14651_prism  # noqa: E402
from frachinge.materials import (Rebar, Steel, compression_annexL, law_linear_softening,  # noqa: E402
                                 law_mc2010_rigid_plastic)
from frachinge.inverse import identify_fct, identify_law  # noqa: E402

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]

st.set_page_config(page_title="FracHinge", layout="wide")
lang = st.sidebar.radio("Idioma / Language", ["Español", "English"], horizontal=True)
ES = lang == "Español"


def t(es, en):
    return es if ES else en


st.title("FracHinge")
st.caption(t("Dimensionamiento de secciones de hormigón armado, con fibras de acero e híbridas "
             "mediante mecánica de la fractura (rótula cohesiva fisurada).",
             "Fracture-mechanics-based design of RC, steel-fibre and hybrid concrete sections "
             "(cohesive cracked hinge)."))

# ----------------------------------------------------------------------------- inputs
sb = st.sidebar
sb.header(t("Sección", "Section"))
b = sb.number_input("b [mm]", 50.0, 5000.0, 300.0, 10.0)
h = sb.number_input("h [mm]", 50.0, 3000.0, 500.0, 10.0)
sb.header(t("Hormigón", "Concrete"))
fck = sb.number_input("fck [MPa]", 12.0, 90.0, 30.0, 1.0)
sb.header(t("Fibras (EN 14651)", "Fibres (EN 14651)"))
use_f = sb.checkbox(t("Hormigón con fibras de acero", "Steel-fibre concrete"), True)
fR1k = sb.number_input("fR,1k [MPa]", 0.0, 15.0, 4.0 if use_f else 0.0, 0.1) if use_f else 0.0
fR3k = sb.number_input("fR,3k [MPa]", 0.0, 15.0, 3.6 if use_f else 0.0, 0.1) if use_f else 0.0
kO = sb.number_input("κO", 0.3, 1.7, 1.0, 0.05) if use_f else 1.0
law_model = sb.selectbox(t("Ley post-fisuración", "Post-cracking law"), ["annexL", "mc2010"]) if use_f else "annexL"
sb.header(t("Armadura", "Reinforcement"))
nbars = sb.number_input(t("Nº de barras", "Number of bars"), 0, 40, 3)
phi = sb.number_input("ϕ [mm]", 6.0, 40.0, 16.0, 1.0)
cover_y = sb.number_input(t("Altura del c.d.g. de la armadura desde la cara traccionada [mm]",
                            "Bar centroid height from tension face [mm]"), 10.0, 500.0, 40.0, 1.0)
fyk = sb.number_input("fyk [MPa]", 200.0, 700.0, 500.0, 10.0)
fuk = sb.number_input("fuk [MPa]", 200.0, 800.0, 540.0, 10.0)
sb.header(t("Acciones", "Actions"))
MEd = sb.number_input("M_Ed [kNm]", 0.0, 1e5, 100.0, 1.0)
NEd = sb.number_input(t("N_Ed [kN] (compresión +)", "N_Ed [kN] (compression +)"), -1e4, 1e5, 0.0, 10.0)
Mser = sb.number_input("M_ser [kNm]", 0.0, 1e5, 70.0, 1.0)
wlim = sb.number_input("w_lim [mm]", 0.05, 1.0, 0.3, 0.05)
VfR = sb.number_input(t("CoV supuesto de fR (ECOV)", "Assumed CoV of fR (ECOV)"), 0.05, 0.5, 0.20, 0.01)

As = nbars * np.pi * phi ** 2 / 4.0
inp = DesignInput(b=b, h=h, fck=fck, fR1k=fR1k, fR3k=fR3k, kO=kO, law_model=law_model,
                  bars=[(cover_y, As, phi)] if nbars > 0 else [], fyk=fyk, fuk=fuk,
                  M_Ed=MEd * 1e6, N_Ed=NEd * 1e3, M_ser=Mser * 1e6, w_lim=wlim, V_fR=VfR)

st.info(t("**Mecánica de fractura (MF)** = resultado principal de FracHinge. "
          "**Ref. normativa** = fórmulas de código mostradas solo como contraste.",
          "**Fracture mechanics (FM)** = FracHinge primary result. "
          "**Code ref.** = code formulas shown for comparison only."))
tabs = st.tabs([t("MF · Materiales", "FM · Materials"), t("MF · Flexión + axil (ELU/ELS)", "FM · Bending + axial (ULS/SLS)"),
                t("MF · Cortante", "FM · Shear"), t("MF · Pilares", "FM · Columns"),
                t("MF · Punzonamiento", "FM · Punching"),
                t("MF · Análisis inverso EN 14651", "FM · EN 14651 inverse analysis"),
                t("MF · Efecto tamaño", "FM · Size effect"), t("Validación experimental", "Experimental validation"),
                t("Ref. normativa", "Code ref."), t("Informe", "Report")])

fcm = fck + 8.0
E = ecm_from_fcm(fcm)
fctm = fctm_from_fck(fck)
GF = gf_mc2010(fcm)

with tabs[0]:
    c1, c2 = st.columns(2)
    with c1:
        st.subheader(t("Propiedades derivadas", "Derived properties"))
        props = {"fcm [MPa]": fcm, "Ecm [MPa]": E, "fctm [MPa]": fctm, "GF [N/m]": GF * 1000,
                 "lch [mm]": characteristic_length(E, GF, fctm), "h/lch": h / characteristic_length(E, GF, fctm)}
        if use_f:
            props[t("Clase MC2010", "MC2010 class")] = codes.classify_mc2010(fR1k, fR3k)
            props[t("Clase EC2 (Tabla L.2)", "EC2 class (Table L.2)")] = codes.classify_ec2(fR1k, fR3k)
            props["(L.1) fR1k/fctk,0.05 ≥ 0.5"] = codes.annexL_ductility_ok(fR1k, 0.7 * fctm)
        st.table(pd.DataFrame({"": [str(round(v, 3)) if isinstance(v, float) else str(v) for v in props.values()]},
                              index=list(props.keys())))
    with c2:
        st.subheader(t("Leyes cohesivas σ(w)", "Cohesive laws σ(w)"))
        w = np.linspace(0, 3.0, 600)
        fig, ax = plt.subplots(figsize=(5, 3.2))
        plain = law_mc2010_plain(fctm, GF)
        ax.plot(w, plain.sigma(w), color=PALETTE[0], label=t("matriz (MC2010)", "matrix (MC2010)"))
        if use_f:
            ax.plot(w, law_mc2010_frc_linear(fR1k, fR3k, w_end=3.0).sigma(w), "--", color=PALETTE[1], label="MC2010")
            ax.plot(w, law_annexL(fR1k, fR3k, fctm, E, wu=3.0).sigma(w), "-.", color=PALETTE[2], label=t("Anejo L", "Annex L"))
        ax.set_xlabel("w [mm]")
        ax.set_ylabel("σ [MPa]")
        ax.legend(frameon=False)
        ax.grid(alpha=0.3)
        st.pyplot(fig)

with tabs[1]:
    if st.button(t("Calcular", "Compute"), type="primary") or "res" not in st.session_state:
        with st.spinner(t("Calculando la rótula cohesiva...", "Running the cohesive hinge...")):
            st.session_state["res"] = run_design(inp, n=140)
    res = st.session_state["res"]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("M_Rd (hinge) [kNm]", f"{res['M_Rd_hinge'] / 1e6:.1f}")
    c2.metric(t("M_Rd Anejo L (ref.) [kNm]", "M_Rd Annex L (ref.) [kNm]"), f"{res['annexL_uls']['M_Rd'] / 1e6:.1f}")
    c3.metric("ECOV R_d [kNm]", f"{res['ecov']['R_d'] / 1e6:.1f}")
    c4.metric(t("Ductilidad M_post/M_cr", "Ductility M_post/M_cr"), f"{res['ductility_ratio']:.2f}")
    if "ULS_ok" in res:
        (st.success if res["ULS_ok"] else st.error)(t("ELU flexión: ", "ULS bending: ") +
                                                      ("OK" if res["ULS_ok"] else t("NO CUMPLE", "NOT OK")))
    if "w_ser" in res:
        ok = res.get("SLS_ok", False)
        (st.success if ok else st.warning)(f"w(M_ser) = {res['w_ser']:.3f} mm  (w_lim = {wlim} mm)")
    if res["brittle_snapback"]:
        st.warning(t("Respuesta frágil: el pico post-fisuración es menor que el de fisuración (snap-back).",
                     "Brittle response: post-cracking peak below cracking peak (snap-back)."))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for i, lev in enumerate(("mean", "char", "design")):
        r = res["responses"][lev]
        ax[0].plot(r["theta"] * 1000, r["M"] / 1e6, color=PALETTE[i], label=lev)
        ax[1].plot(r["M"] / 1e6, r["w_bottom"], color=PALETTE[i], label=lev)
    if MEd > 0:
        ax[0].axhline(MEd, color="k", lw=0.8, ls=":")
    ax[0].set_xlabel("θ [mrad]")
    ax[0].set_ylabel("M [kNm]")
    ax[1].set_xlabel("M [kNm]")
    ax[1].set_ylabel(t("w en la cara traccionada [mm]", "w at tension face [mm]"))
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(frameon=False)
    st.pyplot(fig)
    st.caption(t("Fallo (leyes de cálculo): ", "Failure (design laws): ") + str(res["failure_design"]))

    st.markdown("---")
    st.subheader(t("Diagrama de interacción N–M", "N–M interaction diagram"))
    st.caption(t("Rótula cohesiva con las leyes de cálculo (resultado principal) frente al bloque rígido-plástico "
                 "del Anejo L con N (referencia). Se usa la armadura de la barra lateral.",
                 "Cohesive hinge with design laws (primary result) versus the rigid-plastic Annex L block with N "
                 "(reference). The sidebar reinforcement is used."))
    if st.button(t("Calcular diagrama N–M", "Compute N–M diagram")):
        with st.spinner("N–M..."):
            hgd = build_hinge(inp, "design", res["kappa_G"])
            fcd = codes.eta_cc(fck) * fck / 1.5
            Ns = np.linspace(-0.95 * fcd * b * h, 0.0, 12)
            ia = members.interaction_hinge(hgd.section, hgd.law, hgd.comp, hgd.E, rebars=hgd.rebars,
                                           N_values=Ns, n_layers=150, n=60)
            fFtud = kO * res["kappa_G"] * 0.33 * fR3k / 1.5 if use_f else 0.0
            lay = [(h - y, A, fyk / 1.15) for (y, A, ph) in inp.bars]
            ra = codes.interaction_annexL(b, h, fck, fFtud, lay, gamma_c=1.5, N_values=Ns)
        fig, ax = plt.subplots(figsize=(5.5, 3.6))
        ax.plot(ia["M"] / 1e6, -ia["N"] / 1e3, color=PALETTE[0], lw=2, label=t("MF: rótula cohesiva", "FM: cohesive hinge"))
        ax.plot(ra["M"] / 1e6, -ra["N"] / 1e3, "--", color=PALETTE[1], lw=2, label=t("Ref.: Anejo L", "Ref.: Annex L"))
        if MEd > 0:
            ax.plot([MEd], [NEd], "k*", ms=10, label="(M_Ed, N_Ed)")
        ax.set_xlabel("M_Rd [kNm]")
        ax.set_ylabel(t("N (compresión) [kN]", "N (compression) [kN]"))
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)
        st.pyplot(fig)
        st.dataframe(pd.DataFrame({"N [kN]": -ia["N"] / 1e3, "M_Rd FM [kNm]": ia["M"] / 1e6,
                                   "M_Rd ref [kNm]": ra["M"] / 1e6}).round(1), width="stretch")

with tabs[8]:
    st.subheader(t("Comprobaciones normativas (solo referencia)", "Code checks (reference only)"))
    d = h - cover_y
    rows = []
    if use_f:
        kG = res["kappa_G"]
        fFtuk = kO * kG * 0.33 * fR3k
        tau, eta = codes.tau_RdcF_annexL(As / (b * d) if As > 0 else 1e-4, fck, 0.9 * d, fFtuk / 1.5, fFtuk)
        rows.append([t("τRd,cF Anejo L (L.19) [MPa]", "τRd,cF Annex L (L.19) [MPa]"), f"{tau:.3f}"])
        rows.append(["ηcF (L.20)", f"{eta:.3f}"])
        rows.append([t("Sustitución cercos (L.32): fFtu,ef ≥ 0.08√fck", "Stirrup replacement (L.32): fFtu,ef ≥ 0.08√fck"),
                     str(fFtuk >= 0.08 * np.sqrt(fck))])
    rows.append([t("τRd,c sin fibras (8.16) [MPa]", "τRd,c without fibres (8.16) [MPa]"),
                 f"{codes.tau_Rdc(max(As / (b * d), 1e-4), fck, 0.9 * d):.3f}"])
    Asmin, Mcr = codes.as_min_bending(b, h, d, fctm, fyk, fck, kO * 0.33 * fR3k if use_f else 0.0)
    rows.append([t("As,min (12.1)/(L.29) [mm²]", "As,min (12.1)/(L.29) [mm²]"), f"{Asmin:.0f}"])
    st.table(pd.DataFrame(rows, columns=[t("Magnitud", "Quantity"), t("Valor", "Value")]))
    st.info(t("Recuerde: L.12.2.1 no permite sustituir por fibras la armadura longitudinal mínima en vigas; "
              "en losas, L.12.3.1 limita la sustitución a kAS·As,min (kAS = 0,5).",
              "Note: L.12.2.1 forbids replacing the minimum longitudinal reinforcement of beams by fibres; "
              "in slabs L.12.3.1 limits it to kAS·As,min (kAS = 0.5)."))

with tabs[5]:
    st.subheader(t("Identificación de σ(w) a partir de fR,j", "Identification of σ(w) from fR,j"))
    c = st.columns(6)
    FL = c[0].number_input("F_L [kN]", 1.0, 60.0, 14.6)
    Eexp = c[1].number_input("E [GPa]", 10.0, 60.0, 24.7)
    f1 = c[2].number_input("fR1", 0.0, 15.0, 3.651)
    f2 = c[3].number_input("fR2", 0.0, 15.0, 4.587)
    f3 = c[4].number_input("fR3", 0.0, 15.0, 4.817)
    f4 = c[5].number_input("fR4", 0.0, 15.0, 4.655)
    fcm_t = st.number_input("fcm [MPa]", 10.0, 120.0, 38.6)
    if st.button(t("Identificar", "Identify")):
        Ee = Eexp * 1000
        comp = compression_elastic_plastic(fcm_t, Ee)
        gfm = gf_mc2010(fcm_t)
        fct = identify_fct(FL * 1000, Ee, gfm, comp)
        idr = identify_law({1: f1, 2: f2, 3: f3, 4: f4}, fct, gfm, Ee, comp)
        st.write(f"f_ct = {fct:.2f} MPa; " + t("nodos", "nodes") + ": " +
                 ", ".join(f"σ({a:.2f})={s:.2f}" for a, s in zip(idr["nodes_w"], idr["nodes_s"])))
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.3))
        ax[0].plot(idr["sim"]["CMOD"], idr["sim"]["F"] / 1000, color=PALETTE[0])
        hsp = 125.0
        ax[0].plot([0.5, 1.5, 2.5, 3.5], [f * 2 * 150 * hsp ** 2 / (3 * 500) / 1000 for f in (f1, f2, f3, f4)], "ks")
        ax[0].set_xlabel("CMOD [mm]")
        ax[0].set_ylabel("F [kN]")
        ww = np.linspace(0, 4, 400)
        ax[1].plot(ww, idr["law"].sigma(ww), color=PALETTE[1])
        ax[1].set_xlabel("w [mm]")
        ax[1].set_ylabel("σ [MPa]")
        for a in ax:
            a.grid(alpha=0.3)
        st.pyplot(fig)
        csv = pd.DataFrame({"w_mm": idr["law"].w, "sigma_MPa": idr["law"].s}).to_csv(index=False)
        st.download_button(t("Descargar ley σ(w) (CSV)", "Download σ(w) law (CSV)"), csv, "sigma_w.csv")

with tabs[6]:
    st.subheader(t("Resistencia a flexotracción y tamaño", "Flexural strength and size"))
    hs = st.multiselect("h [mm]", [50, 100, 200, 400, 800, 1600], [100, 200, 400, 800])
    if st.button(t("Calcular efecto tamaño", "Compute size effect")):
        law = law_mc2010_plain(fctm, GF)
        out = []
        for hh in hs:
            hg = CrackedHinge(Section.rectangle(100, hh), law, compression_ec2(fcm, E), E, n_layers=400)
            r = hg.response(theta_max=40 * fctm / E, n=120)
            out.append([hh, hh / characteristic_length(E, GF, fctm), r["M"].max() / (100 * hh ** 2 / 6) / fctm,
                        codes.fct_flexural_mc2010(fctm, hh) / fctm])
        st.table(pd.DataFrame(out, columns=["h", "h/lch", "fct,fl/fct (hinge)", "MC2010 (5.1-8)"]))

with tabs[9]:
    res = st.session_state.get("res")
    if res:
        rep = {k: v for k, v in res.items() if k not in ("responses",)}
        txt = json.dumps(rep, indent=1, default=lambda o: float(o) if np.isscalar(o) else str(o))
        st.code(txt[:4000], language="json")
        st.download_button(t("Descargar informe (JSON)", "Download report (JSON)"), txt, "frachinge_report.json")
        lines = [r"\begin{tabular}{lr}", r"\toprule", r"Quantity & Value\\", r"\midrule",
                 f"$M_{{Rd}}$ hinge [kNm] & {res['M_Rd_hinge'] / 1e6:.2f}\\\\",
                 f"$M_{{Rd}}$ Annex L [kNm] & {res['annexL_uls']['M_Rd'] / 1e6:.2f}\\\\",
                 f"ECOV $R_d$ [kNm] & {res['ecov']['R_d'] / 1e6:.2f}\\\\",
                 f"$\\gamma_R$ (ECOV) & {res['ecov']['gamma_R']:.3f}\\\\",
                 f"$M_{{post}}/M_{{cr}}$ & {res['ductility_ratio']:.2f}\\\\",
                 r"\bottomrule", r"\end{tabular}"]
        st.download_button(t("Descargar tabla LaTeX", "Download LaTeX table"), "\n".join(lines), "frachinge_table.tex")

with tabs[7]:
    st.subheader(t("Validación con ensayos publicados", "Validation against published tests"))
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(ROOT, "examples"))
    st.markdown(t(
        "Resultados precalculados por `examples/generate_all.py` (carpeta `docs/tables`). "
        "Cargas de Carpinteri et al. (2015) digitalizadas de sus figuras.",
        "Results precomputed by `examples/generate_all.py` (folder `docs/tables`). "
        "Loads of Carpinteri et al. (2015) digitised from their figures."))
    ds = st.selectbox(t("Base de datos", "Database"),
                      ["Carpinteri, Cadamuro & Ventura (2015) – R/SFRC", "Minelli et al. (2014) – flexure vs shear",
                       "Cuenca et al. (2018) – shear database", "Bosco, Carpinteri & Debernardi (1990) – RC",
                       "Gao et al. (2025) – EN 14651", "Tiberti et al. (2017) – 81 EN 14651 series"])
    fname = {"Carpinteri": "validation_carpinteri.json", "Minelli": "validation_minelli.json",
             "Cuenca": "validation_cuenca.json", "Bosco": "validation_bosco.json",
             "Gao": "validation_gao.json", "Tiberti": "validation_tiberti.json"}[ds.split(",")[0].split(" ")[0]]
    path = os.path.join(ROOT, "docs", "tables", fname)
    if os.path.exists(path):
        df = pd.DataFrame(json.load(open(path)))
        st.dataframe(df, width="stretch")
        st.download_button("CSV", df.to_csv(index=False), fname.replace(".json", ".csv"))
    else:
        st.info(t("Ejecute examples/generate_all.py", "Run examples/generate_all.py"))
    if ds.startswith("Carpinteri"):
        import data_carpinteri2015 as CARP
        from frachinge.materials import law_linear_softening, Rebar, Steel, tcm_bond
        c1, c2, c3 = st.columns(3)
        vf = c1.selectbox(t("Fibras [kg/m³]", "Fibres [kg/m³]"), [0, 10, 20, 40], index=1)
        bars = c2.selectbox(t("Armadura", "Bars"), list(CARP.BARS.keys()), index=2)
        inp_src = c3.radio(t("Resistencia a tracción", "Tensile strength"),
                       ["B: ft Table 3 (MC90)", "A: ft Table 8 (notched flexural)"])
        if st.button(t("Recalcular la viga", "Recompute beam")):
            m = CARP.MAT[vf]
            ft, GFv = (m["ft3"], m["GF3"]) if inp_src.startswith("B") else (m["ft8"], m["GF8"])
            n, phi, fy, fu = CARP.BARS[bars]
            rb = []
            if n:
                t0, t1 = tcm_bond(fctm_from_fck(m["fc"] - 8.0))
                rb = [Rebar(CARP.C, n * np.pi * phi ** 2 / 4, phi, Steel(0.999 * fu, fu, CARP.ES, 0.10), t0, t1)]
            comp = compression_ec2(m["fc"], CARP.EC)
            hg = CrackedHinge(Section.rectangle(CARP.B, CARP.H), law_linear_softening(ft, GFv), comp, CARP.EC,
                              rebars=rb, n_layers=300)
            r = hg.response(theta_max=0.06, n=150, stop_on_failure=False)
            ok = np.cumprod(-r["eps_top"] <= comp.eps_cu).astype(bool)
            Msw = 24e-6 * CARP.B * CARP.H * CARP.L ** 2 / 8.0
            P = 4.0 * (r["M"][ok] - Msw) / CARP.L / 1000.0
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.plot(r["theta"][ok] * 1000, P, color="#2a78d6")
            pe, pmin, pmax = CARP.PEXP[vf][bars]
            ax.axhspan(pmin, pmax, color="#eb6834", alpha=0.2, label=t("ensayo (intervalo)", "test (range)"))
            ax.axhline(pe, color="#eb6834", lw=1)
            ax.set_xlabel("θ [mrad]")
            ax.set_ylabel("P [kN]")
            ax.legend()
            ax.grid(alpha=0.3)
            st.pyplot(fig)
            st.metric(t("Carga máxima calculada", "Computed peak load"), f"{P.max():.1f} kN",
                      f"{100 * (P.max() / pe - 1):+.0f} % vs test")

# ----------------------------------------------------------------------------- shear
with tabs[2]:
    st.subheader(t("Cortante en vigas (con o sin cercos, con o sin fibras)", "Beam shear (with/without stirrups and fibres)"))
    c = st.columns(4)
    d_s = c[0].number_input("d [mm]", 50.0, 3000.0, float(h - cover_y), 5.0, key="ds")
    a_s = c[1].number_input(t("a (luz de cortante) [mm]", "a (shear span) [mm]"), 100.0, 20000.0, float(3 * (h - cover_y)), 10.0)
    rho_w = c[2].number_input(t("ρw cercos [%]", "ρw stirrups [%]"), 0.0, 2.0, 0.0, 0.01) / 100.0
    fywk = c[3].number_input("fywk [MPa]", 200.0, 700.0, 500.0, 10.0)
    VEd = st.number_input("V_Ed [kN]", 0.0, 1e5, 150.0, 5.0)
    rho_l = max(As / (b * d_s), 1e-4)
    kG = st.session_state["res"]["kappa_G"] if "res" in st.session_state else 1.0
    fFtuk = kO * kG * 0.33 * fR3k if use_f else 0.0
    fFtud = fFtuk / 1.5
    cfm, cref = st.columns(2)
    with cfm:
        st.markdown(t("#### MF · Modelo mecánico (resultado principal)", "#### FM · Mechanical model (primary result)"))
        orient = st.radio(t("Orientación de fibras respecto a la fisura de cortante", "Fibre orientation w.r.t. the shear crack"),
                          ["aligned", "isotropic"], horizontal=True,
                          help=t("aligned: fibres alineadas con el eje de la viga (cota inferior, por defecto); "
                                 "isotropic: orientación aleatoria 3D (cota superior)",
                                 "aligned: fibres aligned with the beam axis (lower bound, default); "
                                 "isotropic: 3D random orientation (upper bound)"))
        r = members.shear_masm_design(b, d_s, a_s, rho_l, fck, fFtud, rho_w=rho_w, fywd=fywk / 1.15,
                                      orientation=orient)
        st.write(r["model"])
        comp_rows = [[k, f"{r[k] / 1e3:.1f}"] for k in ("V_c", "V_w", "V_l", "V_s", "V_u_concrete_fibres") if k in r]
        comp_rows += [[t("Aplastamiento de bielas V_max", "Strut crushing V_max"), f"{r['V_max_strut'] / 1e3:.1f}"],
                      ["V_Rd [kN]", f"{r['V_Rd'] / 1e3:.1f}"], ["x/d", f"{r['x_d']:.3f}"]]
        st.table(pd.DataFrame(comp_rows, columns=[t("Componente", "Component"), "kN"]))
        (st.success if r["V_Rd"] >= VEd * 1e3 else st.error)(f"V_Ed / V_Rd = {VEd * 1e3 / r['V_Rd']:.2f}")
        if use_f:
            st.write(f"λs (Bažant / ACI 318-19) = {r.get('lambda_F', 1.0):.3f};  "
                     f"η_or = {r.get('eta_orientation', 1.0):.3f}")
            st.info(t("Validación (Cuenca et al. 2018, 43 vigas sin cercos, sin ajuste): con efecto tamaño de Bažant y "
                      "fibras alineadas, media Vexp/Vcalc = 1,23 (CoV 26 %, mínimo 0,81, 10 de 43 < 1); con orientación "
                      "isótropa, 0,99 (CoV 29 %, 27 de 43 < 1). Referencia normativa al lado.",
                      "Validation (Cuenca et al. 2018, 43 beams without stirrups, no fitting): with Bažant size effect "
                      "and aligned fibres, mean Vexp/Vcalc = 1.23 (CoV 26 %, minimum 0.81, 10 of 43 < 1); with "
                      "isotropic orientation, 0.99 (CoV 29 %, 27 of 43 < 1). Code reference alongside."))
        else:
            st.info(t("Validación del CCCM (393 vigas de HA sin cercos, base de Sarkhosh 2014 / ACI-DAfStb): "
                      "media 1,05, CoV 14,1 %, sin sesgo con el canto.",
                      "CCCM validation (393 RC beams without stirrups, Sarkhosh 2014 / ACI-DAfStb database): "
                      "mean 1.05, CoV 14.1 %, no bias with depth."))
        st.caption(t("Coeficientes parciales no calibrados por fiabilidad (Marí et al. 2020).",
                     "Partial factors not calibrated by reliability analysis (Marí et al. 2020)."))
    with cref:
        st.markdown(t("#### Ref. normativa (solo contraste)", "#### Code ref. (comparison only)"))
        fctk = 0.7 * fctm
        ref = codes.VRd_frc_with_stirrups(b, d_s, rho_l, fck, fFtuk, fctk, rho_w * b, fywk / 1.15)
        rows = [["MC2010 V_Rd,F (7.7-5) [kN]", f"{ref['V_Rd,F'] / 1e3:.1f}"],
                ["V_Rd,s (cot θ = 2.5) [kN]", f"{ref['V_Rd,s'] / 1e3:.1f}"],
                ["V_Rd,max EC2-2004 (6.9) [kN]", f"{ref['V_Rd,max'] / 1e3:.1f}"],
                ["MC2010 V_Rd = V_Rd,F + V_Rd,s [kN]", f"{ref['V_Rd'] / 1e3:.1f}"]]
        if use_f:
            tau, eta = codes.tau_RdcF_annexL(rho_l, fck, 0.9 * d_s, fFtud, fFtuk)
            rows.append([t("Anejo L τRd,cF·b·z [kN]", "Annex L τRd,cF·b·z [kN]"), f"{tau * b * 0.9 * d_s / 1e3:.1f}"])
            rows.append([t("Sustitución de cercos mínimos", "Minimum-stirrup replacement"),
                         str(codes.annexL_stirrup_replacement_ok(fFtuk, fck))])
        st.table(pd.DataFrame(rows, columns=[t("Magnitud", "Quantity"), t("Valor", "Value")]))

# ----------------------------------------------------------------------------- columns
with tabs[3]:
    st.subheader(t("Pilares: compresión centrada, N–M y esbeltez", "Columns: centred compression, N–M and slenderness"))
    c = st.columns(5)
    Lc = c[0].number_input(t("Longitud de pandeo L [mm]", "Buckling length L [mm]"), 500.0, 30000.0, 4000.0, 100.0)
    e0 = c[1].number_input("e0 [mm]", 0.0, 1000.0, 20.0, 1.0)
    nb_c = c[2].number_input(t("Barras por cara", "Bars per face"), 0, 20, 2)
    phi_c = c[3].number_input("ϕ [mm] ", 6.0, 40.0, 16.0, 1.0)
    cov_c = c[4].number_input(t("Recubr. mecánico [mm]", "Bar cover to centroid [mm]"), 10.0, 200.0, 40.0, 1.0)
    NEd_c = st.number_input("N_Ed [kN] ", 0.0, 1e5, 800.0, 10.0)
    Af = nb_c * np.pi * phi_c ** 2 / 4.0
    inp_c = dc_replace(inp, bars=[(cov_c, Af, phi_c), (h - cov_c, Af, phi_c)] if nb_c > 0 else [])
    level = st.radio(t("Nivel de las leyes", "Law level"), ["design", "mean"], horizontal=True,
                     help=t("design: leyes de cálculo (γ); mean: predicción de ensayo",
                            "design: design laws (γ); mean: test prediction"))
    if st.button(t("Calcular pilar", "Compute column")):
        hgc = build_hinge(inp_c, level, 1.0)
        fc_lev = hgc.comp.fc
        with st.spinner(t("Modelo de fisura única con segundo orden...", "Single-crack second-order model...")):
            cb = members.column_buckling_hinge(hgc.section, hgc.law, hgc.comp, hgc.E, Lc, e0, rebars=hgc.rebars,
                                               n_layers=120, n=50,
                                               N_grid=np.linspace(0.04, 1.0, 25) * fc_lev * b * h)
        cfm, cref = st.columns(2)
        with cfm:
            st.markdown(t("#### MF · Rótula cohesiva + 2º orden (principal)", "#### FM · Cohesive hinge + 2nd order (primary)"))
            st.metric("N_u [kN]", f"{cb['N_u'] / 1e3:.0f}")
            st.write(f"N_Euler = {cb['N_euler'] / 1e3:.0f} kN")
            (st.success if cb["N_u"] >= NEd_c * 1e3 else st.error)(f"N_Ed / N_u = {NEd_c * 1e3 / max(cb['N_u'], 1):.2f}")
            st.caption(t("Hipótesis: una sola fisura a media altura, tramos no fisurados elásticos (en la línea de "
                         "Ruiz Carmona et al. 2009). Modelo sin validar experimentalmente en FracHinge; no sustituye "
                         "el método general de EN 1992.",
                         "Assumptions: single crack at mid-height, elastic uncracked segments (in the line of "
                         "Ruiz Carmona et al. 2009). Not experimentally validated in FracHinge; not a substitute for "
                         "the EN 1992 general method."))
        with cref:
            st.markdown(t("#### Ref. normativa (solo contraste)", "#### Code ref. (comparison only)"))
            Ac = b * h
            NRd = codes.column_axial_capacity(Ac, 2 * Af, fck, fyk)
            nc = codes.nominal_curvature_moment(NEd_c * 1e3, NEd_c * 1e3 * e0, Lc, h - cov_c, h, Ac, 2 * Af, fck, fyk)
            try:
                fFtud = kO * 0.33 * fR3k / 1.5 if use_f else 0.0
                MRd = codes.uls_bending_annexL(b, h, fck, fFtud, [(h - cov_c, Af, fyk / 1.15), (cov_c, Af, fyk / 1.15)],
                                               1.5, N_Ed=-NEd_c * 1e3)["M_Rd"]
            except ValueError:
                MRd = float("nan")
            st.table(pd.DataFrame([[t("N_Rd compresión centrada [kN]", "N_Rd centred compression [kN]"), f"{NRd / 1e3:.0f}"],
                                   ["e2 (EC2-2004 5.8.8) [mm]", f"{nc['e2']:.1f}"],
                                   ["M_Ed,2nd [kNm]", f"{nc['M_Ed'] / 1e6:.1f}"],
                                   [t("M_Rd(N_Ed) Anejo L [kNm]", "M_Rd(N_Ed) Annex L [kNm]"), f"{MRd / 1e6:.1f}"]],
                                  columns=[t("Magnitud", "Quantity"), t("Valor", "Value")]))

# ----------------------------------------------------------------------------- punching
with tabs[4]:
    st.subheader(t("Punzonamiento de losas (pilar interior)", "Slab punching (interior column)"))
    c = st.columns(4)
    d_p = c[0].number_input("d [mm] ", 50.0, 1000.0, 200.0, 5.0)
    rc = c[1].number_input(t("Radio (o semilado) del pilar [mm]", "Column radius (or half side) [mm]"), 50.0, 1000.0, 150.0, 5.0)
    rs = c[2].number_input("r_s [mm]", 200.0, 10000.0, 1100.0, 10.0,
                           help=t("Distancia al punto de momento nulo (≈0,22 L)", "Distance to zero moment (≈0.22 L)"))
    rho_p = c[3].number_input(t("ρ flexión [%]", "ρ flexural [%]"), 0.1, 3.0, 0.8, 0.05) / 100.0
    c2_ = st.columns(3)
    dg = c2_[0].number_input("dg [mm]", 4.0, 32.0, 16.0, 1.0)
    circ = c2_[1].checkbox(t("Pilar circular", "Circular column"), False)
    VEdp = c2_[2].number_input("V_Ed [kN] ", 0.0, 1e5, 500.0, 10.0)
    fy_m = 1.1 * fyk
    mR = rho_p * fy_m * d_p ** 2 * (1.0 - rho_p * fy_m / (2.0 * fcm))
    cfm, cref = st.columns(2)
    lawf = None
    if use_f:
        lawf = law_mc2010_frc_linear(fR1k, fR3k) if law_model == "mc2010" else law_annexL(fR1k, fR3k, fctm, E, include_matrix=False)
    with cfm:
        st.markdown(t("#### MF · CSCT + puenteo cohesivo (principal)", "#### FM · CSCT + cohesive bridging (primary)"))
        pr = members.punching_csct_sfrc(rc, d_p, rs, fcm, fy_m, mR, dg=dg, fibre_law=lawf, circular=circ)
        st.table(pd.DataFrame([["ψ [mrad]", f"{1000 * pr['psi']:.2f}"],
                               ["w = ψd/6 [mm]", f"{pr['psi'] * d_p / 6:.3f}"],
                               ["V_R,c [kN]", f"{pr['V_Rc'] / 1e3:.0f}"], ["V_R,f [kN]", f"{pr['V_Rf'] / 1e3:.0f}"],
                               ["V_flex [kN]", f"{pr['V_flex'] / 1e3:.0f}"], ["V_R [kN]", f"{pr['V_R'] / 1e3:.0f}"],
                               [t("Modo", "Mode"), pr["governing"]]], columns=[t("Magnitud", "Quantity"), t("Valor", "Value")]))
        st.caption(t("Valores medios (γ = 1). Muttoni (2008) + Maya Duque et al. (2012), σ(w) de la ley cohesiva. "
                     "Para cálculo divida por un coeficiente global o use la referencia normativa.",
                     "Mean values (γ = 1). Muttoni (2008) + Maya Duque et al. (2012), σ(w) from the cohesive law. "
                     "For design apply a global factor or use the code reference."))
        (st.success if pr["V_R"] >= VEdp * 1e3 else st.error)(f"V_Ed / V_R = {VEdp * 1e3 / pr['V_R']:.2f}")
        st.markdown(t("##### MF · CCCM para losas de HRFA (Marí et al. 2023)", "##### FM · CCCM for SFRC slabs (Marí et al. 2023)"))
        h_p = st.number_input(t("Canto de la losa h [mm]", "Slab depth h [mm]"), 50.0, 1200.0, float(d_p + 30.0), 5.0)
        st_p = kO * 0.33 * fR3k if use_f else 0.0
        pm = members.punching_cccm_sfrc(rc if circ else 4 * rc / np.pi, d_p, h_p, fcm, rho_p, st_p)
        st.table(pd.DataFrame([[t("σtu = 0,33 fR3 [MPa]", "σtu = 0.33 fR3 [MPa]"), f"{st_p:.2f}"],
                               ["V_cu [kN]", f"{pm['V_cu'] / 1e3:.0f}"], ["F_red", f"{pm['F_red']:.2f}"],
                               ["V_cf [kN]", f"{pm['V_cf'] / 1e3:.0f}"], ["V_u [kN]", f"{pm['V_u'] / 1e3:.0f}"]],
                              columns=[t("Magnitud", "Quantity"), t("Valor", "Value")]))
        st.caption(t("Valores medios. Los autores obtienen 1,105 (CoV 11,9 %) en 94 ensayos con σtu micromecánico; "
                     "con σtu de la ley cohesiva el modelo no está validado (la base no publica fR).",
                     "Mean values. The authors report 1.105 (CoV 11.9 %) on 94 tests with a micromechanical σtu; "
                     "with σtu from the cohesive law the model is not validated (the database lacks fR)."))
    with cref:
        st.markdown(t("#### Ref. normativa (solo contraste)", "#### Code ref. (comparison only)"))
        if circ:
            b0, u1 = 2 * np.pi * (rc + d_p / 2), 2 * np.pi * (rc + 2 * d_p)
        else:
            b0, u1 = 8 * rc + np.pi * d_p, 8 * rc + 4 * np.pi * d_p
        mc = codes.punching_mc2010_level2(b0, d_p, rs, d_p, fck, fyk / 1.15, 1.0, dg)
        st.table(pd.DataFrame([["EC2-2004 (6.47) V_Rd,c [kN]", f"{codes.punching_ec2_2004(u1, d_p, rho_p, fck) / 1e3:.0f}"],
                               ["MC2010 Level II V_Rd,c (m_Ed = m_Rd) [kN]", f"{mc['V_Rd,c'] / 1e3:.0f}"]],
                              columns=[t("Magnitud", "Quantity"), t("Valor", "Value")]))


st.sidebar.markdown("---")
st.sidebar.warning(t("Herramienta de apoyo al cálculo: no sustituye la comprobación del proyectista conforme a EN 1992-1-1:2023.", "Design-support tool: it does not replace the designer's verification according to EN 1992-1-1:2023."))
st.sidebar.caption("FracHinge v1.5.0 · MIT · " + t("Ver docs/ para manual, artículo y curso.",
                                                  "See docs/ for manual, article and course."))
