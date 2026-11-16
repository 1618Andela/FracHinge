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

from frachinge import codes  # noqa: E402
from frachinge.design import DesignInput, run_design  # noqa: E402
from frachinge.materials import (compression_elastic_plastic, ecm_from_fcm, fctm_from_fck,  # noqa: E402
                                 gf_mc2010, law_annexL, law_envelope, law_mc2010_frc_linear,
                                 law_mc2010_plain, characteristic_length, compression_ec2)
from frachinge.hinge import CrackedHinge, Section, en14651_prism  # noqa: E402
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
Mser = sb.number_input("M_ser [kNm]", 0.0, 1e5, 70.0, 1.0)
wlim = sb.number_input("w_lim [mm]", 0.05, 1.0, 0.3, 0.05)
VfR = sb.number_input(t("CoV supuesto de fR (ECOV)", "Assumed CoV of fR (ECOV)"), 0.05, 0.5, 0.20, 0.01)

As = nbars * np.pi * phi ** 2 / 4.0
inp = DesignInput(b=b, h=h, fck=fck, fR1k=fR1k, fR3k=fR3k, kO=kO, law_model=law_model,
                  bars=[(cover_y, As, phi)] if nbars > 0 else [], fyk=fyk, fuk=fuk,
                  M_Ed=MEd * 1e6, M_ser=Mser * 1e6, w_lim=wlim, V_fR=VfR)

tabs = st.tabs([t("Materiales", "Materials"), t("Sección / ELU-ELS", "Section / ULS-SLS"),
                t("Normativa", "Codes"), t("Análisis inverso EN 14651", "EN 14651 inverse analysis"),
                t("Efecto tamaño", "Size effect"), t("Informe", "Report"),
                t("Validación experimental", "Experimental validation")])

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
    c2.metric(t("M_Rd Anejo L [kNm]", "M_Rd Annex L [kNm]"), f"{res['annexL_uls']['M_Rd'] / 1e6:.1f}")
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

with tabs[2]:
    st.subheader(t("Comprobaciones normativas de referencia", "Reference code checks"))
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

with tabs[3]:
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

with tabs[4]:
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

with tabs[5]:
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

with tabs[6]:
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
        inp = c3.radio(t("Resistencia a tracción", "Tensile strength"),
                       ["B: ft Table 3 (MC90)", "A: ft Table 8 (notched flexural)"])
        if st.button(t("Recalcular la viga", "Recompute beam")):
            m = CARP.MAT[vf]
            ft, GFv = (m["ft3"], m["GF3"]) if inp.startswith("B") else (m["ft8"], m["GF8"])
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

st.sidebar.markdown("---")
st.sidebar.warning(t("Herramienta de apoyo al cálculo: no sustituye la comprobación del proyectista conforme a EN 1992-1-1:2023.", "Design-support tool: it does not replace the designer's verification according to EN 1992-1-1:2023."))
st.sidebar.caption("FracHinge v1.0.0 · MIT · " + t("Ver docs/ para manual, artículo y curso.",
                                                  "See docs/ for manual, article and course."))
