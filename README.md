# FracHinge

**Fracture-mechanics-based design of reinforced (RC), steel-fibre-reinforced (SFRC) and hybrid (R/SFRC) concrete sections.**
*Cálculo de secciones de hormigón armado, con fibras de acero e híbridas mediante mecánica de la fractura* — see [README_ES.md](README_ES.md).

FracHinge implements the **FMSD** (Fracture-Mechanics-based Sectional Design) method. It is a cohesive
cracked-hinge model (Hillerborg's fictitious crack + Olesen's non-linear hinge + Tension-Chord-Model bond) linked to
the constitutive laws and safety formats of the **fib Model Code 2010** and **EN 1992-1-1:2023 Annex L**, with code
cross-checks (Annex L, MC2010, ACI 544.4R).

## Features
- Cohesive laws: MC2010 plain bilinear, MC2010 linear / rigid-plastic FRC, Annex L (w-form), envelope (MC2010 5.6.5), polygonal.
- Compression: EC2 (5.6)–(5.10), Annex L (L.12)–(L.14), parabola-rectangle, elastic-plastic.
- Cracked hinge with bars (TCM bond), axial force, T-sections and the c_f,dur layer (L.6).
- EN 14651 simulation and **inverse identification** of σ(w) from F_L and f_R,j.
- ULS/SLS design workflow: partial factors, **ECOV** and global factor (MC2010 7.11), ductility index, size-dependent minimum reinforcement, direct crack width.
- Code formulae: classification (MC2010, Table L.2), κ_G, shear (8.11), (8.16), (L.19)–(L.20), MC2010 (7.7-5), crack spacing (L.26)–(L.28), MC2010 crack width (7.7-22), A_s,min (12.1)/(L.29), Henager–Doherty (ACI 544.4R).
- Interfaces: Python library, **Streamlit app (ES/EN)**, MATLAB/Octave functions.

## Validation (reproducible with `python examples/generate_all.py`)
| Check | Result |
|---|---|
| Flexural strength vs MC2010 Eq. (5.1-8), h = 100–1200 mm | ≤ 8 % difference, no fitting |
| EN 14651 consistency, class 4c | MC2010 linear law +19 % on f_R1; Annex L +7 % |
| 14 SFRC mixes (Gao et al., 2025), blind f_R2, f_R4 | mean 1.03, CoV 7.7 % |
| 81 EN 14651 series, 528 beams (Tiberti et al., 2017) | blind f_L 1.09 (CoV 12.0 %); blind f_R2, f_R4 1.02 (CoV 7.9 %); MC2010 law f_R1 1.19, Annex L 1.04 |
| 29 lightly reinforced beams (Bosco et al., 1990) | P_cr 1.09, P_y 0.84; brittle/ductile 22/29, 0 unsafe |
| 20 R/SFRC beams, 0–40 kg/m³ × ρ = 0–3.14 % (Carpinteri, Cadamuro & Ventura, 2015; loads digitised) | P_max 0.90 (CoV 16.0 %); ρ = 0.25–1.13 %: 0.98 (CoV 9.3 %) |
| 9 deep beams without stirrups (Minelli et al., 2014) | flexure-governed beams M_u/M_max 1.02–1.07; shear failures 0.32–0.84 |
| 43 SFRC beams without stirrups (Cuenca et al., 2018 database), code shear checks | MC2010 1.07 (CoV 19.4 %), Annex L 1.14 (CoV 21.6 %) |

## Quick start
```bash
pip install -r requirements.txt
pytest -q                              # 18 tests
streamlit run app/streamlit_app.py     # graphical application (see DEPLOY.md to publish it online)
python examples/generate_all.py        # all figures, tables and numbers of the documents
cd matlab && octave-cli --eval run_tests
```
```python
from frachinge.design import DesignInput, run_design
res = run_design(DesignInput(b=300, h=500, fck=30, fR1k=4.0, fR3k=3.6,
                             bars=[(40, 603, 16)], M_Ed=120e6, M_ser=70e6))
print(res["M_Rd_hinge"]/1e6, res["annexL_uls"]["M_Rd"]/1e6, res["ecov"]["R_d"]/1e6)
```

## Repository layout
```
frachinge/   materials, hinge, inverse, codes, fracture, safety, design
app/         Streamlit application
matlab/      fh_*.m functions, examples, run_tests.m
examples/    generate_all.py, data_gao2025.py, data_bosco1990.py
tests/       pytest suite
docs/        manual (ES/EN), article (ES/EN), beamer (ES/EN), course (ES/EN), figures, tables, bib
```

## Documents
| | Spanish | English |
|---|---|---|
| Technical manual | `docs/manual/manual_es.pdf` | `docs/manual/manual_en.pdf` |
| Journal article (Elsevier preprint) | `docs/article/article_es.pdf` | `docs/article/article_en.pdf` |
| Conference slides (Beamer) | `docs/beamer/beamer_es.pdf` | `docs/beamer/beamer_en.pdf` |
| Course (PhD / industry) | `docs/course/curso_es.pdf` | `docs/course/course_en.pdf` |

## Limitations
Hinge width s = h/2 (assumption); single dominant crack; rigid-plastic bond; monotonic loading; shear via code
formulae; CoV of residual strengths in ECOV is a user input; validation with fibre-reinforced RC beams pending.
The software is an analysis/cross-check tool and **does not replace** the codes in force.

## Citation and licence
See `CITATION.cff`. Code: MIT. Documents: CC BY 4.0.
