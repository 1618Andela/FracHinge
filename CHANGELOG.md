# Changelog

## 1.5.0 (2026-11-16) – fracture-mechanics-only documentation
- Manual (ES/EN) restricted to a structural design approach based exclusively on fracture mechanics: code, empirical and plasticity checks removed (ACI 544, code shear and punching formulae, Annex L rigid-plastic block, EC2 nominal curvature, code crack-width and minimum-reinforcement formulae, code rules for plastic analysis and reinforcement replacement); the SFRC punching CCCM of Marí et al. (2023) removed from the document (no fracture ingredient beyond σ_tu, not reproducible).
- New validation: maximum moment of 25 slender RC/SFRC columns under axial load (Caballero-Morrison et al. 2012, 2013) with the cohesive hinge, no fitting: mean 1.00, CoV 11.3 % (`examples/data_caballero_columns.py`, `validation_columns_caballero`).
- Example 7 rewritten with the fracture-based shear models; figures and minimum-reinforcement table without code curves; new conclusions, limitations and future work; references to Lantsoght (2019) and Moraes Neto et al. (2013) databases (no fR per test).
- The code library and the app keep their code-reference module/tab, outside the method of the document.

## 1.4.0 (2026-11-16) – SFRC punching CCCM
- New `members.punching_cccm_sfrc` (Marí, Spinella & Rodríguez 2023, Structures 52:854–867) with σ_tu from the cohesive law, switchable thin-slab reduction F_red = 0.005 h ≤ 1; `members.sigma_tu_fibre_factor` (micromechanical reference, Eq. 20) and `members.fctm_mc2010`.
- 94 SFRC punching tests of the paper in `examples/data_mari2023_punching.py`; with the tabulated data the implementation gives 0.93 (CoV 23 %) against the published 1.105 (CoV 11.9 %) – the table lacks fR, slab spans and fibre diameter/strength, so the published predictions cannot be recomputed exactly (documented).
- App punching tab shows both mechanical models; manual section; 1 new test.

## 1.3.0 (2026-11-16) – doctoral theses review
- CCCM (Marí et al. 2014): default fct = 0.30 min(fcm, 60)^(2/3), the convention of the authors' validation. Independent validation on the 393 RC beams of Sarkhosh (2014, Appendix D; ACI-DAfStb database): mean 1.05, CoV 14.1 %, no bias with depth (`examples/data_sarkhosh2014.py`, `validation_cccm_rc`).
- MASM: fibre-orientation option `orientation="aligned"` (default, η_or = sin θ_cr, safe-side bound) or `"isotropic"`; justified with the orientation factors measured by Cunha (2010) and the observations of Grünewald (2004), Özyurt (2006) and García Alberti (2015). On the 43 Cuenca beams: 1.23 (CoV 26 %, min 0.81, 10 < 1) versus 0.99 (CoV 29 %, 27 < 1) isotropic.
- Manual: CCCM validation, orientation section, appendix with the 16 doctoral theses consulted and their use; slides, course exercise 10.5, app option, MATLAB `aligned` argument, 2 new tests.

## 1.2.0 (2026-11-16) – size effect in fibre shear
- MASM (Marí et al. 2020): the fibre term is multiplied by Bažant's energetic size-effect factor λ_s = √(2/(1+d/d0)) ≤ 1, d0 = 254 mm (Bažant & Kim 1984; ACI 318-19), default `size_effect="bazant"`; `size_effect=None` gives the original Eq. (26). New `members.size_factor_bazant`.
- Validation on the 43 Cuenca et al. (2018) shear failures, no fitting: mean 0.91 → 0.99, d ≥ 900 mm 0.75 → 0.99 (CoV 31.1 → 29.0 %). Two other fracture-mechanics variants are reported and rejected: d0 ∝ l_ch,FRC (no size effect, 0.91) and σ(w = ε d) from the critical shear crack theory (0.82).
- Manual section on size effect, updated slides, course exercise, app warning, MATLAB option and test.

## 1.1.0 (2026-11-16) – member checks
- New module `frachinge.members` (fracture mechanics / mechanical models, primary results): N–M interaction with the cohesive hinge; slender-column single-crack second-order model; CCCM shear for RC with or without stirrups (Marí et al. 2014) and MASM for SFRC without stirrups (Marí et al. 2020) with σ_tu from the cohesive law, strut crushing and a design wrapper; punching by the CSCT with cohesive fibre bridging (Muttoni 2008; Maya Duque et al. 2012).
- New code references in `frachinge.codes` (kept separate): V_Rd,s, V_Rd,max, MC2010 V_Rd,F + V_Rd,s, Annex L stirrup replacement, squash load, EC2 nominal curvature, EC2-2004 and MC2010 Level II punching, Annex L N–M interaction.
- `DesignInput.N_Ed` (compression positive) for bending with axial force; the hinge root search is widened when N ≠ 0 and stops gracefully when equilibrium under N is lost.
- Streamlit app reorganised into *FM* tabs (primary) and a *Code ref.* tab; new Shear, Columns and Punching tabs and N–M diagram.
- Validation of the MASM on the Cuenca et al. (2018) database (unconservative, reported as such); MATLAB `fh_masm_sfrc.m`, `fh_punching_csct.m`; 8 new tests.

## 1.0.0 (2026-11-16) – first public release
- Cohesive cracked-hinge sectional model (band s = h/2) for plain, RC, SFRC and hybrid sections, with Tension Chord Model bond, EC2 / Annex L / parabola-rectangle / elastic-plastic compression laws and axial force.
- Constitutive laws: MC2010 plain-concrete bilinear softening, MC2010 linear and rigid-plastic FRC laws, EN 1992-1-1:2023 Annex L laws, envelope rule.
- EN 14651 forward simulation and inverse identification of the stress–crack opening law.
- Code checks (EN 1992-1-1:2023 Annex L and fib MC2010): ULS bending, shear without stirrups, minimum reinforcement, crack spacing and width, classification, ECOV and partial-factor formats (lognormal mean residual strengths, kappa_G at all levels).
- Validation scripts and transcribed experimental databases with provenance: MC2010 size effect, Gao et al. (2025), Tiberti et al. (2017), Bosco et al. (1990), Carpinteri et al. (2015), De La Rosa et al. (2025), Minelli et al. (2014), Cuenca et al. (2018); comparison with the planar crack model.
- Python package, MATLAB/Octave implementation, Streamlit web application (with deployment files) and tests.
