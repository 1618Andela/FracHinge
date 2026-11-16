# Changelog

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
