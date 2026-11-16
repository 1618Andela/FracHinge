# Changelog

## 1.0.0 (2026-11-16) – first public release
- Cohesive cracked-hinge sectional model (band s = h/2) for plain, RC, SFRC and hybrid sections, with Tension Chord Model bond, EC2 / Annex L / parabola-rectangle / elastic-plastic compression laws and axial force.
- Constitutive laws: MC2010 plain-concrete bilinear softening, MC2010 linear and rigid-plastic FRC laws, EN 1992-1-1:2023 Annex L laws, envelope rule.
- EN 14651 forward simulation and inverse identification of the stress–crack opening law.
- Code checks (EN 1992-1-1:2023 Annex L and fib MC2010): ULS bending, shear without stirrups, minimum reinforcement, crack spacing and width, classification, ECOV and partial-factor formats (lognormal mean residual strengths, kappa_G at all levels).
- Validation scripts and transcribed experimental databases with provenance: MC2010 size effect, Gao et al. (2025), Tiberti et al. (2017), Bosco et al. (1990), Carpinteri et al. (2015), De La Rosa et al. (2025), Minelli et al. (2014), Cuenca et al. (2018); comparison with the planar crack model.
- Python package, MATLAB/Octave implementation, Streamlit web application (with deployment files) and tests.
