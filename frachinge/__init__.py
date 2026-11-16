"""FracHinge: fracture-mechanics-based design of RC, SFRC and hybrid (R/SFRC) sections.

Cohesive cracked-hinge model (Hillerborg fictitious crack + Olesen hinge + Tension Chord bond)
with code cross-checks (fib MC2010, EN 1992-1-1:2023 Annex L, ACI 544.4R).
"""
__version__ = "1.2.0"
from . import materials, hinge, codes, fracture, inverse, safety, design  # noqa: F401
