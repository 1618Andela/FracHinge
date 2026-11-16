# FracHinge

**Cálculo de secciones de hormigón armado (HA), reforzado con fibras de acero (HRFA) e híbridas (HA+HRFA) mediante mecánica de la fractura.**

FracHinge implementa el método **FMSD** (cálculo seccional basado en mecánica de la fractura). Es una rótula cohesiva
fisurada (fisura ficticia de Hillerborg, rótula no lineal de Olesen y adherencia del *Tension Chord Model*) vinculada a
las leyes y formatos de seguridad del **Código Modelo 2010** y del **Anejo L de EN 1992-1-1:2023**, con contraste
normativo (Anejo L, MC2010, ACI 544.4R).

## Validación (reproducible con `python examples/generate_all.py`)
| Comprobación | Resultado |
|---|---|
| Flexotracción frente a MC2010 Ec. (5.1-8), h = 100–1200 mm | diferencia ≤ 8 %, sin ajuste |
| Consistencia con EN 14651, clase 4c | ley lineal MC2010 +19 % en f_R1; Anejo L +7 % |
| 14 HRFA (Gao et al., 2025), f_R2 y f_R4 ciegos | media 1,03, CoV 7,7 % |
| 81 series EN 14651, 528 vigas (Tiberti et al., 2017) | f_L ciego 1,09 (CoV 12,0 %); f_R2, f_R4 ciegos 1,02 (CoV 7,9 %); ley MC2010 f_R1 1,19, Anejo L 1,04 |
| 29 vigas de baja cuantía (Bosco et al., 1990) | P_cr 1,09, P_y 0,84; frágil/dúctil 22/29, 0 inseguros |
| 20 vigas HA+HRFA, 0–40 kg/m³ × ρ = 0–3,14 % (Carpinteri, Cadamuro y Ventura, 2015; cargas digitalizadas) | P_max 0,90 (CoV 16,0 %); ρ = 0,25–1,13 %: 0,98 (CoV 9,3 %) |
| 9 vigas de gran canto sin cercos (Minelli et al., 2014) | vigas que agotan a flexión M_u/M_max 1,02–1,07; roturas por cortante 0,32–0,84 |
| 43 vigas de HRFA sin cercos (base de datos de Cuenca et al., 2018), comprobación normativa a cortante | MC2010 1,07 (CoV 19,4 %), Anejo L 1,14 (CoV 21,6 %) |

## Uso rápido
```bash
pip install -r requirements.txt
pytest -q
streamlit run app/streamlit_app.py
python examples/generate_all.py
cd matlab && octave-cli --eval run_tests
```

## Documentos
Manual técnico, artículo, presentaciones y curso, en español e inglés, en `docs/` (fuentes LaTeX y PDF compilados).

## Limitaciones
Ancho de rótula s = h/2 (hipótesis); fisura dominante única; adherencia rígido-plástica; carga monótona; cortante con
fórmulas normativas; el CoV de las resistencias residuales en el ECOV es un dato del usuario; queda pendiente la
validación con vigas HA+HRFA. El programa es una herramienta de análisis y contraste y **no sustituye** a la
normativa vigente.

Licencia del código: MIT. Documentos: CC BY 4.0.
