# FracHinge

**Cálculo de secciones de hormigón armado (HA), reforzado con fibras de acero (HRFA) e híbridas (HA+HRFA), secciones y elementos, mediante mecánica de la fractura.**

FracHinge implementa el método **FMSD** (cálculo seccional basado en mecánica de la fractura). Es una rótula cohesiva
fisurada (fisura ficticia de Hillerborg, rótula no lineal de Olesen y adherencia del *Tension Chord Model*) vinculada a
las leyes y formatos de seguridad del **Código Modelo 2010** y del **Anejo L de EN 1992-1-1:2023**, con contraste
normativo (Anejo L, MC2010, ACI 544.4R).

**Comprobación de elementos (v1.1).** Cada comprobación da un resultado *principal* de mecánica de fractura o modelo
mecánico y, por separado, una *referencia normativa*: flexocompresión (diagrama N–M con la rótula cohesiva), pilares
esbeltos (fisura cohesiva única + segundo orden, ley de compresión del Anejo L de Ruiz et al. 2023), cortante (enfoque unificado de Frosch, Yu, Cusatis y Bažant 2017 con la ley energética de efecto tamaño de Bažant y puenteo cohesivo de las fibras, principal; CCCM de
Marí et al. 2014 con G_f y cercos; MASM de Marí et al. 2020 para HRFA con σ_tu de la ley cohesiva efecto tamaño de Bažant y cota de orientación de fibras en el término de fibras) y punzonamiento
(enfoque unificado bidireccional con efecto tamaño energético y puenteo cohesivo, y CCCM para losas de HRFA de Marí, Spinella y Rodríguez 2023).

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
| Mismas 43 vigas, MASM con σ_tu = 0,33 f_R3 de la ley cohesiva, Ec. (26) original (sin ajuste) | 0,91 (CoV 31,1 %), 0,75 con d ≥ 900 mm: el término de fibras no tiene efecto tamaño |
| Ídem con el factor de efecto tamaño energético de Bažant λ_s = √(2/(1+d/254)) (ACI 318-19) en el término de fibras, orientación isótropa (sin ajuste) | 0,99 (CoV 29,0 %), 0,99 con d ≥ 900 mm; 27/43 cocientes < 1 |
| Ídem, λ_s + fibras alineadas con el eje, η_or = sen θ_cr (por defecto desde 1.3, cota del lado de la seguridad; datos de orientación de Cunha 2010) | 1,23 (CoV 26,0 %), mín. 0,81, 10/43 cocientes < 1 |
| 25 pilares esbeltos de HA y HRFA con axil (Caballero-Morrison et al. 2012, 2013), momento máximo con la rótula cohesiva, sin ajuste | 1,00 (CoV 11,3 %) |
| 393 vigas de HA sin cercos (base de Sarkhosh 2014 / ACI-DAfStb), CCCM de Marí et al. 2014 | 1,05 (CoV 14,1 %), sin sesgo con el canto (EC2 1,14, CoV 15,0 %) |
| 393 vigas de HA, enfoque unificado de Bažant et al. (2017) con efecto tamaño energético | 1,50 (CoV 15,5 %, nominal), sin sesgo de tamaño (1,48 / 1,57 / 1,48), 2 de 393 < 1 |
| 43 vigas de HRFA, enfoque unificado + puenteo cohesivo (alineado) | 1,60 (CoV 22,9 %), ninguna < 1 |
| 94 losas de HRFA a punzonamiento (Marí et al. 2023), Bažant bidireccional + puenteo (σ_tu micromecánico) | 1,79 (CoV 12,3 %), ninguna < 1 |

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
Ancho de rótula s = h/2 (hipótesis); fisura dominante única; adherencia rígido-plástica; carga monótona; el CoV de las
resistencias residuales en el ECOV es un dato del usuario. Elementos: con σ_tu = 0,33 f_R3 el MASM aún sobrestima las vigas con f_R3 alto
(relación σ_tu–f_R no calibrada); el modelo de pilar (fisura única) aún no se ha validado con ensayos con excentricidad, y el punzonamiento con puenteo cohesivo
solo se ha comprobado con σ_tu micromecánico (la base disponible no da f_R); los coeficientes parciales de los modelos mecánicos de cortante no están calibrados. El programa es una herramienta de análisis y contraste y **no sustituye** a la
normativa vigente.

Licencia del código: MIT. Documentos: CC BY 4.0.
