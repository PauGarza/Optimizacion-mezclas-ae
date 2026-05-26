# Plan de Proyecto Final — MAT-34420 Métodos Numéricos y Optimización
**ITAM · Primavera 2026 · Prof. J. Ezequiel Soto**

---

## 1. Identificación del proyecto

| Campo | Detalle |
|---|---|
| **Título** | Formulación óptima de mezclas de aceites esenciales mediante Programación Lineal |
| **Materia** | MAT-34420 — Métodos Numéricos y Optimización |
| **Integrante** | Paulina Garza |
| **Entregables** | Notebook `.ipynb` + Presentación `.pdf` (10–12 min) |
| **Fecha de entrega** | Mayo 2026 |

---

## 2. Planteamiento del problema

### 2.1 Contexto

La aromaterapia y la naturopatía utilizan mezclas de aceites esenciales con fines terapéuticos específicos: relajación, acción antimicrobiana, estimulación cognitiva, entre otros. Cada aceite esencial tiene una **composición química característica**, determinada por cromatografía de gases acoplada a espectrometría de masas (GC-MS), que consiste en un vector de porcentajes de compuestos orgánicos volátiles (linalool, eucaliptol, mentol, limoneno, etc.).

El problema práctico es: dado un objetivo terapéutico, **¿en qué proporciones mezclar un conjunto de aceites esenciales para maximizar la presencia de los compuestos activos deseados**, respetando restricciones de seguridad, costo y formulación?

### 2.2 Por qué es un problema de Programación Lineal

- La mezcla final es una **combinación convexa** de los vectores de composición de cada aceite.
- Maximizar una suma ponderada de concentraciones es una **función objetivo lineal** en las fracciones de mezcla.
- Las restricciones de seguridad, proporciones y presupuesto son todas **restricciones lineales**.
- El problema tiene solución exacta y admite análisis de sensibilidad vía **variables duales (KKT)**.

---

## 3. Técnicas del curso a aplicar

| Técnica | Sección del curso | Rol en el proyecto |
|---|---|---|
| **Programación Lineal — Método Símplex** | §4.1 | Método principal de solución |
| **Condiciones KKT y dualidad** | §4.4 | Interpretación de restricciones activas y variables duales |
| **Método de puntos interiores** | §4.5 | Método alternativo; comparación con Símplex |

---

## 4. Fuente de datos

### 4.1 Base de datos principal: sCentInDB

- **Fuente:** [sCentInDB](https://cb.imsc.res.in/scentindb/) — base de datos de composición química de aceites esenciales de plantas medicinales indias.
- **Publicación:** Samal et al. (2026). *Molecular Diversity*. DOI: 10.1007/s11030-025-11215-5
- **Archivos descargados:**

| Archivo | Contenido | Filas útiles |
|---|---|---|
| `essential_oil_scentindb.csv` | Planta, muestra ESO, compuesto, % GC-MS | 85,341 |
| `therapy_scentindb.csv` | Usos terapéuticos con códigos UMLS/MeSH | 515 |
| `chemical_scentindb.csv` | Diccionario de compuestos (SMILES, PubChem ID) | 3,420 |
| `uses_scentindb.csv` | Usos industriales/comerciales por planta | — |

### 4.2 Aceites de interés confirmados en la BD

| Aceite | Especie | Muestras ESO | Filas GC-MS | Compuesto dominante |
|---|---|:---:|:---:|---|
| Lavanda | *Lavandula angustifolia* | 9 | 414 | Linalyl acetate (27–52%), Linalool (17–31%) |
| Menta | *Mentha piperita* | 12 | 482 | Menthol |
| Romero | *Rosmarinus officinalis* | 12 | 360 | Camphor (24–36%), Eucalyptol (22–24%) |
| Eucalipto | *Eucalyptus radiata* | 1 | 24 | Eucalyptol (74%) |
| Eucalipto | *Eucalyptus globulus* | 2 | 38 | p-Cymene (32%), Eucalyptol (17%) |
| Árbol de té | *Melaleuca alternifolia* | 1 | 26 | 4-Terpineol (48%), Eucalyptol (5%) |
| Incienso† | *Boswellia serrata* | 5 | 168 | α-Thujene (22–61%), α-Pinene (7–11%) |
| Limón | *Citrus limon* | 3 | 83 | Limonene (38–92%), Citral |

> †*Boswellia sacra* y *B. carterii* (especies doTERRA) no están en sCentInDB; se usará *B. serrata* (incienso indio) con nota explicativa.

### 4.3 Columnas clave utilizadas del archivo principal

| Columna | Tipo | Rol en el modelo |
|---|---|---|
| **`Plant name`** | Texto | Identifica el aceite $i$. Cada especie única es una fila de la matriz $A$. |
| **`Chemical Name`** | Texto | Identifica el compuesto $j$. Cada compuesto único es una columna de $A$. |
| **`Percentage`** | Numérico (%) | Es el parámetro $a_{ij}$ — el dato central del modelo. Se promedia entre todas las muestras ESO de esa especie. |
| **`Plant Part`** | Texto | Filtro de calidad: solo usar registros de la parte estándar (flor para lavanda, hoja para eucalipto, etc.) para no mezclar perfiles distintos. |

### 4.4 Dato faltante y cómo resolverlo

- **Naranja dulce (*Citrus sinensis*):** no disponible en sCentInDB. Se sustituirá por *Citrus aurantium* (naranja amarga, disponible en la BD) o se complementará con datos de doTERRA Source-to-You (reportes de lote GC/MS públicos).

---

## 5. Modelo matemático

### 5.1 Variables de decisión

$$x_i \in [0, 1], \quad i = 1, \ldots, n$$

donde $x_i$ es la fracción volumétrica del aceite $i$ en la mezcla final.

### 5.2 Parámetros

| Símbolo | Descripción | Fuente |
|---|---|---|
| $a_{ij}$ | % del compuesto $j$ en el aceite $i$ (promedio entre muestras ESO) | sCentInDB |
| $c_j$ | Peso terapéutico del compuesto $j$ para el escenario elegido | Literatura (Tisserand & Young, 2014) |
| $p_i$ | Precio por mL del aceite $i$ | Catálogo doTERRA 2025 |
| $d_j^{max}$ | Concentración máxima permitida del compuesto $j$ | IFRA / Tisserand & Young |
| $u_i$ | Proporción máxima del aceite $i$ en la mezcla (e.g. 40%) | Buenas prácticas aromaterapia |
| $B$ | Presupuesto máximo por mL de mezcla | Definido por escenario |

### 5.3 Función objetivo

$$\max_{x} \quad \sum_{j=1}^{m} c_j \left( \sum_{i=1}^{n} a_{ij}\, x_i \right) = \max_{x} \quad c^T A^T x$$

### 5.4 Restricciones

$$\sum_{i=1}^{n} x_i = 1 \tag{mezcla completa}$$

$$x_i \geq 0 \quad \forall i \tag{no negatividad}$$

$$\sum_{i=1}^{n} a_{ij}\, x_i \leq d_j^{max} \quad \forall j \tag{seguridad por compuesto}$$

$$x_i \leq u_i \quad \forall i \tag{límite por aceite}$$

$$\sum_{i=1}^{n} p_i\, x_i \leq B \tag{presupuesto (opcional)}$$

### 5.5 Forma estándar para `scipy.optimize.linprog`

Dado que `linprog` minimiza, se convierte el problema:

$$\min_{x} \quad -c^T A^T x \quad \text{s.a.} \quad A_{ub}\,x \leq b_{ub},\; A_{eq}\,x = b_{eq},\; 0 \leq x \leq u$$

---

## 5-Alt. Solución alternativa: Método de Puntos Interiores (Barrera Logarítmica)

> Esta sección es complementaria al modelo de la Sección 5. Se usa el **mismo problema PL**, pero se resuelve con un algoritmo distinto para comparar convergencia, número de iteraciones y valor óptimo obtenido.

### 5-Alt.1 Idea central del método

El método de puntos interiores transforma el problema con restricciones de desigualdad en una **sucesión de problemas sin restricciones** mediante una función de barrera logarítmica que penaliza acercarse a los bordes del conjunto factible.

Para cada restricción de desigualdad $g_k(x) \leq 0$, se añade al objetivo el término:

$$-\mu \sum_{k} \ln(-g_k(x))$$

donde $\mu > 0$ es el **parámetro de barrera**. El algoritmo reduce $\mu \to 0$ iterativamente, de modo que la solución de cada subproblema converge a la solución óptima del problema original.

### 5-Alt.2 Formulación con barrera para nuestro problema

Partiendo de la forma estándar con variables de holgura $s \geq 0$:

$$A_{ub}\,x + s = b_{ub}, \quad s \geq 0$$

El subproblema de barrera en cada iteración $t$ es:

$$\min_{x,\, s} \quad -c^T A^T x - \mu \sum_{k} \ln(s_k)$$

$$\text{s.a.} \quad A_{ub}\,x + s = b_{ub}, \quad A_{eq}\,x = b_{eq}$$

Las condiciones de optimalidad (KKT) del subproblema dan lugar al **sistema de Newton** que se resuelve en cada paso:

$$\begin{pmatrix} 0 & A_{ub}^T & A_{eq}^T \\ A_{ub} & -S^{-1}Z & 0 \\ A_{eq} & 0 & 0 \end{pmatrix} \begin{pmatrix} \Delta x \\ \Delta s \\ \Delta \lambda \end{pmatrix} = -\begin{pmatrix} r_d \\ r_p \\ r_e \end{pmatrix}$$

donde $S = \text{diag}(s)$, $Z = \text{diag}(\lambda_s)$, y $r_d, r_p, r_e$ son los residuos de dualidad, primalidad y complementariedad.

### 5-Alt.3 Criterio de paro

El algoritmo termina cuando se satisfacen simultáneamente:

$$\frac{\|r_p\|}{1 + \|b\|} \leq \varepsilon, \quad \frac{\|r_d\|}{1 + \|c\|} \leq \varepsilon, \quad \frac{x^T s}{n} \leq \varepsilon$$

con tolerancia típica $\varepsilon = 10^{-8}$.

### 5-Alt.4 Implementación en Python

`scipy.optimize.linprog` soporta puntos interiores de forma nativa cambiando un solo parámetro:

```python
from scipy.optimize import linprog

# Símplex revisado (método principal)
res_simplex = linprog(c_obj, A_ub=A_ub, b_ub=b_ub,
                      A_eq=A_eq, b_eq=b_eq,
                      bounds=bounds,
                      method='highs-ds')   # dual simplex

# Puntos interiores (método alternativo)
res_ipm = linprog(c_obj, A_ub=A_ub, b_ub=b_ub,
                  A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds,
                  method='highs-ipm')      # interior point
```

### 5-Alt.5 Tabla de comparación Símplex vs. Puntos Interiores

| Criterio | Símplex (revisado) | Puntos Interiores |
|---|---|---|
| **Trayectoria** | Recorre vértices del poliedro factible | Atraviesa el interior del poliedro |
| **Iteraciones** | Pocas en problemas pequeños | Más iteraciones, cada una más costosa |
| **Convergencia** | Exacta en un vértice | Asintótica (se acerca al óptimo) |
| **Variables duales** | Directamente disponibles al terminar | Requieren recuperación del sistema KKT |
| **Sensibilidad** | Análisis de rango exacto | Aproximado vía perturbación |
| **Complejidad teórica** | Exponencial (peor caso) | Polinomial |
| **Aplicación en este proyecto** | Solución principal + análisis de sensibilidad | Verificación del óptimo y comparación |

### 5-Alt.6 Qué se compara en el notebook

Para cada uno de los 3 escenarios se reporta:

| Métrica | Símplex | Puntos Interiores |
|---|---|---|
| Valor óptimo $z^*$ | — | — |
| Solución $x^*$ (fracciones de mezcla) | — | — |
| Número de iteraciones | — | — |
| Tiempo de cómputo (ms) | — | — |
| Restricciones activas | — | — |

> Se espera que ambos métodos lleguen al **mismo $z^*$** (mismo óptimo global), pero por caminos distintos. La diferencia principal aparecerá en el número de iteraciones y en la "suavidad" de la trayectoria hacia el óptimo.

---

## 6. Los tres escenarios de optimización

### Escenario A — Relajación

- **Objetivo:** maximizar linalool + acetato de linalilo (actividad ansiolítica y sedante)
- **Aceites candidatos:** lavanda, romero, limón, árbol de té, incienso
- **Compuestos objetivo:** linalool, linalyl acetate
- **Restricción especial:** linalool ≤ 10% en mezcla final (límite IFRA piel)

### Escenario B — Antimicrobiano / Respiratorio

- **Objetivo:** maximizar eucaliptol + 4-terpineol + α-pineno
- **Aceites candidatos:** eucalipto, árbol de té, romero, menta, incienso
- **Compuestos objetivo:** eucalyptol (1,8-cineole), 4-terpineol, α-pineno
- **Restricción especial:** eucaliptol ≤ 35% en mezcla (seguridad neurológica en niños)

### Escenario C — Energizante / Mental

- **Objetivo:** maximizar limoneno + mentol + alcanfor
- **Aceites candidatos:** limón, menta, romero, lavanda, naranja/citrus
- **Compuestos objetivo:** limonene, menthol, camphor
- **Restricción especial:** mentol ≤ 5% en mezcla (límite en cosmética facial)

---

## 7. Estructura del notebook (`.ipynb`)

### Sección 1 — Introducción *(~0.5 páginas)*
- Descripción del problema y motivación
- Vínculo con aromaterapia y doTERRA
- Objetivo del proyecto

### Sección 2 — Marco teórico *(~2 páginas)*
- Programación Lineal: forma estándar, geometría del poliedro factible
- Método Símplex: pasos, tableaux, pivoteo
- Condiciones KKT: lagrangianos, variables duales, interpretación económica (precios sombra)
- Comparativa Símplex vs. Puntos Interiores (opcional)

### Sección 3 — Datos *(~2 páginas + código)*
- Descripción de sCentInDB
- Carga y exploración de CSVs (`pandas`)
- Limpieza: normalización de nombres de compuestos, manejo de múltiples muestras ESO
- Construcción de la matriz $A$ (aceites × compuestos): promedio por especie
- Visualización: **heatmap** de composición química (seaborn)

### Sección 4 — Modelado matemático *(~1.5 páginas)*
- Definición formal de $x$, $A$, $c$, $d^{max}$, $p$
- Justificación de los pesos $c_j$ por escenario (referencias bibliográficas)
- Precios doTERRA por mL
- Límites de seguridad IFRA / Tisserand & Young
- Conversión a forma estándar `linprog`

### Sección 5 — Solución *(código + resultados)*
- Implementación con `scipy.optimize.linprog`
- Solución de los 3 escenarios
- Verificación de factibilidad y optimalidad
- Tabla de resultados: $x^*$, valor objetivo $z^*$, restricciones activas

### Sección 6 — Interpretación *(~2 páginas + gráficas)*
- Gráfica de barras apiladas: mezcla óptima por escenario
- Variables duales: ¿cuál restricción "cuesta más" relajar?
- Análisis de sensibilidad: variación de $B$ (presupuesto) y $d_j^{max}$ (límites de seguridad)
- Comparación entre los 3 escenarios

### Sección 7 — Conclusiones *(~0.5 páginas)*
- Resultados principales
- Limitaciones del modelo (un solo lote promedio, precios estimados)
- Extensiones posibles: programación entera mixta para elegir subconjunto de aceites

---

## 8. Librerías y herramientas

```python
pandas          # carga y limpieza de CSVs
numpy           # álgebra lineal, construcción de A, b, c
scipy.optimize  # linprog (Símplex + Interior Point)
matplotlib      # gráficas de resultados
seaborn         # heatmap de composición
pulp            # (opcional) interfaz alternativa más legible para PL
```

**Solver:** `scipy.optimize.linprog` con `method='highs'` (predeterminado, usa HiGHS internamente — soporta Símplex revisado y puntos interiores).

---

## 9. Plan de trabajo y cronograma

| Semana | Fechas | Tarea | Entregable interno |
|:---:|---|---|---|
| 1 | 6–11 may | Exploración de datos, limpieza de CSVs, construcción de matriz $A$ | `data_exploration.ipynb` |
| 2 | 12–18 may | Marco teórico en Markdown/LaTeX, definición formal del modelo | Secciones 1–4 del notebook |
| 3 | 19–25 may | Implementación de los 3 escenarios, verificación de resultados | Sección 5 completa |
| 4 | 26–31 may | Visualizaciones, análisis de sensibilidad, conclusiones | Notebook final `.ipynb` |
| 5 | 1–7 jun | Preparación de presentación PDF, ensayo de 10–12 min | `.pdf` para entrega |

---

## 10. Riesgos y mitigaciones

| Riesgo | Probabilidad | Mitigación |
|---|:---:|---|
| *Citrus sinensis* no está en sCentInDB | **Ocurrido** | Usar *Citrus aurantium* (sí disponible) o datos de doTERRA Source-to-You |
| *Boswellia sacra/carterii* no en BD | **Ocurrido** | Usar *B. serrata* con nota; composición similar en α-pineno y terpenos |
| Problema infactible por restricciones muy estrictas | Media | Relajar $d_j^{max}$ en 10–20% o eliminar aceites con composición atípica |
| Solución degenerada (varios óptimos) | Baja | Añadir criterio de desempate: minimizar costo como objetivo secundario |

---

## 11. Referencias

1. Samal, A. et al. (2026). sCentInDB: a database of essential oil chemical profiles of Indian medicinal plants. *Molecular Diversity*. https://doi.org/10.1007/s11030-025-11215-5
2. Tisserand, R. & Young, R. (2014). *Essential Oil Safety* (2nd ed.). Churchill Livingstone.
3. Soto, J. E. (2026). Notas del curso MAT-34420 — Sección 4.1: Programación Lineal y Método Símplex. ITAM. https://itam-ds.github.io/analisis-numerico-computo-cientifico/
4. IFRA (International Fragrance Association). IFRA Standards — Concentration limits by compound and application category. https://ifrafragrance.org
5. doTERRA International. Source to You — GC/MS batch testing reports. https://sourcetoyou.doterra.com
6. Virtanen, P. et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
7. Leal, W. S. (2013). Odorant reception in insects. *Annual Review of Entomology* (referencia de apoyo para clasificación de compuestos).
