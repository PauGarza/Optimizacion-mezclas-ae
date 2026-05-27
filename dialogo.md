# Diálogo por diapositiva

---

## Slide 1 — Portada


---

## Slide 2 — ¿Qué problema resolvemos?

- Cada aceite tiene una composición química analizada por cromatografía GC-MS — cientos de compuestos por aceite
- La idea es que ciertos compuestos están asociados a efectos terapéuticos; queremos maximizar la presencia de esos compuestos en la mezcla
- La mezcla va diluida al 10% en un roll-on de 25 gotas, así que las proporciones tienen que sumar 1
- Las restricciones son tres: seguridad (límites IFRA por compuesto), proporción máxima por aceite, y presupuesto
- Traducido a LP: variables = fracción de cada aceite, objetivo = beneficio terapéutico ponderado, restricciones = las de la tabla

---

## Slide 3 — Datos — Matriz A

- La matriz A tiene 30 aceites como filas y ~800 compuestos como columnas
- Cada celda es el porcentaje promedio de ese compuesto en ese aceite, promediado sobre todas las muestras GC-MS disponibles
- El 85% de las celdas son cero — la matriz es muy dispersa, porque cada aceite tiene sus compuestos característicos y no comparte muchos con los demás
- Eso en la práctica ayuda al solver: hay menos restricciones activas simultáneamente

---

## Slide 4 — EDA: Precios y Cobertura GC-MS

- El rango de precios por gota es amplio: desde el limoncillo a $0.99/gota hasta el sándalo a $20.65
- Eso hace que el presupuesto sea una restricción activa solo en fórmulas que incluyen aceites caros
- El otro problema que encontramos en los datos: 9 de los 30 aceites tienen suma de fila mayor a 110%
- Eso pasa porque los porcentajes GC-MS se reportan sobre distintas fracciones de la muestra y se acumulan
- Antes de meterlos al LP los normalizamos a 100%, de lo contrario el modelo sobreestima el beneficio

---

## Slide 5 — Construcción del vector c

- El vector c es el que le dice al modelo qué compuestos son importantes para un uso terapéutico específico
- Lo construimos con ScentInDB: una base de datos que asocia plantas a usos terapéuticos
- Para cada escenario tomamos las plantas de referencia, promediamos sus perfiles GC-MS y normalizamos
- Por ejemplo para "ansiolítico" el peso cae fuerte en linalyl acetate y linalool; para "digestivo" en cuminaldehyde y piperitone oxide
- En total construimos 59 vectores c distintos — uno por escenario terapéutico

---

## Slide 6 — Algoritmo Símplex

- El símplex recorre vértices del poliedro factible — y en este problema eso se traduce directamente a mezclas donde la mayoría de los aceites están en exactamente 0% o en su cota máxima
- Por eso las soluciones son dispersas: típicamente 5 o 6 aceites activos de los 30 posibles
- Usamos la variante dual símplex porque tenemos muchas restricciones de desigualdad (los límites IFRA)
- Un ajuste que tuvimos que hacer: con U_MAX_CAP = 0.40 el solver llenaba solo 3 aceites al tope y la mezcla quedaba muy concentrada; al bajarlo a 0.20 forzamos más diversidad
- Converge en 15 a 35 iteraciones para nuestros 30 aceites

---

## Slide 7 — Puntos Interiores

- A diferencia del símplex, el IPM no salta entre vértices: empieza en el interior del poliedro y se mueve por un camino central hacia el óptimo
- Para mantenerse adentro usa una función barrera logarítmica que penaliza acercarse a las fronteras
- El parámetro μ controla qué tanto pesa esa barrera; al irlo reduciendo a cero, el camino converge al vértice óptimo
- Durante las iteraciones intermedias todos los aceites tienen fracciones positivas — es una mezcla "densa" que después colapsa
- Converge en 8 a 15 iteraciones de Newton, pero cada iteración es más costosa que en símplex para este tamaño de problema

---

## Slide 8 — Resultados — Mezclas Óptimas

- Veamos tres escenarios concretos
- Para ansiolítico la mezcla usa 5 aceites: pimienta negra, geranio, albahaca, lavanda y limón — todos al 20% con limón completando el resto
- Para soporte de sueño aparecen pachulí, albahaca, lavanda, nardo, pimienta negra y ylang ylang — es la fórmula más cara, $3.71/gota
- Para digestivo los aceites son hierbabuena, hinojo, limón, cilantro, cardamomo y limoncillo — la más barata a $1.94/gota
- Lo interesante es que la composición refleja la bioquímica: aceites conocidos por esos usos aparecen de forma natural en la solución

---

## Slide 9 — Transición (¿Qué pasa si simulamos todos los escenarios?)

- Eso fue con tres escenarios específicos; pero tenemos 59
- ¿Hay patrones? ¿Qué aceites son los más versátiles? ¿Cuánto varía el costo?

---

## Slide 10 — Los 59 escenarios — Panorama general

- El scatter muestra beneficio vs costo para los 59 escenarios — el color indica cuántos aceites activos tiene la mezcla
- No hay una correlación fuerte: pagar más no garantiza mayor beneficio terapéutico
- Hay un outlier con beneficio muy alto pero costo bajo — esos son escenarios donde los compuestos clave están en aceites baratos
- En la gráfica de frecuencia: pimienta negra aparece en 33 de 59 escenarios, limón en 28, romero en 23
- Son aceites "versátiles" — sus compuestos aparecen en muchos vectores c distintos
- Algo que detectamos: eugenol y timol son las restricciones de seguridad que se activan más seguido — la seguridad limita más que el presupuesto en la mayoría de los casos
- También encontramos 12 pares de escenarios con el mismo vector c — nombres distintos en ScentInDB pero mismas plantas de referencia

---

## Slide 11 — Dualidad y Comparación de Métodos

- Las variables duales nos dan los precios sombra: cuánto mejoraría el óptimo si relajamos cada restricción en una unidad
- Para ansiolítico con B = $3.5/gota, la restricción de eugenol tiene λ = 0.284 — relajar el límite de eugenol en 1 punto porcentual mejora el beneficio en 0.284
- La restricción de timol tiene λ = 0 — no está activa, relajarla no cambia nada
- Eso confirma que eugenol es el cuello de botella de seguridad, no el presupuesto
- En cuanto a los métodos: símplex e IPM llegan exactamente al mismo óptimo — la diferencia es menor a 10⁻¹⁰
- La diferencia de complejidad teórica (exponencial vs polinomial) no importa aquí con solo 30 variables; importaría a escala de miles

---

## Slide 12 — Conclusiones y Limitaciones

- El modelo funciona: la mezcla óptima refleja conocimiento real de aromaterapia sin que se lo dijéramos explícitamente
- La seguridad IFRA es el límite más frecuentemente activo, no el costo
- Símplex e IPM son equivalentes para este tamaño — elegir uno u otro es indiferente en la práctica
- Hay un punto de saturación de presupuesto alrededor de $4.5/gota; después de ahí pagar más no mejora la fórmula
- En limitaciones: el modelo es un proxy de co-ocurrencia, no prueba eficacia clínica
- No captura sinergias — dos compuestos juntos pueden tener un efecto distinto al de cada uno por separado
- 6 aceites del catálogo no tienen datos GC-MS en ScentInDB y quedaron fuera

---

## Slide 13 — Gracias / Preguntas


---
