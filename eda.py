"""
eda.py  —  Análisis Exploratorio de Datos
Proyecto : Optimización de mezclas de aceites esenciales (MAT-34420 ITAM)
Autoras  : Felipe Castro, Paulina Garza

Genera figuras en figures/ y estadísticas de resumen para el README.

Uso:
    python eda.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

CLEAN = Path("data/clean")
FIGS  = Path("figures")
FIGS.mkdir(exist_ok=True)

sns.set_theme(style="whitegrid", palette="muted", font_scale=1.05)
plt.rcParams["figure.dpi"] = 150

# ─── Cargar artefactos limpios ────────────────────────────────────────────────
A         = pd.read_csv(CLEAN / "matrix_A.csv",    index_col=0)
prices    = pd.read_csv(CLEAN / "prices_clean.csv")
scenarios = pd.read_csv(CLEAN / "scenarios_c.csv", index_col=0)
agg       = pd.read_csv(CLEAN / "eo_aggregated.csv")
therapy   = pd.read_csv(CLEAN / "therapy_clean.csv")

# Agregar precio inscrito por gota (16 gotas/mL, estándar doTERRA)
# precio_por_mL ya usa el precio de mayoreo (Costo Inscritos)
if "precio_por_gota" not in prices.columns:
    prices.insert(prices.columns.get_loc("precio_por_mL") + 1,
                  "precio_por_gota",
                  (prices["precio_por_mL"] / 16).round(4))
    prices.to_csv(CLEAN / "prices_clean.csv", index=False)
    print("prices_clean.csv: columna precio_por_gota agregada")

print(f"Matriz A:    {A.shape[0]} aceites × {A.shape[1]} compuestos")
print(f"Escenarios:  {len(scenarios)}")
print(f"Precios:     {len(prices)} aceites")
print()

# ─── Fig 1 — Precio inscrito por gota ────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 7))
p = prices.sort_values("precio_por_gota")
colors = sns.color_palette("Blues_d", len(p))
bars = ax.barh(p["ACEITES"], p["precio_por_gota"], color=colors)
ax.bar_label(bars, labels=[f"${v:.2f}" for v in p["precio_por_gota"]],
             padding=4, fontsize=7.5)
ax.set_xlabel("Precio inscrito por gota (MXN)")
ax.set_title("Costo por gota — aceites esenciales doTERRA\n(precio de mayoreo · 16 gotas/mL)")
ax.set_xlim(0, p["precio_por_gota"].max() * 1.22)
plt.tight_layout()
plt.savefig(FIGS / "fig1_precio_por_gota.png", bbox_inches="tight")
plt.close()
print("fig1_precio_por_gota.png")

# ─── Fig 2 — Heatmap matriz A (top 25 compuestos por abundancia media) ───────
top25 = A.mean(axis=0).nlargest(25).index
A_top = A[top25]

fig, ax = plt.subplots(figsize=(17, 8))
sns.heatmap(A_top, annot=True, fmt=".1f", annot_kws={"size": 7},
            cmap="YlOrRd", linewidths=0.3,
            cbar_kws={"label": "% GC-MS (promedio entre muestras)"},
            ax=ax)
ax.set_title("Matriz A — Composición química (% GC-MS)\nTop 25 compuestos por abundancia promedio", fontsize=12)
ax.set_xlabel("")
ax.set_ylabel("Aceite doTERRA")
plt.xticks(rotation=45, ha="right", fontsize=7.5)
plt.tight_layout()
plt.savefig(FIGS / "fig2_heatmap_matriz_A.png", bbox_inches="tight")
plt.close()
print("fig2_heatmap_matriz_A.png")

# ─── Fig 3 — Cobertura GC-MS por aceite (sumas de fila) ──────────────────────
row_sums = A.sum(axis=1).sort_values(ascending=False)
colors   = ["#d73027" if v > 110 else "#4575b4" for v in row_sums]

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(range(len(row_sums)), row_sums.values, color=colors, width=0.7)
ax.axhline(100, color="black", linestyle="--", linewidth=1.5,
           label="100 % (cobertura completa esperada)")
ax.set_xticks(range(len(row_sums)))
ax.set_xticklabels(row_sums.index, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Suma de porcentajes GC-MS (%)")
ax.set_title("Cobertura GC-MS por aceite\n"
             "Rojo = suma > 110 % (normalizar antes del LP); azul = dentro del rango normal")
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(FIGS / "fig3_cobertura_gcms.png", bbox_inches="tight")
plt.close()
print("fig3_cobertura_gcms.png")

# ─── Fig 4 — Top 20 compuestos más ubicuos (presentes en más aceites) ─────────
presence = (A > 0).sum(axis=0).sort_values(ascending=False).head(20)

fig, ax = plt.subplots(figsize=(11, 5))
ax.bar(range(len(presence)), presence.values,
       color=sns.color_palette("Blues_d", len(presence)))
ax.set_xticks(range(len(presence)))
ax.set_xticklabels(presence.index, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Número de aceites (de 30 totales)")
ax.set_title("Top 20 compuestos más ubicuos\n"
             "(presentes en el mayor número de aceites doTERRA)")
ax.set_ylim(0, 30)
for i, v in enumerate(presence.values):
    ax.text(i, v + 0.3, str(v), ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(FIGS / "fig4_compuestos_ubicuos.png", bbox_inches="tight")
plt.close()
print("fig4_compuestos_ubicuos.png")

# ─── Fig 5 — Vector c para escenarios representativos ────────────────────────
key_scenarios = {
    "Anxiolytics":    "Anxiolytics [TC-therapeutic category]",
    "Sleep Support":  "Sleep Support",
    "Digestive":      "Digestive Support",
    "Purifying":      "Purifying and Cleansing",
}
# filtrar a los que existen en scenarios_c
key_scenarios = {k: v for k, v in key_scenarios.items() if v in scenarios.index}

fig, axes = plt.subplots(1, len(key_scenarios), figsize=(5 * len(key_scenarios), 6))
if len(key_scenarios) == 1:
    axes = [axes]

palette = sns.color_palette("viridis", 10)
for ax, (short, full) in zip(axes, key_scenarios.items()):
    c = scenarios.loc[full]
    top10 = c.nlargest(10)
    ax.barh(range(len(top10)), top10.values, color=palette)
    ax.set_yticks(range(len(top10)))
    ax.set_yticklabels([t[:35] for t in top10.index], fontsize=8)
    ax.set_xlabel("Peso $c_j$ (normalizado)")
    ax.set_title(short, fontsize=10)
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.3f}"))

plt.suptitle("Vector c — top 10 compuestos por escenario terapéutico",
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(FIGS / "fig5_vectores_c.png", bbox_inches="tight")
plt.close()
print("fig5_vectores_c.png")

# ─── Fig 6 — Distribución de n_samples en la tabla agregada ──────────────────
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(agg["n_samples"], bins=30, color="#4575b4", edgecolor="white", linewidth=0.5)
ax.set_xlabel("Número de mediciones promediadas por par (especie, compuesto)")
ax.set_ylabel("Frecuencia")
ax.set_title("Calidad del promedio GC-MS\n"
             "Distribución de n_samples en eo_aggregated.csv")
pct_1 = (agg["n_samples"] == 1).mean()
ax.axvline(1, color="red", linestyle="--", label=f"n=1 ({pct_1:.0%} de los pares)")
ax.legend()
plt.tight_layout()
plt.savefig(FIGS / "fig6_nsamples.png", bbox_inches="tight")
plt.close()
print("fig6_nsamples.png")

# ─── Estadísticas de resumen ──────────────────────────────────────────────────
sparsity    = (A == 0).sum().sum() / A.size
row_sums    = A.sum(axis=1)
over100     = (row_sums > 110).sum()

print("\n" + "=" * 55)
print("ESTADISTICAS PARA README")
print("=" * 55)
print(f"Matriz A:          {A.shape[0]} aceites × {A.shape[1]} compuestos")
print(f"Sparsity:          {sparsity:.1%}  (celdas con 0)")
print(f"Aceites > 110%:    {over100} de {len(A)}  (requieren normalización)")
print(f"Compuesto con más presencia: "
      f"{presence.index[0]} ({presence.iloc[0]} aceites)")
print(f"Escenarios totales: {len(scenarios)}")
print(f"  scentindb:        {sum(1 for i in scenarios.index if 'Support' not in i and 'Care' not in i and 'Relief' not in i and 'Health' not in i and 'Repellent' not in i and 'Uplifting' not in i and 'Meditation' not in i and 'Oral' not in i and 'Insect' not in i and 'Purif' not in i and 'Hormonal' not in i and 'Focus' not in i and 'Grounding' not in i and 'Sleep' not in i and 'Hair' not in i and 'Skin' not in i)}")
print(f"Rango precio/gota: "
      f"${prices['precio_por_gota'].min():.2f} – ${prices['precio_por_gota'].max():.2f} MXN")
print(f"Aceite más caro/gota:   "
      f"{prices.loc[prices['precio_por_gota'].idxmax(), 'ACEITES']}  "
      f"(${prices['precio_por_gota'].max():.2f})")
print(f"Aceite más barato/gota: "
      f"{prices.loc[prices['precio_por_gota'].idxmin(), 'ACEITES']}  "
      f"(${prices['precio_por_gota'].min():.2f})")
print(f"\nTop 5 compuestos más ubicuos:")
for comp, n in presence.head(5).items():
    print(f"  {comp:<40s}  en {n}/30 aceites")

print("\nFiguras guardadas en figures/")
