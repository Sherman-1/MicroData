import polars as pl 
import matplotlib.pyplot as plt 
import seaborn as sns 
import numpy as np
from tqdm import tqdm

print(f"Polars version: {pl.__version__}")
print(f"NumPy version: {np.__version__}")
print(f"Seaborn version: {sns.__version__}")
print(f"Matplotlib version: {plt.matplotlib.__version__}")


globular_vs_afdb = pl.scan_csv("Globular_vs_afbd_micro.tsv", separator = "\t").with_columns(
    query = pl.col("query").str.split(".").list.get(0),
    target = pl.col("target").str.replace(r"AF-(.*)-F1.*", r"${1}"),
    evalue = pl.col("evalue").cast(pl.Float64),
    qcov = pl.col("qcov").cast(pl.Float64),
    tcov = pl.col("tcov").cast(pl.Float64)
).filter(
    (pl.col("qcov") >= 0.7) & (pl.col("tcov") >= 0.7)
)

molten_vs_afdb = pl.scan_csv("Molten_vs_afbd_micro.tsv", separator = "\t").with_columns(
    query = pl.col("query").str.split(".").list.get(0),
    target = pl.col("target").str.replace(r"AF-(.*)-F1.*", r"${1}"),
    evalue = pl.col("evalue").cast(pl.Float64),
    qcov = pl.col("qcov").cast(pl.Float64),
    tcov = pl.col("tcov").cast(pl.Float64)
).filter(
    (pl.col("qcov") >= 0.7) & (pl.col("tcov") >= 0.7)
)

evalue_thresholds = [1, 1e-1, 1e-2, 1e-3, 1e-6, 1e-9]
overlapping_targets = []

for threshold in tqdm(evalue_thresholds, desc="Thresholds", total = len(evalue_thresholds)):
    globular_filtered = globular_vs_afdb.filter(pl.col("evalue") <= threshold)
    molten_filtered = molten_vs_afdb.filter(pl.col("evalue") <= threshold)
    overlapping = globular_filtered.join(molten_filtered, on="target", how="inner").select(pl.len()).collect().item()
    overlapping_targets.append(overlapping)

print(overlapping_targets)

print("Saving figure")
plt.figure(figsize=(10, 6))
plt.plot(evalue_thresholds, overlapping_targets, marker='o')
plt.xscale('log')
plt.xlabel('E-value threshold')
plt.ylabel('Number of overlapping targets')
plt.title('Evolution of overlapping targets with increasing E-value thresholds')
plt.savefig("overlapping_targets_evolution.png")
plt.close()

micro_plddt_stats = pl.scan_csv("afdb_microproteins_plddt_stats.tsv", separator = "\t").with_columns(
    id = pl.col("id").str.split("-").list.get(1)
)


globular = globular_vs_afdb.join(micro_plddt_stats, left_on = "target", right_on = "id").with_columns(
    category = pl.lit("globular")
)
molten = molten_vs_afdb.join(micro_plddt_stats, left_on = "target", right_on = "id").with_columns(
    category = pl.lit("molten")
)
macro_plddt_stats = pl.scan_csv("afdb_macro_sample_plddt_stats.tsv", separator = "\t").with_columns(
    id = pl.col("id").str.split("-").list.get(1),
    category = pl.lit("macro")
)

plddt_data = pl.concat([
    globular.select(["category","median_plDDT","fraction_res_above_70","protein_length"]),
    molten.select(["category","median_plDDT","fraction_res_above_70","protein_length"]),
    macro_plddt_stats.select(["category","median_plDDT","fraction_res_above_70","protein_length"])
]).filter(
    pl.int_range(pl.len()).shuffle().over("category") < 100000
).collect().to_pandas()

plt.figure(figsize = (10,6))
sns.boxplot(data = plddt_data, x = "category", y = "median_plDDT")
plt.xlabel("Category")
plt.ylabel("Median plDDT")
plt.title("Median plDDT distribution")
plt.savefig("median_plDDT_distribution.png")
plt.close()

plt.figure(figsize = (10,6))
sns.boxplot(data = plddt_data, x = "category", y = "fraction_res_above_70")
plt.xlabel("Category")
plt.ylabel("Median plDDT")
plt.title("Median plDDT distribution")
plt.savefig("fraction_res_above_70_distribution.png")
plt.close()

plddt_data = pl.concat(
    [
        micro_plddt_stats.with_columns(category = pl.lit("micro")).filter(
            pl.int_range(pl.len()).shuffle().over("category") < 100000
        ),
        macro_plddt_stats.with_columns(category = pl.lit("macro"))
    ]
).collect().to_pandas()

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sns.kdeplot(data=plddt_data, hue="category", x="fraction_res_above_70", ax=axes[0])
axes[0].set_xlabel("Fraction of residues above 70")
axes[0].set_ylabel("Density")
axes[0].set_title("Fraction of residues above 70 distribution : Micro vs Macro")

sns.kdeplot(data=plddt_data, hue="category", x="median_plDDT", ax=axes[1])
axes[1].set_xlabel("Median plDDT")
axes[0].set_ylabel("Density")
axes[1].set_title("Median plDDT distribution : Micro vs Macro")

plt.tight_layout()
plt.savefig("Micro_vs_Macro.png")
plt.close()

