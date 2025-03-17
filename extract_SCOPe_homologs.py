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
    target = pl.col("target").str.replace(r"AF-(.*)-F1.*", r"${1}")
).filter(
    (pl.col("qcov") >= 0.7) & (pl.col("tcov") >= 0.7) & (pl.col("evalue") <= 1e-3)
).collect()

print("Globular vs AFDB microproteins collected")

molten_vs_afdb = pl.scan_csv("Molten_vs_afbd_micro.tsv", separator = "\t").with_columns(
    query = pl.col("query").str.split(".").list.get(0),
    target = pl.col("target").str.replace(r"AF-(.*)-F1.*", r"${1}")
).filter(
    (pl.col("qcov") >= 0.7) & (pl.col("tcov") >= 0.7) & (pl.col("evalue") <= 1e-3)
).collect()

print("Molten vs AFDB microproteins collected")

print("Computing overlapping targets")

overlapping = globular_vs_afdb.join(molten_vs_afdb, on="target", how="inner").select("target").to_series().unique().to_list()


print(f"Number of globular targets before removing overlaps : {globular_vs_afdb.select('target').n_unique('target')}")
print(f"Number of molten targets before removing overlaps : {molten_vs_afdb.select('target').n_unique('target')}")
print(f"Number of overlapping targets : {len(overlapping)}")


molten_vs_afdb = molten_vs_afdb.filter(~pl.col("target").is_in(overlapping))
globular_vs_afdb = globular_vs_afdb.filter(~pl.col("target").is_in(overlapping))

print(f"Number of globular targets after removing overlaps : {globular_vs_afdb.select('target').n_unique('target')}")
print(f"Number of molten targets after removing overlaps : {molten_vs_afdb.select('target').n_unique('target')}")

micro_plddt_stats = pl.scan_csv("afdb_microproteins_plddt_stats.tsv", separator = "\t").with_columns(
    id = pl.col("id").str.split("-").list.get(1)
).filter(
    (pl.col("fraction_res_above_70") >= 0.5) & 
    (pl.col("median_plDDT") >= 70) 
).collect()

globular = globular_vs_afdb.join(micro_plddt_stats, left_on = "target", right_on = "id").with_columns(
    category = pl.lit("globular")
)
molten = molten_vs_afdb.join(micro_plddt_stats, left_on = "target", right_on = "id").with_columns(
    category = pl.lit("molten")
)

globular_structural_homologs = globular.unique("target").select("target").to_series().to_list()
molten_structural_homologs = molten.unique("target").select("target").to_series().to_list()

print(f"Number of globular structural homologs : {len(globular_structural_homologs)}")
print(f"Number of molten structural homologs : {len(molten_structural_homologs)}")

print("Saving structural homologs")
with open("globular_structural_homologs.txt", "w") as f:
    f.write("\n".join(globular_structural_homologs))

with open("molten_structural_homologs.txt", "w") as f:
    f.write("\n".join(molten_structural_homologs))