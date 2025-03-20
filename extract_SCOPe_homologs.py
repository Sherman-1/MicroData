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
)

molten_vs_afdb = pl.scan_csv("Molten_vs_afbd_micro.tsv", separator = "\t").with_columns(
    query = pl.col("query").str.split(".").list.get(0),
    target = pl.col("target").str.replace(r"AF-(.*)-F1.*", r"${1}")
).filter(
    (pl.col("qcov") >= 0.7) & (pl.col("tcov") >= 0.7) & (pl.col("evalue") <= 1e-3)
)

micro_plddt_stats = pl.scan_csv("afdb_microproteins_plddt_stats.tsv", separator = "\t").with_columns(
    short_id = pl.col("id").str.split("-").list.get(1)
)

moltens_hit = molten_vs_afdb.select("target").collect().unique("target").to_series().to_list()
globulars_hit = globular_vs_afdb.select("target").collect().unique("target").to_series().to_list()

micro_plddt_stats = micro_plddt_stats.with_columns(
    hit_molten = pl.col("short_id").is_in(moltens_hit),
    hit_globular = pl.col("short_id").is_in(globulars_hit)
).with_columns(
    hit = (pl.col("hit_molten") | pl.col("hit_globular")),
    median_plDDT = pl.col("median_plDDT").cast(pl.Float32), 
    fraction_res_above_70 = pl.col("fraction_res_above_70").cast(pl.Float32)
).collect()

print(micro_plddt_stats.head)

plt.figure(figsize=(20,10))
sns.boxenplot(data=micro_plddt_stats.to_pandas(), y="median_plDDT", hue="hit")
plt.title("Distribution of median pLDDT scores for AFDB Microproteins")
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), 
            title = "AFDB MicroProtein hit \nagainst SCOPe")

plt.savefig("afdb_microproteins_plddt_stats.png")

plt.figure(figsize=(20,10))
sns.boxenplot(data=micro_plddt_stats.to_pandas(), y="fraction_res_above_70", hue="hit")
plt.title("Distribution of fraction of residues with pLDDT > 70 for AFDB Microproteins")
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), 
            title = "AFDB MicroProtein hit \nagainst SCOPe")

plt.savefig("afdb_microproteins_plddt_stats.png")

# Extract globular hits
globular_hits = micro_plddt_stats.with_columns(
    id = pl.concat_str([pl.lit("/datas/SIMON/AFDB_20_to_100/pdbs_microproteins/"), pl.col("id"), pl.lit(".pdb")])
).filter(
    (pl.col("hit_globular") == True) &
    (pl.col("hit_molten") == False) &
    (pl.col("median_plDDT") > 70) &
    (pl.col("fraction_res_above_70") > 0.5)
).select("id").unique().write_csv("globular_structural_homologs.txt", include_header = False)

# Extract molten hits
micro_plddt_stats.with_columns(
    id = pl.concat_str([pl.lit("/datas/SIMON/AFDB_20_to_100/pdbs_microproteins/"), pl.col("id"), pl.lit(".pdb")])
).filter(
    (pl.col("hit_molten") == True) &
    (pl.col("hit_globular") == False) &
    (pl.col("median_plDDT") > 70) &
    (pl.col("fraction_res_above_70") > 0.5)
).select("id").unique().write_csv("molten_structural_homologs.txt", include_header = False)


macro_plddt_stats = pl.scan_csv("Macroproteins/afdb_macro_sample_plddt_stats.tsv", separator = "\t").with_columns(category = pl.lit("MacroProtein"))
micro_plddt_stats = pl.scan_csv("afdb_microproteins_plddt_stats.tsv", separator = "\t").with_columns(category = pl.lit("MicroProtein"))

data = pl.concat([macro_plddt_stats, micro_plddt_stats]).select(["category","median_plDDT", "fraction_res_above_70"]).collect()

plt.figure(figsize=(20,10))
sns.set(style="darkgrid")
sns.violinplot(data=data.to_pandas(), hue="category", y="median_plDDT")
plt.title("Distribution of median pLDDT scores for AFDB Proteins")
plt.savefig("median_plDDT_distribution.png")

plt.figure(figsize=(20,10))
sns.set(style="darkgrid")
sns.violinplot(data=data.to_pandas(), hue="category", y="fraction_res_above_70")
plt.title("Distribution of fraction of residues with pLDDT > 70 for AFDB Proteins")
plt.savefig("fraction_res_above_70_distribution.png")