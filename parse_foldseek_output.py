import polars as pl 


lf = pl.scan_csv("Globular_vs_afbd_micro.tsv", separator = "\t").select(
    ["query","target","evalue","qcov","tcov"]
)


filter = (pl.col("evalue") < 1e-3) & (pl.col("qcov") > 0.7) & (pl.col("tcov") > 0.7) 

lf = lf.with_columns(
    target = pl.lit("/datas/SIMON/AFDB_20_to_100/pdbs_microproteins/") + pl.col("target"), 
    query = pl.lit("/datas/SIMON/AFDB_20_to_100/globular/") + pl.col("query")
)

lf = lf.filter(filter)

paths = lf.select("query").unique().collect().to_series().to_list() + lf.select("target").unique().collect().to_series().to_list()

with open("globular_paths.txt", "w") as f:
    for path in paths:
        f.write(path + "\n")

lf = pl.scan_csv("Molten_vs_afbd_micro.tsv", separator = "\t").select(
    ["query","target","evalue","qcov","tcov"]
)

lf = lf.with_columns(
    target = pl.lit("/datas/SIMON/AFDB_20_to_100/pdbs_microproteins/") + pl.col("target"), 
    query = pl.lit("/datas/SIMON/AFDB_20_to_100/molten/") + pl.col("query")
)

lf = lf.filter(filter)

paths = lf.select("query").unique().collect().to_series().to_list() + lf.select("target").unique().collect().to_series().to_list()

with open("molten_paths.txt", "w") as f:
    for path in paths: 
        f.write(path + "\n")