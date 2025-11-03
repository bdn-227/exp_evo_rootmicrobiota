
# script to test for an intersection between the pacbio sequencing and shotgun sequencing-based variant calling


# modules
import polars as pl
import numpy as np

from shotgun_analysis_functions import load_mapping

# mapping
mapping = load_mapping()


# mapping table of evolved strains
evolved_strains_df = (pl.read_csv("mapping/evolved_strains.tsv", separator="\t")
                        .with_columns((pl.col("strain")+"_"+pl.col("ID")).alias("evolved_strain"))
                        .rename({"ID": "reiso_id", "syncom": "syncom_id", "evolved_strain": "target_strain"})
                        .drop("strain", "family", "condition", "plant")
                        .with_columns(pl.col("line").str.strip_chars_start("l"))
                     )

# shotgun data
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                            separator="\t",
                            schema_overrides={"line": str}
                            )
                            .filter(pl.col("line").is_in(["1","2","3","4","5",
                                                          "6.2","7.2","8.2","9","10",
                                                          "11","12","13","14","15",
                                                          "16","17","18","19","20"]))
                )




# MAPPING USING MINIMAP
pacbio_deepvariant = (pl.read_csv("raw_data/pacbio_variants_deepvariant.tsv", separator="\t", schema_overrides={"line": str})
                              .with_columns(pl.col("evolved_strain").str.split("_").list[1].alias("reiso_id"))
                              .with_columns(pl.col("reiso_id").replace("12-45", "2-45"))
                              .join(evolved_strains_df, how="left", on=["reiso_id"])
                              .with_columns(pl.col("line").str.strip_chars("l"))
                              .with_columns((pl.col("chrom") + "_" + pl.col("pos").cast(str) + "_" + pl.col("ref") + "_" + pl.col("alt")).alias("ID"))
                     )


# concatenated metadata table
variant_long_meta_data = variants_long.select("ID", "geneID", "alias", "desc", "ko").unique()
pacbio_deepvariant_meta_data = pacbio_deepvariant.select("ID", "geneID", "alias", "desc", "ko").unique()
metadata = pl.concat([variant_long_meta_data, pacbio_deepvariant_meta_data]).unique()


var_ls = []
var_intersect_ls = []
for line in np.intersect1d(pacbio_deepvariant["line"].drop_nulls(), variants_long["line"]):
    pacbio_deepvariant_line = pacbio_deepvariant.filter(pl.col("line")==line)
    seq_strains = pacbio_deepvariant_line["strain"].unique().to_list()
    variants_long_line = (variants_long
                                    .filter(pl.col("line")==line)
                                    .filter(pl.col("generation")==16)
                                    .filter(pl.col("strain").is_in(seq_strains))
                         )
    

    # intersection to match variants
    variants_long_line = variants_long_line.select("syncom", "plant", "chrom", "pos", "ref", "alt", "ID", "line", "generation", "AF").unique()
    pacbio_deepvariant_line = pacbio_deepvariant_line.select("chrom", "pos", "ref", "alt", "ID", "line", "af_deepvariant").unique()
    variants_intersect = variants_long_line.join(pacbio_deepvariant_line, how="inner", on=["chrom", "pos", "ref", "alt", "ID", "line"])
    variants_long_lost = variants_long_line.filter( ~pl.col("ID").is_in(variants_intersect["ID"].to_list()) )
    syncom=variants_long.filter(pl.col("line")==line)["syncom"].unique()[0]
    plant=variants_long.filter(pl.col("line")==line)["plant"].unique()[0]
    variants_deepvar_lost = pacbio_deepvariant_line.filter( ~pl.col("ID").is_in(variants_intersect["ID"].to_list()) ).with_columns(pl.lit(syncom).alias("syncom")).with_columns(pl.lit(plant).alias("plant"))
    variants_all = pl.concat([variants_intersect, variants_long_lost, variants_deepvar_lost], how="diagonal")
    var_ls.append(variants_all)
    var_intersect_ls.append(variants_intersect)


variants_all = (pl.concat(var_ls)
                    .rename({"AF": "af_shotgun"})
                    .join(metadata, how="left", on="ID")
                    .sort("plant", "syncom", "chrom", "pos", "ref", "alt")
               )
variants_all.drop("chrom", "pos", "ref", "alt").write_csv("results/shotgun_pacbio_all.tsv", separator="\t")



variants_intersect = (pl.concat(var_intersect_ls)
                            .rename({"AF": "af_shotgun"})
                            .join(metadata, how="left", on="ID")
                            .sort("plant", "syncom", "chrom", "pos", "ref", "alt")
                    )
variants_intersect.drop("chrom", "pos", "ref", "alt").write_csv("results/shotgun_pacbio_intersect1d.tsv", separator="\t")


# fixerd ones
fixed_variants = variants_intersect.filter(pl.col("af_shotgun")>=0.95)
fixed_variants.write_csv("results/shotgun_pacbio_fixed.tsv", separator="\t")
