
# import modules
import polars as pl
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import false_discovery_control



# ~~~~~~~~~~~~ configs ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *
import warnings
warnings.filterwarnings("ignore")


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()
kegg = load_kegg()
gff = loadGff()



# and import the inputs
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                            separator="\t",
                            schema_overrides={"line": str}
                            )
                            .filter(pl.col("af_type")=="original")
                            .filter(pl.col("type")!="upstream")
                )


# reformat variants for the glm
index_cols = ["plant", "syncom", "condition", "treatment", "line", "generation"]
variants_glm = (variants_long
                    .pivot(on="ID", values="AF", index=index_cols)
                    .fill_null(0)
                    .unpivot(index=index_cols, variable_name="ID", value_name="AF"))


# get the lowest non-zero val
min_val = variants_glm.filter(pl.col("AF")!= 0)["AF"].min()
max_val = variants_glm.filter(pl.col("AF")!= 0)["AF"].max()
padding_val = (min_val/max_val)*min_val

# parameters for testing
formula = "AF ~ plant+generation"
family = sm.families.Gaussian()


# test each variant, whether it is enriched for one of the two hosts
glm_result_ls = []
for syncom in ["at", "lj"]:
    
    # subset the variant table
    variants_sc_subset = variants_glm.filter(pl.col("syncom") == syncom)

    sc_ls = []
    for variant_id in variants_sc_subset["ID"].unique():

        # subsetting the data
        variants = variants_sc_subset.filter(pl.col("ID") == variant_id)

        # only if non-zeros
        if (variants["AF"] != 0).any() and (not (variants["AF"] == 1).all()):
    
            # perform fitting of the glm
            glm_fit = smf.glm(formula=formula, data=variants.to_pandas(), family=family).fit()
            params = glm_fit.params
            pvalue = glm_fit.pvalues

            # get the fold change
            padded_af = variants.with_columns( pl.when(pl.col("AF") == 0).then(padding_val).otherwise("AF").alias("AF") )
            plant_ratio = np.log2(padded_af.filter(pl.col("plant") == "gifu")["AF"].mean() / padded_af.filter(pl.col("plant") == "col0")["AF"].mean())

            glm_res = (pl.DataFrame({"syncom": syncom, "ID": variant_id, "ratio_plant": plant_ratio})
                                .with_columns(pl.lit(value).alias(f"coef_{name}") for name, value in params.items())
                                .with_columns(pl.lit(value).alias(f"pval_{name}") for name, value in pvalue.items()))
            sc_ls.append(glm_res)

    # concat the syncom specific data
    glm_res_sc = pl.concat(sc_ls)

    # perform FDR correction of pvalues
    pval_cols = [col for col in glm_res_sc.columns if "pval" in col]
    glm_res_sc = (glm_res_sc
                        .with_columns(pl.col(col).map_batches(function = lambda x: false_discovery_control(x)).alias(f"padj_{col}") for col in pval_cols)
                        .with_columns( pl.when( (pl.col('padj_pval_plant[T.gifu]') < 0.05) & (pl.col("ratio_plant").abs() > 1) ).then(True).otherwise(False).alias("significant") )
                        .with_columns( pl.when(pl.col("coef_plant[T.gifu]") > 0).then(pl.lit("gifu")).otherwise(pl.lit("col0")).alias("enriched") )
                        .with_columns( pl.when(pl.col("significant")).then("enriched").otherwise(pl.lit("NS")).alias("enriched") )
                        )
    
    glm_result_ls.append(glm_res_sc)


# concat the corrected data 
glm_result_df = pl.concat(glm_result_ls)
print(f"fitted glms: " + str(glm_result_df["ID"].unique().len()))

# aggregate and plot significant hits
glm_result_df = (glm_result_df.with_columns( pl.when((pl.col("enriched") == "col0") & (pl.col("syncom") == "at")).then(pl.lit("native"))
                                               .when((pl.col("enriched") == "gifu") & (pl.col("syncom") == "lj")).then(pl.lit("native"))
                                               .otherwise(pl.lit("non-native")).alias("condition") )
                            .with_columns( pl.when(pl.col("significant")).then("condition").otherwise(pl.lit("NS")).alias("condition") )
                )




# native and non-native variants
condi_specific = glm_result_df.filter(pl.col("significant")).group_by("condition").len()
increase = round(((condi_specific.filter(pl.col("condition")=="non-native")["len"][0] - condi_specific.filter(pl.col("condition")=="native")["len"][0])/condi_specific.filter(pl.col("condition")=="native")["len"][0])*100, 2)
print(f'host-specific: {condi_specific["len"].sum()} out of {glm_result_df["ID"].n_unique()}: {round(condi_specific["len"].sum()/glm_result_df["ID"].n_unique()*100, 2)}%')
print(condi_specific)
print(f"increase in significant hits in non-native {increase}%")


# add the variant meta data
glm_result_df_meta = (glm_result_df
                        .join(variants_long.select("ID", "alias", "ko", "type", "geneID").unique(), how="left", on="ID")
                        .with_columns((-pl.col('padj_pval_plant[T.gifu]').log(base=10)).alias("-log10(pval)"))
                        .sort("syncom")
                        .with_columns( (pl.when( (pl.col("ko").is_null()) | (pl.col("alias") == "hypothetical protein") ).then(False).otherwise(True)).alias("annotated") )
                        .with_columns( pl.when(pl.col("ko") == pl.col("alias")).then(False).otherwise("annotated").alias("annotated") )
                        .with_columns( (pl.when(pl.col("alias") == "hypothetical protein").then("geneID").otherwise("alias")).alias("alias") )
                        .with_columns( (pl.when(pl.col("alias").is_null()).then("ID").otherwise("alias")).alias("alias") )
                        .with_columns((pl.col("ID") + "~" + pl.col("alias")).alias("ID_alias"))
                        )


# save the table
glm_result_df_meta.write_csv("results/volcano.tsv", separator="\t")


# reformat the final table with the mutations
cols = ["ID", "strain", "chrom", "start", 
        "plant", "syncom", "padj", 
        "type",
        "ko", "alias", "desc", 
        "a", "b", "c", ]

mapped_genes = (glm_result_df_meta
                .join(kegg.group_by("ko").agg(pl.col("a").str.join(" / "), pl.col("b").str.join(" / "), pl.col("c").str.join(" / ")).unique(), how="left", on="ko")
                .join(gff.select("geneID", "start").unique(), how="left", on="geneID")
                .filter(pl.col("significant"))
                .join(kegg.select("ko", "c", "desc").group_by("ko", "desc").agg(pl.col("c").str.join(",")), how="left", on="ko")
                .with_columns(pl.col("ID").str.split("_").list[0].alias("strain"))
                .with_columns((pl.col("strain") + "_contig_" + pl.col("ID").str.split("_").list[2]).alias("chrom"))
                .with_columns(pl.col("ID").str.split("_").list[3].alias("location"))
                .rename({"enriched": "plant", 'padj_pval_plant[T.gifu]': "padj"})
                .select(cols)
                .sort("plant", "syncom", "ID", "type")
                )

# save the results to drive
mapped_genes.write_csv("results/S_significant_variants.tsv", separator="\t")











