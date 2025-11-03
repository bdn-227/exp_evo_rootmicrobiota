

# rarefaction of variants

# imports
import polars as pl
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt



# dummy dataset
variants_long = pl.read_csv(f"raw_data/variants.raw.long.tsv", separator="\t")


# ~~~~~~~~~~~~ LOAD METADATA ~~~~~~~~~~~~ #
mapping = pl.read_csv("mapping/mapping.txt", separator="\t").sort("generation")



# we keep the original columns ref_fw, ref_rv, alt_fw and alt_rv
# and calcluate their frequencies from the whole pool of variants
# then we "subsample" the reads by defining a lower set of total 
# reads we might get and calculate the expected number of
# reads we are getting. we pass these through the filtering and
# and also keep the raw data and count them both

def calculate_results(variants_raw):
    return variants_raw.group_by("sampleID", "sampling_depth", "rare_frac").len().rename({"len": "#variants"})

# step to use
alignment_step = "trimmed_reads"

# get the information about the sequencing depth
# get the sequencing depth
depth_df = (pl.read_csv(f"processedData/coverage.tsv", separator="\t")
                    .pivot(index = "sampleID", columns="processing_step", values="#Reads")
                    .select("sampleID", "trimmed_reads", "aligned_to_syncom_analysis_ready"))


# now get the support for each of the variants
variant_support = (variants_long
                    .select(["sampleID", "ID_line", "ref_fw", "ref_rv", "alt_fw", "alt_rv"])
                    .join(depth_df, how="left", on=["sampleID"])
                    .with_columns( (pl.col( "ref_fw", "ref_rv", "alt_fw", "alt_rv")/pl.col(alignment_step)) )
                    )


# define parameters for subsampling
depth_min = variant_support[alignment_step].min()
depth_max = variant_support[alignment_step].max()
depth_range = np.linspace(0, depth_max+1, 100)

# collect results
results_filtered = []
results_raw = []



for idx, sampling_depth in enumerate(depth_range):

    # print progress
    print(f"{idx/depth_range.shape[0]}")

    # cacluate rarefaction fractions and detection frequencies
    variants_raw = (variant_support
                        .with_columns(pl.lit(sampling_depth).alias("sampling_depth"))
                        .with_columns( ((pl.col("sampling_depth") / pl.col(alignment_step)).alias("rare_frac") ))
                        .with_columns((pl.col(["ref_fw", "ref_rv", "alt_fw", "alt_rv"])*pl.col("sampling_depth")))
                        .with_columns( (pl.col("ref_fw") + pl.col("ref_rv") + pl.col("alt_fw") + pl.col("alt_rv")).alias("DP"))
                        .filter(pl.col("rare_frac") <= 1)
                        .filter(pl.col("DP") >= 10)
                    )


    # break the loop of we are over 1
    if variants_raw.shape[0] == 0:
        continue

    
    # filter the variants
    variants_filtered = filterVariants(variants_raw)

    

    # get the counts for variants
    variants_raw = calculate_results(variants_raw)
    variants_filtered = calculate_results(variants_filtered)
    
    
    # collect the results
    results_filtered.append(variants_filtered)
    results_raw.append(variants_raw)


# concat the results
rare_results_filtered = pl.concat(results_filtered, how="vertical").join(mapping, how="left", on="sampleID")
rare_results_raw = pl.concat(results_raw, how="vertical").join(mapping, how="left", on="sampleID")


# save the rarefaction values
rare_results_filtered.write_csv(f"processedData/rarefied_total_freq.filtered.tsv", separator="\t")
rare_results_raw.write_csv(f"processedData/rarefied_total_freq.raw.tsv", separator="\t")


# convert to pandas for plotting
rare_results_filtered   = rare_results_filtered.to_pandas()
rare_results_raw        = rare_results_raw.to_pandas()


# 
rare_results_raw.loc[(rare_results_raw["generation"] == 12) & (rare_results_raw["plant"] == "col0") & (rare_results_raw["syncom"] == "lj")].groupby("sampleID").max()

rare_results_raw["sampling_depth"].max()


# Plot the lines on two facets
ax = sns.relplot(
                data=rare_results_filtered, 
                x="sampling_depth", 
                y="#variants", 
                hue="sampleID",
                legend=False,
                col="plant",
                row = "syncom",
                kind="line", 
                height=5, 
                aspect=1.5, 
                facet_kws=dict(sharex=True),
                )

plt.xlim([0, depth_max*1.1])
plt.tight_layout()
plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.filtered.pdf")
plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.filtered.png")
plt.clf()



# Plot the lines on two facets
ax = sns.relplot(
                data=rare_results_raw, 
                x="sampling_depth", 
                y="#variants", 
                hue="sampleID",
                legend=False,
                col="plant",
                row = "syncom",
                kind="line", 
                height=5, 
                aspect=1.5, 
                facet_kws=dict(sharex=True),
                )

plt.xlim([0, depth_max*1.1])
plt.tight_layout()
plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.raw.pdf")
plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.raw.png")
plt.clf()



# generation specific plot
for generation in rare_results_raw["generation"].unique():
    rare_results_raw_susbet = rare_results_raw.loc[rare_results_raw["generation"] == generation]
    # Plot the lines on two facets
    ax = sns.relplot(
                    data=rare_results_raw_susbet, 
                    x="sampling_depth", 
                    y="#variants", 
                    hue="sampleID",
                    legend=True,
                    col="plant",
                    row = "syncom",
                    kind="line", 
                    height=5, 
                    aspect=1.5, 
                    facet_kws=dict(sharex=True),
                    )

    plt.xlim([0, depth_max*1.1])
    plt.tight_layout()
    plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.g{generation}.raw.pdf")
    plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.g{generation}.raw.png")
    plt.clf()



for generation in rare_results_filtered["generation"].unique():
    rare_results_filtered_subset = rare_results_filtered.loc[rare_results_filtered["generation"] == generation]
    # Plot the lines on two facets
    ax = sns.relplot(
                    data=rare_results_filtered_subset, 
                    x="sampling_depth", 
                    y="#variants", 
                    hue="sampleID",
                    legend=True,
                    col="plant",
                    row = "syncom",
                    kind="line", 
                    height=5, 
                    aspect=1.5, 
                    facet_kws=dict(sharex=True),
                    )

    plt.xlim([0, depth_max*1.1])
    plt.tight_layout()
    plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.g{generation}.filtered.pdf")
    plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.g{generation}.filtered.png")
    plt.clf()


cmap =          [
                 [(0,0,0)],
                 sns.color_palette("Greys", 5)[1:],
                 sns.color_palette("Reds", 5)[1:],
                 sns.color_palette("Blues", 5)[1:],
                 sns.color_palette("Greens", 5)[1:],
                ]
from itertools import chain
cmap = list(chain(*cmap))



for line in rare_results_raw["line"].unique():

    # subset and get correct colors
    rare_results_raw_subset = rare_results_raw.loc[rare_results_raw["line"] == line]
    rare_results_raw_subset = rare_results_raw_subset.sort_values("generation")

    # get colors
    cmap_ = {f"g{gen}_{line}":color for gen, color in enumerate(cmap)}

    # Plot the lines on two facets
    ax = sns.relplot(
                    data=rare_results_raw_subset, 
                    x="sampling_depth", 
                    y="#variants", 
                    hue="sampleID",
                    legend=True,
                    col="plant",
                    palette=cmap_,
                    row = "syncom",
                    kind="line", 
                    height=5, 
                    aspect=1.5, 
                    facet_kws=dict(sharex=True),
                    )

    plt.xlim([0, depth_max*1.1])
    plt.tight_layout()
    plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.{line}.raw.pdf")
    plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.{line}.raw.png")
    plt.clf()




# plot the variants that need to be resequenced
reseq = pl.read_ods("../reseq_samples.ods")
reseq = rare_results_raw.loc[rare_results_raw["sampleID"].isin(reseq["reseq"])]


# Plot the lines on two facets
ax = sns.relplot(
                data=reseq, 
                x="sampling_depth", 
                y="#variants", 
                hue="sampleID",
                legend=True,
                col="plant",
                row = "syncom",
                kind="line", 
                height=5, 
                aspect=1.5, 
                facet_kws=dict(sharex=True),
                )

plt.xlim([0, depth_max*1.1])
plt.tight_layout()
plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.reseq.raw.pdf")
plt.savefig(f"figures/rarefaction/rare_total_freq/rarefied_freq.reseq.raw.png")
plt.clf()
