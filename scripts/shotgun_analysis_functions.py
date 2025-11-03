
# ~~~~~~~~~~~~ IMPORTS ~~~~~~~~~~~~ #
import polars as pl
import numpy as np
import warnings
warnings.filterwarnings("ignore")




def load_kegg():
    return pl.read_csv(f"mapping/kegg.txt", separator="\t")



def loadGff():
    gff = pl.read_csv("mapping/mixed_sphere.gff", separator="\t", comment_prefix="#", has_header=False, null_values= ".")
    gff.columns = ["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"]

    # extend attribute column
    attr_cols = ["geneID", "alias", "ko", "desc"]
    gff = gff.with_columns([pl.col("attribute").str.split_exact(";", 4).struct.rename_fields(attr_cols).alias("fields")]).unnest("fields")
    gff = gff.drop("attribute")

    for col in attr_cols:
        gff = gff.with_columns(pl.col(col).str.split_exact("=", n=1).struct[-1].alias(col))


    # extract the strain information
    gff = gff.with_columns(pl.col("chrom").str.split_exact("_", n=1).struct[0].alias("strain"))
    gff = gff.select(["chrom", "feature", "start", "end", "strand", "phase", "strain", "geneID", "ko", "alias", "desc"])
    return gff


def load_mapping():
    return pl.read_csv("mapping/shotgun_mapping.txt", separator="\t").sort("generation")

def load_colors():
    return pl.read_csv(f"mapping/syncom_colors.tsv", separator="\t")

def color_dict_family():
    colors = load_colors()
    return {family: color for family, color in zip(colors["family"].to_numpy(), colors["color"].to_numpy())}

def color_dict_strain():
    colors = load_colors()
    return {strain: color for strain, color in zip(colors["strain"].to_numpy(), colors["color"].to_numpy())}






def number_formatter(number, precision=3):
    """
    function to reformat numbers into nice scientific notation
    """
    
    if "e" in str(number):
        return np.format_float_scientific(number, precision=precision).replace("e", "×$10^{") + "}$"
    else:
        return np.round(number, precision)
    
def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "NS"

def get_plant_label(plant):
    if plant == "col0":
        return "$A. thaliana$ (Col-0)"
    if plant == "gifu":
        return "$L. japonicus$ (Gifu)"
    else:
        return "Input (ancestral)"
    
def get_syncom_label(syncom):
    if syncom == "at":
        return "$At$-SPHERE strains"
    if syncom == "lj":
        return "$Lj$-SPHERE strains"
    
def get_label_enriched(enriched):
    if enriched == "col0":
        return "$A. thaliana$ (Col-0)"
    if enriched == "gifu":
        return "$L. japonicus$ (Gifu)"
    else:
        return "n.s."
    
def get_label_condition(condition):
    if condition == "native":
        return "Native"
    if condition == "non-native":
        return "Non-native"
    else:
        return "n.s."




def pretty_exponent(number, digits=2):

    # first, convert to scientific notation
    exp_number = "{:e}".format(number)

    # apply digits
    splitted_num = exp_number.split("e")

    if len(splitted_num) > 1:
        exp_number = "".join([str(round(float(splitted_num[0]), digits)), "e", splitted_num[1]])

    # now convert to true exponential
    if "e+00" in exp_number:
        return exp_number.split("e")[0]
    else:
        return exp_number.replace("e", "×10$^{") + "}$"



def hex_to_rgb(hex_color):
    # 1. Clean and validate the input string
    hex_color = hex_color.lstrip('#')
    hex_len = len(hex_color)

    if hex_len not in (3, 6, 8):
        raise ValueError("Invalid hex string: must be 3, 6, or 8 characters long.")

    # 2. Parse the hex string based on its length
    if hex_len == 3:
        # Shorthand hex (#RGB -> #RRGGBB)
        r = int(hex_color[0] * 2, 16)
        g = int(hex_color[1] * 2, 16)
        b = int(hex_color[2] * 2, 16)
    elif hex_len == 6:
        # Standard hex (#RRGGBB)
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)
    else: # hex_len == 8
        # Hex with alpha (#RRGGBBAA)
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)

    return (r, g, b)



def percent_change(old, new):
    if old == 0:
        if new == 0:
            return 0.0
        return float('inf') if new > 0 else -float('inf')
    return (new - old) / old * 100.0
