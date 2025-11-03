
import polars as pl


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
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



def load_colors():
    return pl.read_csv(f"mapping/syncom_colors.tsv", separator="\t")

def color_dict_family():
    colors = load_colors()
    return {family: color for family, color in zip(colors["family"].to_numpy(), colors["color"].to_numpy())}

def color_dict_strain():
    colors = load_colors()
    return {strain: color for strain, color in zip(colors["strain"].to_numpy(), colors["color"].to_numpy())}


def get_plant(plant):
    if plant == "col0":
        return "$A. thaliana$ (Col-0)"
    if plant == "gifu":
        return "$L. japonicus$ (Gifu)"
    else:
        return "*"


def get_marker(syncom, plant):
    if syncom == "at" and plant == "col0":
        return "s"
    if syncom == "at" and plant == "gifu":
        return "D"
    if syncom == "lj" and plant == "col0":
        return "D"
    if syncom == "lj" and plant == "gifu":
        return "s"
