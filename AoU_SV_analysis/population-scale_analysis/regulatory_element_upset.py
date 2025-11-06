import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

populations = ["EAS", "SAS", "AFR", "AMR", "EUR"]

# Input files
variants_df = pd.read_csv("SV_summary.txt", sep="\t")
sv_elements_df = pd.read_csv(
    "SV_regulatory_element_ID.txt",
    sep="\t",
    header=None,
    names=["Variant_ID", "Regulatory_Element"]
)

elements_data = {}

for pop in populations:
    selected_variants = variants_df[
        (variants_df["Sample_Count"] >= 1) &
        (variants_df["Populations"].str.contains(pop))
    ]
    selected_variant_ids = selected_variants["Variant_ID"].tolist()

    selected_sv_elements = sv_elements_df[
        sv_elements_df["Variant_ID"].isin(selected_variant_ids)
    ]

    final_output = selected_sv_elements[["Regulatory_Element"]].drop_duplicates()

    output_filename = f"{pop}_regulatory_elements.txt"
    final_output.to_csv(output_filename, sep="\t", index=False, header=False)

    elements_data[pop] = set(final_output["Regulatory_Element"])

all_elements = sorted(set().union(*elements_data.values()))

data = pd.DataFrame(
    {
        pop: [element in elements_data[pop] for element in all_elements]
        for pop in populations
    },
    index=all_elements
)

upset_input = from_indicators(populations, data[populations])

upset = UpSet(upset_input, subset_size='count', sort_by='cardinality')

fig = plt.figure(figsize=(11, 4.7), dpi=600)
upset.plot(fig=fig)

plt.ylabel("Number of regulatory elements")
plt.tight_layout()
plt.savefig("regulatory_elements_upset.png", dpi=600)
plt.savefig("regulatory_elements_upset.pdf", format="pdf", dpi=600)
plt.close(fig)
