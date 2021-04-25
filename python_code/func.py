import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


def make_a_plot(file, x, y):
    rho, spear_pvalue = stats.spearmanr(x, y)
    rho_check, spear_pvalue_check = stats.spearmanr(x, y)

    assert rho == rho_check

    # Step 3: Plot the data
    # You can present the rho value in various ways. Table, scatterplot, etc.
    # Here you will make a scatterplot

    plt.figure(figsize=(10, 10))

    # Replace x and y with appropriate variables
    plt.scatter(x, y)

    title = "rho: {} for ESR1".format(rho)  # This is string formatting. The variable in the () will print in the {}
    plt.title(title)

    # Fill in informative x and y labels
    plt.xlabel("RNA expression levels")
    plt.ylabel("Protein expression levels")

    # plt.show()  # Comment out when running in script
    plt.savefig(f"figs/{file}.png", bbox_inches="tight")  # Use this when saving figure in script


if __name__ == "__main__":
    protein_data = pd.read_csv("data/cptac_protein.csv", index_col=0).loc[:, "ESR1"]
    rna_data = pd.read_csv("data/cptac_rna.csv", index_col=0).loc[:, "ESR1"]
    make_a_plot("scatter_func_thing", protein_data, rna_data)
