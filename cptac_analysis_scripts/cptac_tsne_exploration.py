import cptac
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


def mod_clinical(file="./data/cptac_clinical.csv"):
    clinical_data = pd.read_csv("./data/cptac_clinical.csv")

    clinical_data["age_in_years"] = clinical_data["Age.in.Month"] // 12

    age_cat = []
    for age in clinical_data["age_in_years"]:
        if age < 54:
            age_cat.append("Age < 54 Years")
        else:
            age_cat.append("Age > 54 Years")

    stage_simplified = []
    for stage in clinical_data["Stage"]:
        stage_simplified.append(str(stage).strip(" ABC"))

    clinical_data["Patient Age Category"] = age_cat
    clinical_data["Stage (Simplified)"] = stage_simplified

    clinical_data.to_csv("./data/cptac_clinical_mod.csv")


# creates tsne plots
def plot_tsne(colorby="Patient Age", title="Patient Age", save=True):
    # read in data
    protein_data = pd.read_csv("./data/cptac_protein.csv", index_col=0)
    rna_data = pd.read_csv("./data/cptac_rna.csv", index_col=0)
    clinical_data = pd.read_csv("./data/cptac_clinical_mod.csv")

    # check to make sure labels match
    assert all(protein_data.axes[0] == clinical_data["Patient_ID"])
    assert all(rna_data.axes[0] == clinical_data["Patient_ID"])

    # fill in NA values with 0
    protein_data = protein_data.fillna(0)
    rna_data = rna_data.fillna(0)

    # performs the TSNE clustering, seed is set to 0
    tsne_protein_dims = TSNE(n_components=2, random_state=0).fit_transform(protein_data)
    tsne_rna_dims = TSNE(n_components=2, random_state=0).fit_transform(rna_data)

    # create two side-by-side scatter plots
    # first is protein expression; second is mRNA
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    sns.scatterplot(
        x=tsne_protein_dims[:, 0],
        y=tsne_protein_dims[:, 1],
        hue=clinical_data[colorby],
        data=protein_data,
        legend="full",
        alpha=0.7,
        ax=ax[0])
    ax[0].set(xlabel="Protein TSNE 1", ylabel="Protein TSNE 2")

    sns.scatterplot(
        x=tsne_rna_dims[:, 0],
        y=tsne_rna_dims[:, 1],
        hue=clinical_data[colorby],
        data=rna_data,
        legend="full",
        alpha=0.7,
        ax=ax[1])
    ax[1].set(xlabel="RNA TSNE 1", ylabel="RNA TSNE 2")

    fig.suptitle(f"Clustering CPTAC Patients by {title}")

    # save to a file
    if save:
        plt.savefig(f"figs/tsne_clusterby_{title}.png", bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    plot_tsne()
