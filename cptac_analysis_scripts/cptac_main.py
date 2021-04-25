import cptac
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

from cptac_tsne_exploration import plot_tsne


def plot_ages():
    clinical_data = pd.read_csv("./data/cptac_clinical.csv", index_col=0)
    clinical_data["age_in_years"] = clinical_data["Age.in.Month"] // 12
    fig, ax = plt.subplots(1, 1)

    sns.distplot(clinical_data["age_in_years"], hist=True, kde=True,
                 bins=20)

    median_age = clinical_data["age_in_years"][len(clinical_data["age_in_years"]) // 2]
    ax.set(xlabel="Age (years)", ylabel="")
    plt.suptitle(f"Median age: {median_age} years")
    plt.show()


def main():
    # plot_tsne(colorby="PAM50", title="PAM50 Type")
    plot_tsne(colorby="Patient Age Category", title="Age")
    plot_tsne(colorby="Stage (Simplified)", title="Stage")


if __name__ == "__main__":
    main()
