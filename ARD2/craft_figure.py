import pandas as pd
import datashader as ds
import colorcet as cc
import matplotlib.pyplot as plt




def plot_umap():
    """ """

    # load data
    df_p = pd.read_csv("data/umap/precisesads.csv")
    df_g = pd.read_csv("data/umap/geo.csv")
    df_all = pd.concat([df_p, df_g], ignore_index=True)

    # plot
    plt.figure(figsize=(8, 6))
    for label, group in df_all.groupby('LABEL'):
        plt.scatter(group['X'], group['Y'], label=str(label), s=20, alpha=0.7)
    plt.title("UMAP - trained on PRECISESADS clusters")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend(title="LABEL", loc='best')
    plt.grid(True, alpha=0.2)
    plt.tight_layout()

    # save
    plt.savefig('images/umap.png')
    plt.close()




if __name__ =="__main__":

    plot_umap()
