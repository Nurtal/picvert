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


def plot_umap_visit():
    """ """

    # load data
    df_p = pd.read_csv("data/umap/precisesads.csv")
    df_g = pd.read_csv("data/umap/geo_visit.csv")
    df_all = pd.concat([df_p, df_g], ignore_index=True)

    # plot
    plt.figure(figsize=(8, 6))
    for label, group in df_all.groupby('LABEL'):
        plt.scatter(group['X'], group['Y'], label=str(label), s=20, alpha=0.7)
    plt.title("UMAP - Visit")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend(title="LABEL", loc='best')
    plt.grid(True, alpha=0.2)
    plt.tight_layout()

    # save
    plt.savefig('images/umap_visit.png')
    plt.close()


def plot_umap_treatment():
    """ """

    # load data
    df_p = pd.read_csv("data/umap/precisesads.csv")
    df_g = pd.read_csv("data/umap/geo_treatment.csv")
    df_all = pd.concat([df_p, df_g], ignore_index=True)

    # plot
    plt.figure(figsize=(8, 6))
    for label, group in df_all.groupby('LABEL'):
        plt.scatter(group['X'], group['Y'], label=str(label), s=20, alpha=0.7)
    plt.title("UMAP - Treatment")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend(title="LABEL", loc='best')
    plt.grid(True, alpha=0.2)
    plt.tight_layout()

    # save
    plt.savefig('images/umap_treatment.png')
    plt.close()


def plot_trajectory():
    """ """

    df_p = pd.read_csv("data/umap/precisesads.csv")
    df_t = pd.read_csv("data/umap/geo_treatment.csv")
    df_v = pd.read_csv("data/umap/geo_visit.csv")
    df_t['VISIT'] = df_v['LABEL']
    df_p['VISIT'] = "NA"
    df_all = pd.concat([df_p, df_t], ignore_index=True)

    time_point = set(list(df_v['LABEL']))
    for t in time_point:
        df = df_all[df_all['VISIT'].isin([t, "NA"])]
        df = df.drop(columns=['VISIT'])
    
        plt.figure(figsize=(8, 6))
        for label, group in df.groupby('LABEL'):
            plt.scatter(group['X'], group['Y'], label=str(label), s=20, alpha=0.7)
        plt.title(f"UMAP - Treatment at {t}")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend(title="LABEL", loc='best')
        plt.grid(True, alpha=0.2)
        plt.tight_layout()

        plt.savefig(f"images/umap_treatment_at_{t.replace(' ', '_')}.png")
        plt.close()



def plot_trajectory_clean():
    """ """

    df_p = pd.read_csv("data/umap/precisesads.csv")
    df_t = pd.read_csv("data/umap/geo_treatment.csv")
    df_v = pd.read_csv("data/umap/geo_visit.csv")
    df_t['VISIT'] = df_v['LABEL']
    df_p['VISIT'] = "NA"
    df_all = pd.concat([df_p, df_t], ignore_index=True)

    time_point = set(list(df_v['LABEL']))
    for t in time_point:
        df = df_all[df_all['VISIT'].isin([t, "NA"])]
        df = df.drop(columns=['VISIT'])
        df = df[df['LABEL'].isin(['treatment group: Cyclophosphamide', 'treatment group: Transplant'])]
    
        plt.figure(figsize=(8, 6))
        for label, group in df.groupby('LABEL'):
            plt.scatter(group['X'], group['Y'], label=str(label), s=20, alpha=0.7)
        plt.title(f"UMAP - Treatment at {t}")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend(title="LABEL", loc='best')
        plt.grid(True, alpha=0.2)
        plt.tight_layout()

        plt.savefig(f"images/umap_treatment_at_{t.replace(' ', '_')}_clean.png")
        plt.close()


if __name__ =="__main__":

    plot_umap()
    plot_umap_visit()
    plot_umap_treatment()
    plot_trajectory()
    plot_trajectory_clean()
