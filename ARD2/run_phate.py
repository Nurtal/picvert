import phate
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.preprocessing import StandardScaler


def run(precisesads_file, geo_file, nb_voisins, nb_steps):

    # init output
    if not os.path.isdir("data/phate"):
        os.mkdir("data/phate")
        os.mkdir("data/phate/all")
        os.mkdir("data/phate/precisesads")
        os.mkdir("data/phate/geo")

    # load data
    df_unlabeled = pd.read_csv("data/data_normalized2_lda.csv").drop(columns=['ID'])
    df_unlabeled['LABEL'] = 4
    df_labeled = pd.read_csv("data/precisesads/input_normalized2_lda.csv").drop(columns=['ID'])
    df_labeled['LABEL'] = df_labeled['LABEL'].replace({'clust1':1, 'clust2':2, 'clust3':3})
    target = list(df_labeled['LABEL'])
    classes = ['cluster1', 'cluster2', 'cluster3']

    # process data
    X = df_labeled.dropna()
    y = list(X['LABEL'])
    print(len(y))
    X = X.drop('LABEL', axis=1)
    df_unlabeled['LABEL'] = 'geo'
    X_new = df_unlabeled.dropna()
    y_new = list(X_new['LABEL'])
    X_new = X_new.drop('LABEL', axis=1)
    y_all = y + y_new
    
    # Standardisation
    X_scaled = StandardScaler().fit_transform(X)
    X_new_scaled = StandardScaler().fit_transform(X_new)

    # aggregate
    X_total = np.vstack([X_scaled, X_new_scaled])

    # Appliquer PHATE - all
    phate_op = phate.PHATE(n_components=2, knn=nb_voisins, t=nb_steps)
    X_phate = phate_op.fit_transform(X_total)
    
    # rearange data for plot
    df = pd.DataFrame({
        'PHATE1': X_phate[:,0],
        'PHATE2': X_phate[:,1],
        'Label': y_all
    })

    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=df, x='PHATE1', y='PHATE2',
        hue='Label',       # couleur selon label
        palette='tab10',   # choix de palette, tu peux changer
        s=50,              # taille des points
        alpha=0.8
    )
    plt.title("PHATE embedding - ALL")
    plt.xlabel("PHATE1")
    plt.ylabel("PHATE2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f"data/phate/all/knn_{nb_voisins}_steps_{nb_steps}.png")

    # Appliquer PHATE - Precisesads
    phate_op = phate.PHATE(n_components=2, knn=nb_voisins, t=nb_steps)
    X_phate = phate_op.fit_transform(X_scaled)

    # rearange data for plot
    df = pd.DataFrame({
        'PHATE1': X_phate[:,0],
        'PHATE2': X_phate[:,1],
        'Label': y
    })

    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=df, x='PHATE1', y='PHATE2',
        hue='Label',       # couleur selon label
        palette='tab10',   # choix de palette, tu peux changer
        s=50,              # taille des points
        alpha=0.8
    )
    plt.title("PHATE embedding - PRECISESADS")
    plt.xlabel("PHATE1")
    plt.ylabel("PHATE2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f"data/phate/precisesads/knn_{nb_voisins}_steps_{nb_steps}.png")

    # Appliquer PHATE - GEO
    phate_op = phate.PHATE(n_components=2, knn=nb_voisins, t=nb_steps)
    X_phate = phate_op.fit_transform(X_new_scaled)

    # rearange data for plot
    df = pd.DataFrame({
        'PHATE1': X_phate[:,0],
        'PHATE2': X_phate[:,1],
        'Label': y_new
    })

    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=df, x='PHATE1', y='PHATE2',
        hue='Label',       # couleur selon label
        palette='tab10',   # choix de palette, tu peux changer
        s=50,              # taille des points
        alpha=0.8
    )
    plt.title("PHATE embedding - GEO")
    plt.xlabel("PHATE1")
    plt.ylabel("PHATE2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f"data/phate/geo/knn_{nb_voisins}_steps_{nb_steps}.png")



def explore():
    """ """

    for nb_voisins in range(3,30):
        for nb_steps in range(10, 150, 5):
            run("data/precisesads/lda_input.csv", "data/data_vst_lda.csv", nb_voisins, nb_steps)

    
if __name__ == "__main__":

    
    # run("data/precisesads/input_normalized2_lda.csv", "data/data_normalized2_lda.csv")
    # run("data/precisesads/lda_input.csv", "data/data_vst_lda.csv", 25, 25)
    explore()
