import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap
import matplotlib.pyplot as plt
import numpy as np



def run_semi_supervised_umap():
    """ """

    # load data
    df_unlabeled = pd.read_csv("data/data_normalized2_lda.csv").drop(columns=['ID'])
    df_labeled = pd.read_csv("data/precisesads/input_normalized2_lda.csv").drop(columns=['ID'])
    df_labeled['LABEL'] = df_labeled['LABEL'].replace({'clust1':1, 'clust2':2, 'clust3':3})
    target = list(df_labeled['LABEL'])
    classes = ['cluster1', 'cluster2', 'cluster3']

    # combine
    X_all = pd.concat([df_labeled.drop('LABEL', axis=1), df_unlabeled], ignore_index=True)
    X_all = X_all.dropna()
    print(X_all)
    
    # Extract labels
    y_all = pd.concat([
        df_labeled['LABEL'],
        pd.Series([-1] * len(df_unlabeled))
        
    ], ignore_index=True)

    # Standardisation
    X_scaled = StandardScaler().fit_transform(X_all)

    # UMAP semi-supervised
    umap_model = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        target_metric='categorical',
        target_weight=0.8
    )

    # Compute
    embedding = umap_model.fit_transform(X_scaled, y=y_all)

    fig, ax = plt.subplots(1, figsize=(14, 10))
    plt.scatter(*embedding.T, s=0.1, c=y_all, cmap='Spectral', alpha=1.0)
    plt.setp(ax, xticks=[], yticks=[])
    cbar = plt.colorbar(boundaries=np.arange(11)-0.5)
    cbar.set_ticks(np.arange(3))
    cbar.set_ticklabels(classes)
    plt.title('Tadam')

    plt.show()



def run_projection_umap():
    """Train a umap, project new data in it and save embeddings in a data/umap subfolder"""

    # load data
    df_unlabeled = pd.read_csv("data/data_normalized2_lda.csv").drop(columns=['ID'])
    df_unlabeled['LABEL'] = 4
    df_labeled = pd.read_csv("data/precisesads/input_normalized2_lda.csv").drop(columns=['ID'])
    df_labeled['LABEL'] = df_labeled['LABEL'].replace({'clust1':1, 'clust2':2, 'clust3':3})
    target = list(df_labeled['LABEL'])
    classes = ['cluster1', 'cluster2', 'cluster3']


    # process data
    X = df_labeled.drop('LABEL', axis=1)
    X = X.dropna()
    X_new = df_unlabeled.drop('LABEL', axis=1)
    X_new = X_new.dropna()
    
    # Extract labels
    y = df_labeled['LABEL']

    # Standardisation
    X_scaled = StandardScaler().fit_transform(X)
    X_new_scaled = StandardScaler().fit_transform(X_new)

    # create mapper
    mapper = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        target_metric='categorical',
        target_weight=0.8
    ).fit(X_scaled, y)

    # compute new emnedding
    new_embedding = mapper.transform(X_new_scaled) 

    # save PRECISESADS
    data = []
    cmpt = 0
    for c in mapper.embedding_:
        data.append({"X":c[0], "Y":c[1], "LABEL":y[cmpt]})
        cmpt+=1
    df = pd.DataFrame(data)
    df.to_csv("data/umap/precisesads.csv")

    # save GEO
    data = []
    for c in new_embedding:
        data.append({"X":c[0], "Y":c[1], "LABEL":"GEO"})
    df = pd.DataFrame(data)
    df.to_csv("data/umap/geo.csv")



if __name__ == "__main__":


    # run_semi_supervised_umap()
    run_projection_umap()

