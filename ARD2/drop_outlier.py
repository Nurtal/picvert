import pandas as pd
from sklearn.ensemble import IsolationForest
from scipy.spatial import distance
import numpy as np


def isolation_forest():
    """ """

    # load data
    df = pd.read_csv("data/data_vst_lda.csv")
    df = df.drop(columns=['Unnamed: 0'])
    X = df.select_dtypes(include="number")

    # Modèle d'isolation forest
    iso = IsolationForest(contamination=0.1, random_state=42)  # contamination = % attendu d'outliers
    df["outlier"] = iso.fit_predict(X)

    # Dans sklearn, outlier = -1, normal = 1
    print(df["outlier"].value_counts())

    # Supprimer les outliers
    df_clean = df[df["outlier"] == 1].drop(columns="outlier")

    # save
    df_clean.to_csv("data/data_vst_outliers_droped_lda.csv", index=False)


def dist_outlier():
    """ """

    # load data
    df = pd.read_csv("data/data_vst_lda.csv")
    df = df.drop(columns=['Unnamed: 0'])
    X = df.select_dtypes(include="number")

    # Moyenne et matrice de covariance
    mean = np.mean(X, axis=0)
    cov = np.cov(X, rowvar=False)
    inv_cov = np.linalg.inv(cov)

    # Distance de Mahalanobis pour chaque patient
    dists = [distance.mahalanobis(x, mean, inv_cov) for x in X]

    df["mahalanobis"] = dists

    # Seuil = par ex. au-dessus du 97.5e percentile
    threshold = np.percentile(dists, 97.5)
    df["outlier"] = df["mahalanobis"] > threshold

    df_clean = df[~df["outlier"]].drop(columns=["mahalanobis", "outlier"])
    print(df_clean)
    
    # save
    df_clean.to_csv("data/data_vst_outliers_droped_lda.csv", index=False)


if __name__ == "__main__":

    isolation_forest()
    # dist_outlier()
