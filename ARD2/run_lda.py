import joblib
import pandas as pd
import matplotlib.pyplot as plt




def run_lda_normalize():
    """ """

    # load model
    lda = joblib.load("model/lda_normalized2.joblib")

    # load new data
    df = pd.read_csv("data/data_normalized2_lda.csv")
    df['LABEL'] = "machin"
    df['ID'] = "truc"
    df = df.drop(columns=['ID', 'LABEL'])
    X = df
    X1_proj = lda.transform(X)

    # load precisesads data
    df = pd.read_csv("data/precisesads/input_normalized2_lda.csv")
    df['LABEL'] = "machin"
    df['ID'] = "truc"
    df = df.drop(columns=['ID', 'LABEL'])
    X = df
    X2_proj = lda.transform(X)

    plt.figure(figsize=(8,6))

    # Jeu 1 en bleu
    plt.scatter(X1_proj[:,0], X1_proj[:,1], alpha=0.7, label="Jeu 1", color="blue")

    # Jeu 2 en rouge
    plt.scatter(X2_proj[:,0], X2_proj[:,1], alpha=0.7, label="Jeu 2", color="red")

    plt.xlabel("LD1")
    plt.ylabel("LD2")
    plt.title("LDA new and old")
    plt.legend()
    plt.show()



def run_lda():
    """ """

    # load model
    lda = joblib.load("model/lda_basic.joblib")

    # load new data
    df = pd.read_csv("data/data_vst_lda.csv")
    df = df.drop(columns=["Unnamed: 0"])
    X = df
    X1_proj = lda.transform(X)

    # load precisesads data
    df = pd.read_csv("data/precisesads/lda_input.csv")
    df = df.drop(columns=['ID', 'LABEL'])
    X = df
    X2_proj = lda.transform(X)

    plt.figure(figsize=(8,6))

    # Jeu 1 en bleu
    plt.scatter(X1_proj[:,0], X1_proj[:,1], alpha=0.7, label="Jeu 1", color="blue")

    # Jeu 2 en rouge
    plt.scatter(X2_proj[:,0], X2_proj[:,1], alpha=0.7, label="Jeu 2", color="red")

    plt.xlabel("LD1")
    plt.ylabel("LD2")
    plt.title("LDA new and old")
    plt.legend()
    plt.show()


if __name__ == "__main__":


    # run_lda()
    run_lda_normalize()
