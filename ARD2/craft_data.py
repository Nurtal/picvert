import pandas as pd
import re
import numpy as np

def craft_precisesads_lda_data():
    """Regenarate input dataset for LDA"""

    # load data
    df = pd.read_csv("data/precisesads/vst_rnaseq_SSc_genes_transposed.csv")
    df = df.rename(columns={'Unnamed: 0':'ID'})

    # filter on authorized genes
    legal_gene = list(pd.read_csv('data/gene_intersection.csv')['SYMBOL'])
    var_to_keep = ['ID']
    for g in legal_gene:
        var_to_keep.append(g)
    df = df[var_to_keep]

    # load manifest
    id_to_label = {}
    manifest = pd.read_csv("data/manifest2.csv")
    for index, row in manifest.iterrows():
        id_to_label[row['ID']] = row['Cluster']

    # add label
    df['LABEL'] = df['ID'].replace(id_to_label)

    # clean
    df = df[df['LABEL'].isin(['clust1', 'clust2', 'clust3'])]

    # save
    df.to_csv("data/precisesads/lda_input.csv", index=False)


def compute_vst(mu, stderr, n, min_alpha=1e-8):
    # variance empirique
    var = (stderr * np.sqrt(n))**2  # = stderr^2 * n
    
    # dispersion alpha
    alpha = (var - mu) / (mu**2)
    alpha = np.where(alpha > min_alpha, alpha, min_alpha)
    
    # VST
    vst = (1.0 / np.sqrt(alpha)) * np.arcsinh(np.sqrt(alpha * mu))
    vst = np.where(alpha < min_alpha * 10, np.sqrt(mu), vst)  # approx Poisson
    return vst



def create_vst_new_file():
    """Take downloaded txt file and turn it into a vst file"""
    
    # load data
    df = pd.read_csv("data/data.txt", sep="\t")

    # compute VST
    for col in df.columns:
        if re.search("AVG_Signal", col):
            sample = col.replace(".AVG_Signal", "")
        
            mu = df[f"{sample}.AVG_Signal"]
            stderr = df[f"{sample}.BEAD_STDERR"]
            n = df[f"{sample}.Avg_NBEADS"]
        
            df[f"{sample}.VST"] = compute_vst(mu, stderr, n)

    # keep VST signal
    var_to_keep = ["SYMBOL"]
    for k in list(df.keys()):
        if re.search("VST", k):
            var_to_keep.append(k)
    df = df[var_to_keep]
    df = df.rename(columns={'SYMBOL':'ID'})
    df = df.set_index('ID')
    df = df.T

    # save
    df.to_csv("data/data_vst.csv")



def compute_intersection():

    # load data preciesads
    df = pd.read_csv("data/precisesads/lda_input.csv")

    # load new data
    df_new = pd.read_csv("data/data_vst.csv")

    # compute intersection
    intersection = []
    for x in list(df.keys()):
        if x in list(df_new.keys()):
            intersection.append(x)

    # save
    intersection_file = open("data/gene_intersection.csv", "w")
    intersection_file.write("SYMBOL\n")
    for s in intersection:
        intersection_file.write(s+"\n")
    intersection_file.close()


def craft_new_lda_data():
    """Craft lda input file from new data"""

    # load new vst
    df = pd.read_csv("data/data_vst.csv")

    # filter on authorized genes
    var_to_keep = []
    legal_gene = list(pd.read_csv('data/gene_intersection.csv')['SYMBOL'])
    for g in legal_gene:
        var_to_keep.append(g)
    df = df[var_to_keep]

    # sort accordingly to lda file use as input
    df_sort = pd.read_csv("data/precisesads/lda_input.csv")
    order = []
    for var in list(df_sort.keys()):
        if var in list(df.keys()):
            order.append(var)
    df = df[order]

    # save
    df.to_csv("data/data_vst_lda.csv")




def create_normalize_data():
    """ """

    df = pd.read_csv("data/precisesads/lda_input.csv")
    df_new = pd.read_csv("data/data_vst_lda.csv")
    df_new = df_new.rename(columns={"Unnamed: 0":"ID"})

    df['GROUP'] = "P"
    df_new['LABEL'] = "NA"
    df_new['GROUP'] = "N"
    df_all = pd.concat([df, df_new])
    df_all['ID'] = 'not used'

    # Sélectionne uniquement les colonnes numériques
    num_cols = df_all.select_dtypes(include="number").columns

    # Applique la normalisation min–max sur [0, 100] à ces colonnes
    df_all[num_cols] = (df_all[num_cols] - df_all[num_cols].min()) / (df_all[num_cols].max() - df_all[num_cols].min()) * 100

    # save precisesads
    df_a = df_all[df_all['GROUP'] == "P"]
    df_a = df_a.drop(columns=['GROUP'])
    df_a.to_csv("data/precisesads/input_normalized_lda.csv", index=False)

    # save new
    df_b = df_all[df_all['GROUP'] == "N"]
    df_b = df_b.drop(columns=['GROUP'])
    df_b.to_csv("data/data_normalized_lda.csv", index=False)


def create_normalize_data2():
    """ """

    df = pd.read_csv("data/precisesads/lda_input.csv")
    df = df.drop(columns=['ID'])
    num_cols = df.select_dtypes(include="number").columns
    df[num_cols] = (df[num_cols] - df[num_cols].min()) / (df[num_cols].max() - df[num_cols].min()) * 100
    df['ID'] = 'not used'
    df.to_csv("data/precisesads/input_normalized2_lda.csv", index=False)

    df = pd.read_csv("data/data_vst_lda.csv")
    df = df.rename(columns={"Unnamed: 0":"ID"})
    num_cols = df.select_dtypes(include="number").columns
    df[num_cols] = (df[num_cols] - df[num_cols].min()) / (df[num_cols].max() - df[num_cols].min()) * 100
    df['ID'] = 'not used'
    df.to_csv("data/data_normalized2_lda.csv", index=False)

if __name__ == "__main__":


    # compute_intersection()
    # craft_precisesads_lda_data()
    # create_vst_new_file()
    # craft_new_lda_data()
    # create_normalize_data()
    create_normalize_data2()
