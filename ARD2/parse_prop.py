import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def spot_high_range(n):
    """ """

    # load data
    dfn = pd.read_csv("data/data_new_properties.csv")
    dfp = pd.read_csv("data/precisesads/data_precisesads_properties.csv")

    dfn = dfn.nlargest(n, 'RANGE')
    dfp = dfp[dfp['SYMBOL'].isin(list(dfn['SYMBOL']))]

    sn_to_range = {}
    for index, row in dfn.iterrows():
        sn_to_range[row['SYMBOL']] = row['RANGE']
    sp_to_range = {}
    for index, row in dfp.iterrows():
        sp_to_range[row['SYMBOL']] = row['RANGE']

    keys = list(sn_to_range.keys())
    vn = []
    vp = []
    for k in keys:
        vn.append(sn_to_range[k])
        vp.append(sp_to_range[k])
    
    # save keys
    key_file = open("data/gene_to_drop.csv", "w")
    key_file.write("SYMBOL\n")
    for k in keys:
        key_file.write(f"{k}\n")
    key_file.close()
    
    # Position des barres
    x = np.arange(len(keys))  
    width = 0.35  # largeur des barres

    # Création du graphique
    fig, ax = plt.subplots(figsize=(8,6))
    ax.bar(x - width/2, vn, width, label="New")
    ax.bar(x + width/2, vp, width, label="PREC")

    # Personnalisation
    ax.set_xticks(x)
    ax.set_xticklabels(keys)
    ax.set_ylabel("Valeurs")
    ax.set_title("Comparaison")
    ax.legend()

    plt.show()


if __name__ == "__main__":

    spot_high_range(10)
