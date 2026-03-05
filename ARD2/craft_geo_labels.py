"""Absolute shit show"""

import pandas as pd

def generate_visit_label():
    """ """

    # load data
    df = pd.read_csv("data/umap/geo.csv")
    dfl = pd.read_csv("data/GEO/labels_formated.csv")

    # assign label
    df['LABEL'] = dfl['VISIT']

    # save
    df.to_csv("data/umap/geo_visit.csv", index=False)

def generate_treatment_label():
    """ """

    # load data
    df = pd.read_csv("data/umap/geo.csv")
    dfl = pd.read_csv("data/GEO/labels_formated.csv")

    # assign label
    df['LABEL'] = dfl['TREATMENT']

    # save
    df.to_csv("data/umap/geo_treatment.csv", index=False)

if __name__ == "__main__":

    generate_visit_label()
    generate_treatment_label()
    
