import os
import urllib.request


def dl_data():
    """DL data to parse"""

    target_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130953&format=file&file=GSE130953%5FNon%2Dnormalized%5Fdata%2Etxt%2Egz"

    # create sub data directory if not exist
    if not os.path.isdir("data"):
        os.mkdir("data")

    # run download
    urllib.request.urlretrieve(target_url, "data/data.txt.gz")
    

    


if __name__ == "__main__":

    dl_data()
