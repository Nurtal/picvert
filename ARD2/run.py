import craft_data
import parse_prop
import run_lda
import os


def run(n):
    """Zbleh"""

    pfile = "data/precisesads/lda_input.csv"
    nfile = "data/data_vst_lda.csv"

    # extrac prop
    craft_data.extract_properties()    

    # pick gene to delete
    parse_prop.spot_high_range(n)

    # craft reduced dataset
    craft_data.craft_reduce_dataset()

    # normalize dataset
    craft_data.create_normalize_data2()

    # train lda
    os.system("python /home/drfox/workspace/SPELLBOOK/spellbook.py lda /home/drfox/workspace/picvert/ARD2/data/precisesads/input_normalized2_lda.csv /home/drfox/workspace/picvert/ARD2/data/lda_auto")
    os.system("cp /home/drfox/workspace/picvert/ARD2/data/lda_auto/model.joblib /home/drfox/workspace/picvert/ARD2/model/lda_auto.joblib")

    # run lda
    run_lda.run_lda_normalize()
    


    


if __name__ == "__main__":


    run(150)
