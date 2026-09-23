
# Polyreativity Prediction baased on sequences embedding data
# Embedding method: Ablang2 with  VH and VL sequences
# ML method: Deep learning NN, DL RNN, Random Forest, XGBoost, SVM
# ML trainningset: IPI data and PeterTessierLab Dataset1
# Input data:  excel file includes cdr3_aa,  IPI vh_scaffold and vl_scaffold        
# Output: five machine learning prediction scores


# python library import 
import tensorflow as tf
import keras
import scipy 
from interpolation import interp
from itertools import chain
from itertools import product as iterproduct
import warnings
import pickle
from sklearn.preprocessing import LabelEncoder
import joblib
from pathlib import Path

import pandas as pd
import numpy as np
import ablang2
np.random.seed(42) # fix random seed for reproducibility

import sys
import os
import glob
from anarci import anarci
from abnumber import Chain
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from pandarallel import pandarallel
pandarallel.initialize(25)


# PIPELINE CONFIGUATION HERE #
pipeline_path="/Users/Hoan.Nguyen/ComBio/MachineLearning/SEC"
sys.path.append(pipeline_path)
os.chdir(pipeline_path)

# IPI antibody library info 
constant="ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
fr4="WGQGTLVTVSS"
IPI_LIB="IPI_VLVH_LIB_ALL.csv"





def load_previous_db(db_path: Path) -> tuple[set, set, set, dict, dict, dict, dict]:
    """Load previous antibodies and create lookup sets + dicts for BARCODEs"""
    if not db_path.exists():
        print("Previous antibodies DB not found — skipping repeat flags")
        empty_set = set()
        empty_dict = {}
        return empty_set, empty_set, empty_set, empty_dict, empty_dict, empty_dict

    df = pd.read_excel(db_path)
    
    # Assume columns: "BARCODE", "CDR3", "heavy", "light" (adjust if different)
    df = df[["BARCODE", "CDR3", "heavy", "light","antigen"]].dropna()
    #df = df[["CDR3", "heavy", "light"]].dropna()
    df['CDR3'] = df['CDR3'].str[1:]
    df['heavy'] = df['heavy'].str[1:]
    df['light'] = df['light'].str[1:]

    #df['CDR3'] = df['CDR3'].str[1:] if df['CDR3'].str.startswith('C').all() else df['CDR3']  # optional strip "C"
    #df['heavy'] = df['heavy'].str[1:] if df['heavy'].str.startswith('V').all() else df['heavy']
    #df['light'] = df['light'].str[1:] if df['light'].str.startswith('V').all() else df['light']

    # Sets for boolean flags
    cdr3_set = set(df["CDR3"])
    vh_set = set(df["heavy"] + "|" + df["CDR3"])
    ab_set = set(df["heavy"] + "|" + df["light"] + "|" + df["CDR3"])

    df.rename(columns={"CDR3": "cdr3_aa", "heavy": "vh_scaffold", "light": "vl_scaffold"}, inplace=True)
    cdr3_set = set(df["cdr3_aa"])
    vh_set = set(df["vh_scaffold"] + "|" + df["cdr3_aa"])
    ab_set = set(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])

    # Dicts for BARCODE lists (sequence → ";" joined BARCODEs)
    cdr3_dict = df.groupby("cdr3_aa")["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    vh_dict = df.groupby(df["vh_scaffold"] + "|" + df["cdr3_aa"])["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    ab_dict = df.groupby(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    ab_dict_ag = df.groupby(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])["antigen"].apply(lambda x: ";".join(x.astype(str))).to_dict()

    print(f"Loaded {len(df)} previous antibodies for repeat checking (with BARCODEs)")
    return cdr3_set, vh_set, ab_set, cdr3_dict, vh_dict, ab_dict,ab_dict_ag


# Generating the full AA Sequence from IPI NGS Report or Miseq report 
def Generate_FullVHVL(NGS_LEADS,heavy_scaffold='vh_scaffold',light_scaffold='vl_scaffold',cdr3='cdr3_aa', liabilities=False):
    IPI_VHVLLIB=pd.read_csv(IPI_LIB)
    NGS_LEADS['HSEQ']=''
    NGS_LEADS['LSEQ']=''
    NGS_LEADS['BARCODE']= NGS_LEADS.index
    for i in range(0,len(NGS_LEADS)):
        CDR3=NGS_LEADS.loc[i][cdr3]
        if (CDR3.startswith("CAR")):
            CDR3=CDR3[1:]
        Heavy=NGS_LEADS.loc[i][heavy_scaffold]
        Light=NGS_LEADS.loc[i][light_scaffold]
        #NGS_LEADS['BARCODE'][i]=NGS_LEADS.loc[i]['target']+'|'+NGS_LEADS.loc[i][heavy_scaffold]+'|'+NGS_LEADS.loc[i][light_scaffold]+'|'+NGS_LEADS.loc[i][cdr3]
        # ADD H and L into vl/vh_scaffold
        if (not Heavy.startswith("V")):
            Heavy=''.join(['V',Heavy])
        if (not Light.startswith("V")):
            Light=''.join(['V',Light])

        VH=IPI_VHVLLIB[IPI_VHVLLIB['Name']==Heavy]
        VL=IPI_VHVLLIB[IPI_VHVLLIB['Name']==Light]
        if (len(VH)>0):
            NGS_LEADS['HSEQ'][i]=''.join((VH['AA']+CDR3+fr4))
        else:
            NGS_LEADS['HSEQ'][i]= CDR3
        if (len(VL)>0):
            NGS_LEADS['LSEQ'][i]=''.join((VL['AA']))
    if (liabilities):
        import liabilities
        NGS_LEADS=liabilities.annotate_liabilities(NGS_LEADS,cdr3)
    return NGS_LEADS



def PSR_MULTIPLE_PRED(INPUT_LEADS,heavy_scaffold='vh_scaffold',light_scaffold='vl_scaffold',liab=False,cdr3='cdr3_aa',lang='ablang'):
    #OUTPUT_NAME=INPUT_LEADS+".PSR_Prediction.xlsx"
    OUTPUT_NAME=INPUT_LEADS
    ## read excel and csv file
    if INPUT_LEADS.endswith('.xlsx'):
         LEADS= pd.read_excel(INPUT_LEADS)
    elif INPUT_LEADS.endswith('.csv'):      
         LEADS=pd.read_csv(INPUT_LEADS)




    LEADS=Generate_FullVHVL(LEADS,heavy_scaffold=heavy_scaffold,light_scaffold=light_scaffold,cdr3=cdr3,liabilities=liab)
    
    # remove vl_scaffold and vh_scaffold with UNK
    LEADS=LEADS[~LEADS[heavy_scaffold].str.contains('UNK')]
    LEADS=LEADS[~LEADS[light_scaffold].str.contains('UNK')]


    # remove BARCODE column
    if 'BARCODE' in LEADS.columns:
        LEADS=LEADS.drop(columns=['BARCODE'])
    
    # remove cdr3_functional as FALSE
    if 'cdr3_functional' in LEADS.columns:
        LEADS=LEADS[LEADS['cdr3_functional']==True]

    LEADS=LEADS[["cdr3_aa","vh_scaffold","vl_scaffold","count","freq","HSEQ","LSEQ"]]

    print(f"Total unique antibodies before collapsing: {len(LEADS)}")
    
    # now, LEADS dataframe has columns: cdr3_aa, vh_scaffold, vl_scaffold, count, freq, HSEQ, LSEQ
    # collapse dataset base on cdr3_aa, vh_scaffold, vl_scaffold and sum on count, and recalate freq=sum/total_count
    total_count = LEADS['count'].sum()
    LEADS = LEADS.groupby(['cdr3_aa', 'vh_scaffold', 'vl_scaffold','HSEQ','LSEQ'], as_index=False).agg({'count': 'sum'})
    LEADS['freq'] = LEADS['count'] / total_count
    LEADS = LEADS.sort_values(by='count', ascending=False).reset_index(drop=True)
    print(f"Total unique antibodies after collapsing: {len(LEADS)}")
    
    # Load previous DB for repeat flags

    cdr3_prev, vh_prev, ab_prev, cdr3_barcode_dict, vh_barcode_dict, ab_barcode_dict,ab_ag_dict = load_previous_db(Path("/Users/Hoan.Nguyen/ComBio/IPIAbDiscov/data/All_mAb_20251216_FACS_BLI.xlsx"))

    LEADS["TAB_ID"] = (LEADS["vh_scaffold"] + "|" + LEADS["vl_scaffold"] + "|" + LEADS['cdr3_aa']).map(ab_barcode_dict).fillna("")



    if INPUT_LEADS.endswith('.xlsx'):
         LEADS.to_excel(OUTPUT_NAME ,index=False)
    elif INPUT_LEADS.endswith('.csv'):    
         LEADS.to_csv(OUTPUT_NAME ,index=False) 

    return LEADS



def get_full_anarci_anno(df: pd.DataFrame) -> pd.DataFrame:
    def get_ANARCI(seq: str):
        if len(seq) <= 50:
            return "Not Fully Annotated||||||"
        try:
            chain = Chain(seq, scheme="imgt")
            return f"{chain}|{chain.fr1_seq}|{chain.cdr1_seq}|{chain.fr2_seq}|{chain.cdr2_seq}|{chain.fr3_seq}|{chain.cdr3_seq}"
        except Exception:
            return "Not Fully Annotated||||||"

    df["anarci"] = df["aa"].parallel_apply(get_ANARCI)
    df[["CHAIN", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3"]] = df["anarci"].str.split("|", expand=True)
    return df


def NT2AA(seq: str) -> str:
    candidates = [Seq(seq[i:]).translate(to_stop=True) for i in range(3)]
    candidates += [Seq(str(Seq(seq).reverse_complement())[i:]).translate(to_stop=True) for i in range(3)]
    return max(candidates, key=len) if candidates else ""

def run_anarci():
    df=pd.read_excel("/Users/Hoan.Nguyen/ComBio/MachineLearning/IPIPred/data/kmab_a_2593055_sm5212_sequence.xlsx")
    df["aa"] = df["HSEQ"]
    df = get_full_anarci_anno(df)
    df.to_excel("/Users/Hoan.Nguyen/ComBio/MachineLearning/IPIPred/data/kmab_a_2593055_sm5212_sequence.xlsx",index=False)
        


if __name__ == "__main__":
    LEADS=pd.read_excel("/Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq112/leads.xlsx")
    LEADS=Generate_FullVHVL(LEADS,heavy_scaffold="vh_scaffold",light_scaffold="vl_scaffold",cdr3="cdr3_aa")
    LEADS.to_excel("/Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq112/leads.xlsx",index=False)
    #run_anarci()
    #sys.exit(0)
     

    # DATA FOLDER INIPUT HERE
    FASTQ_FOLDER="/Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq112/"
    lang='antiberta2-cssp'
    column_names = ['target', 'pass', 'fail']
    psr_table = pd.DataFrame(columns=column_names)
    for INPUT_LEADS in glob.glob(FASTQ_FOLDER+ '*.xlsx'):
        print(INPUT_LEADS)
        LEADS=PSR_MULTIPLE_PRED(INPUT_LEADS,liab=False,cdr3='cdr3_aa',lang=lang)

 





