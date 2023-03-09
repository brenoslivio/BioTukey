import streamlit as st
import utils
import pandas as pd
import re

def kmer_feature(seq, feature):
    k = int(re.findall(r'\d+', feature)[0])

    _, kmer_df = seq.kmer_count(k)

    return kmer_df

def feature_extraction(seqs, features):

    df_feat = pd.DataFrame()

    for feature in features:
        df = pd.DataFrame()
        for type in seqs:

            if 'k-mer' in feature:
                new_df = kmer_feature(seqs[type], feature)

            df = pd.concat([df, new_df]).reset_index(drop=True)
        
        df_feat = pd.concat([df_feat, df], axis = 1).T.drop_duplicates().T

    return df_feat

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = 'imgs/biotukey_icon.png', layout="wide")

    utils.inject_css()

if __name__ == '__main__':
    runUI()