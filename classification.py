import streamlit as st
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, cross_val_predict, StratifiedKFold
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


def classify(seqs):
    st.markdown("You can classify between " + str(len(seqs)) + " RNA types: " + ', '.join(seqs) + '. ')

    col1, col2 = st.columns(2)

    df_show = pd.DataFrame()

    with col1:
        with st.form("classify_form"):
            features = st.multiselect('Select features you want for classification:', ["",
                                                                                        "Nucleotide Composition (k-mer = 1)", 
                                                                                        "Dinucleotide Composition (k-mer = 2)", 
                                                                                        "Trinucleotide Composition (k-mer = 3)", 
                                                                                        "Tetranucleotide Composition (k-mer = 4)", 
                                                                                        "Pentanucleotide Composition (k-mer = 5)", 
                                                                                        "kGap", "ORF", "Fickett Score"])

            submitted = st.form_submit_button("Submit")

        if submitted and features:
            df = feature_extraction(seqs, features)
            df_show = df.copy().set_index('type')

            df = df.drop('nameseq', axis = 1)

            y = df.pop('type')

            X = df.values
    
    with col2:
        if features:
            st.markdown("##### Dataset with concatenated features")
            st.dataframe(df_show)
        
    ml = st.selectbox("Select a classifier:", ["", "Random Forest (RF)", "Multilayer Perceptron (MLP)", "Support Vector Machines (SVM)"])

    if ml and features:
        clf = RandomForestClassifier()

        if ml == "Random Forest (RF)":
            st.markdown("Select the number of trees in the random forest.")
        elif ml == "Multilayer Perceptron (MLP)":
            st.markdown("MLP")
        elif ml == "Support Vector Machines (SVM)":
            C = float(st.text_input("C", 10))
            gamma = float(st.text_input("gamma", 0.01))
            kernel = st.text_input("kernel", "rbf")

            clf = SVC(C = C, gamma = gamma, kernel = kernel, random_state = 42)

        steps = [('scaler', StandardScaler()), ('model', clf)]

        model = Pipeline(steps = steps)

        avg_score = cross_val_score(model, X, y, cv = 10, scoring = 'accuracy', n_jobs = -1).mean()

        st.markdown('ACC:' + str(avg_score))


    #st.checkbox("Dataset for testing the model")