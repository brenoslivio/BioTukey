import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
import utils, setup
import os, shutil

def runUI():
    st.set_page_config(page_title = "BioTukey", initial_sidebar_state = "expanded", page_icon = 'imgs/biotukey_icon.png', layout="wide")

    utils.inject_css()

    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey')

    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

    os.makedirs(dir_path)

    for i in range(10):
        st.sidebar.markdown("")
    st.sidebar.markdown("---")

    uploaded_files = st.sidebar.file_uploader("Select your FASTA files by sequence class",
                                                accept_multiple_files=True, type=["fasta", "fa", "faa"], 
                                                help="Each file must be named according to its class (e.g., sRNA.fasta). FASTA, FA, FAA files only.")

    study_example = st.sidebar.selectbox("Or select study example", ['', "ncRNAs", "secreted proteins"])

    option = st.sidebar.radio("Select option to load", ["Manual", "Example"], horizontal=True)

    page = option_menu(None, ["Home", "Overview", "Feature Engineering", "Classification"], 
    icons=['house', 'search', "gear", 'diagram-2'], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    if page == "Home":
        st.title("BioTukey")
        st.markdown('''##### <span style="color:gray">Comprehensive toolkit for interactive data analysis, engineering, and classification of biological sequences.</span>
                    ''', unsafe_allow_html=True)
    elif page == "Overview" and not (uploaded_files or study_example):
        st.warning("Select files or study example for Overview.")
    elif page == "Feature Engineering" and not (uploaded_files or study_example):
        st.warning("Select files or study example for Feature Engineering.")
    elif page == "Classification" and not (uploaded_files or study_example):
        st.warning("Select files or study example for Classification.")

    match option:
            case "Manual":
                if not uploaded_files:
                    st.sidebar.warning("Please select files.")
                else:
                    st.sidebar.success("Files submitted with success.")
                    files, seq_type = utils.processing.process_files(uploaded_files)
            case "Example":
                if not study_example:
                    st.sidebar.warning("Please select study example.")
                else:
                    st.sidebar.success("Example submitted with success.")
                    files, seq_type = utils.processing.load_study(study_example)

    if (option == "Manual" and uploaded_files) or (option == "Example" and study_example):

        if page == "Overview":
            if option == "Example":
                match study_example:
                    case "ncRNAs":
                        st.info("**Dataset from the following published paper:** \
                                Robson P Bonidia, Anderson P Avila Santos, Breno L S de Almeida, \
                                Peter F Stadler, Ulisses N da Rocha, Danilo S Sanches, \
                                André C P L F de Carvalho, BioAutoML: automated feature engineering \
                                and metalearning to predict noncoding RNAs in bacteria, \
                                Briefings in Bioinformatics, Volume 23, Issue 4, July 2022, \
                                bbac218, https://doi.org/10.1093/bib/bbac218")
                    case "secreted proteins":
                        st.info("**Dataset from the following published paper:** \
                                 Yanju Zhang, Sha Yu, Ruopeng Xie, Jiahui Li, André Leier, Tatiana T Marquez-Lago, \
                                Tatsuya Akutsu, A Ian Smith, Zongyuan Ge, Jiawei Wang, Trevor Lithgow, Jiangning Song, \
                                PeNGaRoo, a combined gradient boosting and ensemble learning framework for predicting \
                                non-classical secreted proteins, Bioinformatics, Volume 36, \
                                Issue 3, February 2020, Pages 704–712, \
                                https://doi.org/10.1093/bioinformatics/btz629")
                        
            setup.overview.load(files, seq_type)
            
        #elif page == "Feature Engineering":

if __name__ == '__main__':
    runUI()