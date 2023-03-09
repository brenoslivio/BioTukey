import streamlit as st
import utils.css_injection as css_injection
import pandas as pd
from PIL import Image
import utils, setup

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = 'imgs/biotukey_icon.png', layout="wide")

    utils.inject_css()
    
    st.sidebar.markdown("---")

    uploaded_files = st.sidebar.file_uploader("Select your FASTA files by sequence class", 
                                                accept_multiple_files=True, type=["fasta", "fa", "faa"], 
                                                help="Each file must be named according to its class (e.g., sRNA.fasta). FASTA, FA, FAA files only.")

    study_example = st.sidebar.selectbox("Or select study example", ['', "ncRNAs"])

    option = st.sidebar.radio("Select option to load", ["Manual", "Example"], horizontal=True)
                
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


    st.sidebar.markdown("---")

    if (option == "Manual" and uploaded_files) or (option == "Example" and study_example):

        page = st.sidebar.radio("Select page:", ["Home", "Overview", "Feature Engineering", "Classification"])

        if page == "Home":
            st.title("BioTukey")
            st.markdown('''##### <span style="color:gray">Comprehensive toolkit for interactive data analysis, engineering, and classification of biological sequences.</span>
                        ''', unsafe_allow_html=True)
        elif page == "Overview":
            match study_example:
                case "ncRNAs":
                    st.info("**Dataset from the following published paper:**\n \
                            Robson P Bonidia, Anderson P Avila Santos, Breno L S de Almeida, \
                            Peter F Stadler, Ulisses N da Rocha, Danilo S Sanches, \
                            André C P L F de Carvalho, BioAutoML: automated feature engineering \
                            and metalearning to predict noncoding RNAs in bacteria, \
                            Briefings in Bioinformatics, Volume 23, Issue 4, July 2022, \
                            bbac218, https://doi.org/10.1093/bib/bbac218")
            setup.overview.load(files, seq_type)
        #elif page == "Feature Engineering":

if __name__ == '__main__':
    runUI()