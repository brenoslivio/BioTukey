import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
import utils, setup
import os, shutil

def clear_session():
    if 'test' in st.session_state:
        del st.session_state['test']

def runUI():
    st.set_page_config(page_title = "BioTukey", initial_sidebar_state = "expanded", page_icon = 'imgs/biotukey_icon.png', layout="wide")
    
    utils.inject_css()

    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey')

    if 'biotukey_dir' not in st.session_state:
        try:
            shutil.rmtree(dir_path)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
            print('Creating Directory...')
    
        os.makedirs(dir_path, exist_ok=True)
        os.makedirs(os.path.join(dir_path, 'train'), exist_ok=True)
        os.makedirs(os.path.join(dir_path, 'test'), exist_ok=True)
        os.makedirs(os.path.join(dir_path, 'feat_engineering/train'), exist_ok=True)
        os.makedirs(os.path.join(dir_path, 'feat_engineering/test'), exist_ok=True)

        st.session_state['biotukey_dir'] = True

    for _ in range(10):
        st.sidebar.markdown("")
    st.sidebar.markdown("---")

    uploaded_files = st.sidebar.file_uploader("Select your FASTA files by sequence class",
                                                accept_multiple_files=True, type=["fasta", "fa", "faa"], 
                                                help="Each file must be named according to its class (e.g., sRNA.fasta). FASTA, FA, FAA files only.")
    type_selection = st.sidebar.selectbox("Select sequence type:", ["DNA/RNA", "Protein"])
    study_example = st.sidebar.selectbox("Or select study example", ['', "ncRNAs", "secreted proteins"])

    option = st.sidebar.radio("Select option to load", ["Manual", "Example"], horizontal=True)

    page = option_menu(None, ["Home", "Overview", "Feature Engineering", "Classification"], 
    icons=['house', 'search', "gear", 'diagram-2'], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    if page == "Home":
        st.title("BioTukey")
        st.markdown('''##### <span style="color:gray">Comprehensive toolkit for interactive data analysis, engineering, and classification of biological sequences.</span>
                    ''', unsafe_allow_html=True)
        st.markdown("""BioTukey was developed using Python programming, with the Streamlit framework as the foundation for its GUI. Streamlit is a web application framework that simplifies creating interactive data-driven applications.
                    In this context, BioTukey leverages the capabilities of Streamlit to function as a web application, enabling users to interact with the application through its GUI. The application is designed as a responsive 
                    single-page application (SPA) to provide users with a user experience akin to a desktop application. It is organized into three main modules/pages: Overview, Feature Engineering, and Classification. 
                    Each module offers distinct functionalities to support analyzing and exploring biological sequences.""")

        st.markdown("""**Overview** module provides comprehensive data visualization and summary statistics. It presents users with diverse statistics on the provided sequences, including sequence alignment visualization, k-mer distribution, 
                    and structure visualization for DNA/RNA and protein sequences. This module aims to provide users with a holistic understanding of the characteristics and properties of the input sequences.""")
        
        st.markdown("""**Feature Engineering** module empowers users to extract different types of descriptors from the sequences, generating feature vectors that can be concatenated. These feature vectors serve as input for visualizations 
                    such as distribution plots, dimensionality reduction plots using techniques like PCA, t-SNE, or UMAP, correlation heatmaps, and feature importance bar charts. By allowing users to manipulate and explore the extracted 
                    features, this module facilitates a deeper understanding of the data and enables informed decision-making in subsequent analysis steps.""")

        st.markdown("""**Classification** module allows users to perform classification tasks on the provided sequences. Users can select between two popular machine learning models, Random Forest, and XGBoost, as the algorithms for classification. 
                    They can further choose to apply 10-fold cross-validation, either with or without a separate test set. Additionally, users have the option to engage with the model manually by creating a custom pipeline. This includes 
                    selecting hyperparameters, filtering for the most important features, and applying dimensionality reduction techniques. The module also provides performance scores of the model, such as accuracy, precision, recall, and F1-score, along with a confusion matrix to evaluate the model's predictive capabilities.""")
        
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
                    clear_session()
                    st.sidebar.success("Files submitted with success.")
                    files, seq_type = utils.processing.process_files(uploaded_files, type_selection, True)
            case "Example":
                if not study_example:
                    st.sidebar.warning("Please select study example.")
                else:
                    clear_session()
                    st.sidebar.success("Example submitted with success.")
                    files, seq_type = utils.processing.load_study(study_example, True)

    if (option == "Manual" and uploaded_files and files) or (option == "Example" and study_example and files):
        if page == "Overview":
            if option == "Example":
                match study_example:
                    case "ncRNAs":
                        st.info("**Data set from the following published paper:** \
                                Robson P Bonidia, Anderson P Avila Santos, Breno L S de Almeida, \
                                Peter F Stadler, Ulisses N da Rocha, Danilo S Sanches, \
                                André C P L F de Carvalho, BioAutoML: automated feature engineering \
                                and metalearning to predict noncoding RNAs in bacteria, \
                                Briefings in Bioinformatics, Volume 23, Issue 4, July 2022, \
                                bbac218, https://doi.org/10.1093/bib/bbac218")
                    case "secreted proteins":
                        st.info("**Data set from the following published paper:** \
                                    Yanju Zhang, Sha Yu, Ruopeng Xie, Jiahui Li, André Leier, Tatiana T Marquez-Lago, \
                                Tatsuya Akutsu, A Ian Smith, Zongyuan Ge, Jiawei Wang, Trevor Lithgow, Jiangning Song, \
                                PeNGaRoo, a combined gradient boosting and ensemble learning framework for predicting \
                                non-classical secreted proteins, Bioinformatics, Volume 36, \
                                Issue 3, February 2020, Pages 704–712, \
                                https://doi.org/10.1093/bioinformatics/btz629")
            setup.overview.load(files, seq_type)
        elif page == "Feature Engineering":
            setup.feature_engineering.load(files, seq_type)
        elif page == "Classification":
            setup.classification.load(files, seq_type, option, study_example)

if __name__ == '__main__':
    runUI()