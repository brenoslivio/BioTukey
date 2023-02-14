import streamlit as st
import zipfile
import os
import shutil
import pandas as pd
import stats, classification, faq, seqdata
import subprocess

def process_files(uploaded_file):
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    else:
        os.makedirs('tmp')

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall('tmp')

    files = {os.path.splitext(f)[0] : 'tmp/' + f for f in os.listdir('tmp') if os.path.isfile(os.path.join('tmp', f))}

    for type in files:
        subprocess.run(['python', 'MathFeature/preprocessing/preprocessing.py', '-i', files[type], '-o', 'tmp/pre_' + type + '.fasta'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        files[type] = 'tmp/pre_' + type + '.fasta'

    return files

def runUI():
    st.set_page_config(page_title = "RNAseq Analysis Tool", page_icon = ":microscope:", layout="wide")

    # CSS to inject contained in a string
    hide_table_row_index = """
                <style>
                thead tr th:first-child {display:none}
                tbody th {display:none}
                </style>
                """

    # Inject CSS with Markdown
    st.markdown(hide_table_row_index, unsafe_allow_html=True)

    st.title("Analysis tool for RNA sequences")
    st.markdown('''##### <span style="color:gray">Discover unique information about RNA sequences with an overall interactive environment and Machine Learning approaches.</span>
                ''', unsafe_allow_html=True)
                    
    col1, col2, col3 = st.sidebar.columns([1,8,1])

    with col2:
        st.markdown('# RNAseq Analysis Tool')

    st.sidebar.markdown('---')

    uploaded_file = st.sidebar.file_uploader("Upload your ZIP file with FASTA files by RNA type", type=["zip"])

    st.sidebar.markdown('---')

    tab1, tab2, tab3, tab4 = st.tabs(['General Statistics', 'Classification', 'Feature Visualization', 'FAQ'])

    if uploaded_file:

        files = process_files(uploaded_file)

        seqs = {}

        for type in files:
            seq = seqdata.Seq(files[type], type)
            seqs[type] = seq
    
        with tab1:
            stats.general_stats(seqs)

        with tab2:
            classification.classify(seqs)

    with tab4:
        faq.help()

if __name__ == '__main__':
    runUI()