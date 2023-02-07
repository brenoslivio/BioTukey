import streamlit as st
import zipfile
import os
import shutil

st.set_page_config(page_title = "RNAseq Analysis Tool", page_icon = ":microscope:", layout="wide")

st.title("Analysis tool for RNA sequences")
st.markdown('''##### <span style="color:gray">Discover unique information about RNA sequences with Machine Learning approaches.</span>
            ''', unsafe_allow_html=True)
                
col1, col2, col3 = st.sidebar.columns([1,8,1])

with col2:
    st.markdown('# RNAseq Analysis Tool')

st.sidebar.markdown('---')

with col2:
    uploaded_file = st.sidebar.file_uploader("Upload your ZIP file with FASTA files by RNA type", type=["zip"])

tab1, tab2, tab3, tab4 = st.tabs(['General Statistics', 'Classification', 'Visualization', 'FAQ'])


if uploaded_file:

    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    else:
        os.makedirs('tmp')

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall('tmp')

    filenames = [f for f in os.listdir('tmp') if os.path.isfile(os.path.join('tmp', f))]
    classes = [os.path.splitext(f)[0] for f in filenames]

    with tab1:
        st.markdown("You have " + str(len(classes)) + " RNA types: " + ', '.join(classes) + '.')

    

with tab4:
    st.markdown("**1. How my ZIP file has to look like for me to upload it?**")

    st.markdown("Your ZIP file need to have FASTA files by RNA type with their respective names.")

    st.image('imgs/zipfile.png')
