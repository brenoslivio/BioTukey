import streamlit as st
import zipfile
import os
import shutil
import stats
import pandas as pd

st.set_page_config(page_title = "RNAseq Analysis Tool", page_icon = ":microscope:", layout="wide")

st.title("Analysis tool for RNA sequences")
st.markdown('''##### <span style="color:gray">Discover unique information about RNA sequences with an overall interactive environment and Machine Learning approaches.</span>
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
        st.markdown("You have " + str(len(classes)) + " RNA types: " + ', '.join(classes) + '. ' \
                    + 'We have the following summary statistics for the sequences:')

        df = pd.DataFrame()

        for i in range(len(classes)):

            seq_desc = stats.SeqData('tmp/' + filenames[i]).desc()

            stats_df = pd.DataFrame({"type": classes[i], 
                                "num_seqs": seq_desc['count'], 
                                "min_len": seq_desc['min'],
                                "avg_len": seq_desc['mean'],  
                                "max_len": seq_desc['max'],
                                "std_len": seq_desc['std'],
                                "25%": seq_desc['25%'],
                                "50%": seq_desc['50%'],
                                "75%": seq_desc['75%'],}, 
                                index = [0])
            
            df = pd.concat([df, stats_df]).reset_index(drop=True)

        st.table(df)
    

with tab4:
    st.markdown("**1. How my ZIP file has to look like for me to upload it?**")

    st.markdown("Your ZIP file need to have FASTA files by RNA type with their respective names.")

    st.image('imgs/zipfile.png')
