import streamlit as st
import zipfile
import os
import shutil
import stats
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import collections

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

st.sidebar.markdown('TODO:')

st.sidebar.markdown('- Histograms exploring nucleotides frequency')

st.sidebar.markdown('- Violinplots for nucleotides frequency')

st.sidebar.markdown('- Hypothesis tests between n-mer averages')


tab1, tab2, tab3, tab4 = st.tabs(['General Statistics', 'Classification', 'Feature Visualization', 'FAQ'])


if uploaded_file:

    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    else:
        os.makedirs('tmp')

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall('tmp')

    files = {os.path.splitext(f)[0] : f for f in os.listdir('tmp') if os.path.isfile(os.path.join('tmp', f))}
 
    with tab1:
        st.markdown("You have " + str(len(files.keys())) + " RNA types: " + ', '.join(files.keys()) + '. ' \
                    + 'We have the following summary statistics for the sequences:')

        df = pd.DataFrame()
        seqs = {}

        for type in files.keys():

            seq_df = stats.SeqData('tmp/' + files[type])
            seqs[type] = seq_df

            seq_desc = seq_df.desc()

            stats_df = pd.DataFrame({"type": type, 
                                "num_seqs": seq_desc['count'], 
                                "min_len (bp)": seq_desc['min'],
                                "avg_len (bp)": seq_desc['mean'],  
                                "max_len (bp)": seq_desc['max'],
                                "std_len (bp)": seq_desc['std'],
                                "Q1 (bp)": seq_desc['25%'],
                                "Q2 (bp)": seq_desc['50%'],
                                "Q3 (bp)": seq_desc['75%'],
                                "gc_content (%)": seq_df.gc_content()}, 
                                index = [0])
            
            df = pd.concat([df, stats_df]).reset_index(drop=True)

        st.table(df.style.format({col: "{:.2f}" for col in df.columns if col != 'type'}))

        tab1_tab1, tab1_tab2 = st.tabs(['Distribution', 'Statistical Significance'])

        with tab1_tab1:
            
            df_plot = pd.DataFrame()

            for type in seqs.keys():
                new_df = pd.DataFrame(seqs[type].df['seq'])
                new_df['type'] = type
                
                df_plot = pd.concat([df_plot, new_df]).reset_index(drop=True)

            figures = {}

            for N in ['A', 'G', 'C', 'T']:
                df_plot[N] = df_plot['seq'].apply(lambda x : x.count(N) / len(x))
                figures[N] = px.box(df_plot, x='type', y=N, color='type')

            figures_traces = collections.defaultdict(list)

            for N in figures.keys():
                for trace in range(len(figures[N]["data"])):
                    figures_traces[N].append(figures[N]["data"][trace])

            fig = make_subplots(rows=2, cols=2,
                subplot_titles = ['Adenine', 'Guanine', 'Cytosine', 'Thymine'])

            for traces in figures_traces['A']:
                fig.append_trace(traces, row=1, col=1)

            for traces in figures_traces['G']:
                fig.append_trace(traces, row=1, col=2)
            
            for traces in figures_traces['C']:
                fig.append_trace(traces, row=2, col=1)

            for traces in figures_traces['T']:
                fig.append_trace(traces, row=2, col=2)

            fig.update_layout(height=800, width=1000, showlegend=False, title_text="Boxplot for nucleotide ratio by RNA type")

            st.plotly_chart(fig)

with tab4:
    st.markdown("**1. How my ZIP file has to look like for me to upload it?**")

    st.markdown("Your ZIP file need to have FASTA files by RNA type with their respective names.")

    st.image('imgs/zipfile.png')
