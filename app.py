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
from scipy.stats import ttest_ind
import statsmodels.stats.proportion as prop

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

            seq = stats.SeqData('tmp/' + files[type])
            seqs[type] = seq

            seq_desc = seq.desc()

            stats_df = pd.DataFrame({"type": type, 
                                "num_seqs": seq_desc['count'], 
                                "min_len (bp)": seq_desc['min'],
                                "avg_len (bp)": seq_desc['mean'],  
                                "max_len (bp)": seq_desc['max'],
                                "std_len (bp)": seq_desc['std'],
                                "Q1 (bp)": seq_desc['25%'],
                                "Q2 (bp)": seq_desc['50%'],
                                "Q3 (bp)": seq_desc['75%'],
                                "gc_content (%)": seq.gc_content()[0]}, 
                                index = [0])

            df = pd.concat([df, stats_df]).reset_index(drop=True)

        th_props = [
            ('text-align', 'center'),
            ('font-weight', 'bold')
            ]
        
        styles = [
            dict(selector="th", props=th_props)
        ]

        df = df.style.format({col: "{:.2f}" for col in df.columns if col != 'type'}).set_table_styles(styles)

        st.table(df)

        st.markdown('---')
            
        col1, col2 = st.columns(2)

        with col1:
            st.markdown('### Nucleotide distribution insights')

            df_plot = pd.DataFrame()

            for type in seqs.keys():
                new_df = pd.DataFrame(seqs[type].df['seq'])
                new_df['type'] = type
                
                df_plot = pd.concat([df_plot, new_df]).reset_index(drop=True)

            figures = {}

            for N in ['A', 'G', 'C', 'T']:
                df_plot[N] = df_plot['seq'].apply(lambda x : x.count(N) / len(x))
                figures[N] = px.box(df_plot, x='type', y=N, color='type', color_discrete_sequence=px.colors.cyclical.Twilight)

            figures_traces = collections.defaultdict(list)

            for N in figures.keys():
                for trace in range(len(figures[N]["data"])):
                    figures_traces[N].append(figures[N]["data"][trace])

            fig = make_subplots(rows=2, cols=2,
                subplot_titles = ['Adenine', 'Guanine', 'Cytosine', 'Thymine'])

            for i, N in enumerate(['A', 'G', 'C', 'T']):
                for traces in figures_traces[N]:
                    fig.append_trace(traces, row=(i//2) + 1, col=(i%2) + 1)

            fig.update_layout(height=500, width=800, showlegend=False, title_text="Boxplot for nucleotide proportion by RNA type")

            st.plotly_chart(fig)

            st.markdown('###### Hypothesis Testing')

            st.markdown('Compare two sequence types to check for statistical significance in proportions related to GC% Content and nucleotides.')

            st.markdown('Null Hypothesis is:')

            st.markdown('$$H_0: p_1 = p_2$$')

            st.markdown('You can select the following Alternative Hypotheses:')

            st.markdown('$$H_1: p_1 \\neq p_2;\; p_1 > p_2;\; p_1 < p_2$$')

            st.markdown('Select the types of sequences to be compared:')

            type1 = st.selectbox('1st type:', [''] + list(seqs.keys()))

            type2 = st.selectbox('2nd type:', [''] + list(seqs.keys()))

            st.markdown('Select the alternative hypothesis:')

            alternative = st.selectbox('Alternative Hypothesis:', ['two-sided', 'larger', 'smaller'])

            if type1 and type2:
                count1, nobs1 = seqs[type1].gc_content()[1], seqs[type1].gc_content()[2]
                count2, nobs2 = seqs[type2].gc_content()[1], seqs[type2].gc_content()[2]

                _, pvalue = prop.test_proportions_2indep(count1, nobs1, count2, nobs2, compare='diff', alternative=alternative)
            
                st.markdown(pvalue)

            st.markdown('LET USER SELECT OPTIONS')


        with col2:
            st.markdown('### K-mers distribution')

            #https://immunarch.com/articles/web_only/v9_kmers.html
            st.markdown('TODO: Let user select which k-mer and if we want the average or sum considering all sequences')
            st.markdown('Color by the type of the sequence')

with tab4:
    st.markdown("**1. How my ZIP file has to look like for me to upload it?**")

    st.markdown("Your ZIP file need to have FASTA files by RNA type with their respective names.")

    st.image('imgs/zipfile.png')
