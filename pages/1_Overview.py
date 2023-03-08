import streamlit as st
import statsmodels.stats.proportion as prop
from plotly.subplots import make_subplots
import plotly.express as px
import collections
import pandas as pd
import seqdata, css_injection
import os
import shutil
import subprocess
from stmol import showmol,render_pdb
import py3Dmol
import sys
import urllib.request
from pypdb import *
from stmol import showmol, render_pdb_resn
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from pyfamsa import Aligner, Sequence
import plotly.graph_objects as go
from Bio import SeqIO
import re

def show_protein():
    q = Query("MKELQTVLKNHFEIEFADKKLLETAFTHTSYANEHRLLKISHNERLEFLGDAVLQLLISEYLYKKYPKKPEGDLSKLRAMIVREESLAGFARDCQFDQFIKLGKGEEKSGGRNRDTILGDAFEAFLGALLLDKDVAKVKEFIYQVMIPKVEAGEFEMITDYKTHLQELLQVNGDVAIRYQVISETGPAHDKVFDVEVLVEGKSIGQGQGRSKKLAEQEAAKNAVEKGLDSCI", 
    query_type="sequence", 
    return_type="polymer_entity")

    id_pdb = q.search()["result_set"][0]["identifier"].split('_')[0]

    xyzview = py3Dmol.view(query=f'pdb:{id_pdb}')
    xyzview.setStyle({'cartoon':{'color':'spectrum'}})
    
    showmol(render_pdb_resn(viewer = xyzview, resn_lst = ['']), height = 500,width=800)# = ['ALA',]))

def process_files(uploaded_files):
    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey')

    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

    os.makedirs(dir_path)

    for file in uploaded_files:
        save_path = os.path.join(dir_path, file.name)
        with open(save_path, mode='wb') as f:
            f.write(file.getvalue())

    dir_path += '/'

    files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}

    seq_type = None

    alphabets = {'nt': re.compile('^[acgtu]*$', re.I), 
                'aa': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}

    for seq_class in files:

        pre_file = dir_path + 'pre_' + seq_class + '.fasta'

        with open(pre_file, mode='a') as f:
            for record in SeqIO.parse(files[seq_class], 'fasta'):

                if alphabets['nt'].search(str(record.seq)) is not None:
                    if seq_type == 'Protein':
                        st.error("Error: Dataset contains both DNA/RNA and Protein sequences.")
                        break
                    else:
                        if seq_type is None:
                            seq_type = 'DNA/RNA'
                        
                        f.write(f">{record.id}\n")
                        f.write(f"{record.seq}\n")
                elif alphabets['aa'].search(str(record.seq)) is not None:
                    if seq_type == 'DNA/RNA':
                        st.error("Error: Dataset contains both DNA/RNA and Protein sequences.")
                        break
                    else:
                        if seq_type is None:
                            seq_type = 'Protein'
                        
                        f.write(f">{record.id}\n")
                        f.write(f"{record.seq}\n")

        files[seq_class] = pre_file
            
        #subprocess.run(['python', 'MathFeature/preprocessing/preprocessing.py', '-i', files[seq_class], '-o', dir_path + 'pre_' + seq_class + '.fasta'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return files, seq_type

def load_study(study):
    files = {}
    seq_type = ''

    match study:
        case "ncRNAs":
            dir_path = "examples/example2/train/"
            files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}
            seq_type = "DNA/RNA"

    return files, seq_type

def seq_alignment(seqs):
    df_sequences = pd.DataFrame()

    for seq_class in seqs:
        df_seq = seqs[seq_class].df.copy()
        df_seq['list_name'] = seq_class + ' - ' + df_seq['name']
        df_sequences = pd.concat([df_sequences, df_seq]).reset_index(drop=True)

    seq_select = st.multiselect("Select sequences to view alignment:", df_sequences['list_name'], default =df_sequences['list_name'][0])

    if seq_select:
        with st.spinner('Loading...'):
            nseqs_selected = len(seq_select)

            df_select = df_sequences[df_sequences['list_name'].isin(seq_select)]

            sequences = [Sequence(name.encode(), 
                        df_select.loc[df_select['list_name'] == name]['seq'].item().encode()) 
                        for name in df_select['list_name']]

            aligner = Aligner(guide_tree="upgma")
            msa = aligner.align(sequences)
            
            msa_dict = {sequence.id.decode():sequence.sequence.decode() for sequence in msa}
            
            nucleotide_color = {"A": 0, "G": 0.25, "T": 0.5, "C": 0.75, "-": 1}

            msa_seqs = [[*msa_dict[id]] for id in msa_dict]

            seqs_num = [[nucleotide_color[N] for N in msa_seq] for msa_seq in msa_seqs]

            colorscale= [[0, 'rgb(217,95,2)'], [0.25, 'rgb(230,171,2)'], [0.5, 'rgb(27,158,119)'], [0.75, 'rgb(117,112,179)'], [1, 'white']]
            fig = make_subplots(rows=2, shared_xaxes=True, vertical_spacing=0.065, row_heights=[0.2, 0.8])

            np_seqs = np.array(msa_seqs)

            len_seqs = len(np_seqs[0])

            df = pd.DataFrame({'Position': list(range(len_seqs)),
                'A': [list(np_seqs[:,i]).count('A') for i in range(len_seqs)],
                'G': [list(np_seqs[:,i]).count('G') for i in range(len_seqs)],
                'T': [list(np_seqs[:,i]).count('T') for i in range(len_seqs)],
                'C': [list(np_seqs[:,i]).count('C') for i in range(len_seqs)]})

            colors = {'A': 'rgb(217,95,2)',
                    'G': 'rgb(230,171,2)',
                    'T': 'rgb(27,158,119)',
                    'C': 'rgb(117,112,179)'}

            data = []

            for i in range(df.shape[0]):
                ordered_columns = df.columns[1:][np.argsort(df.iloc[i, 1:].values)]

                text_col = []
                for col_index, column in enumerate(ordered_columns):
                    if col_index == 3:
                        for col in ordered_columns:
                            if df[col][i] == df[column][i]:
                                text_col.append(col)
                    
                    data.append(go.Bar(x=[df['Position'][i]],
                                    y=[df[column][i]],
                                    text = "<br>".join(text_col),
                                    textfont=dict(
                                        family="Trebuchet MS",
                                        size = 10,
                                    ),
                                    cliponaxis = False,
                                    marker=dict(color=colors[column]),
                                    name = column,
                                    textposition="outside",
                                    legendgroup=column,
                                    showlegend=i == 0,
                                    hovertemplate = "<b>Character: </b>" + column + " "
                                            "<br><b>Frequency:</b> %{y}" +
                                            "<br><b>Position:</b> %{x} <extra></extra>")) 

            fig_hm = px.imshow(seqs_num, aspect = 'auto')
            fig_hm.update_traces(text=msa_seqs,
                            hovertemplate = "<b>Sequence name:</b> %{y}" +
                                            "<br> <b>Character:</b> %{text}" +
                                            "<br> <b>Position:</b> %{x} <extra></extra>",
                            
                            textfont=dict(
                                    family="Trebuchet MS",
                                    size = 6,
                                ), texttemplate="<b>%{text}</b>")
            
            for bar in data:
                fig.add_trace(bar, row=1, col=1)

            fig.add_trace(fig_hm['data'][0], row = 2, col = 1)

            fig.update_coloraxes(colorscale = colorscale, cmin = 0, cmax = 1, showscale=False)
            fig.update_layout(height = 1000, barmode = 'stack', dragmode='pan',
                            legend_title_text='Characters',
                            yaxis2 = dict(tickmode = 'array',
                                tickvals = list(range(nseqs_selected)),
                                ticktext = list(msa_dict.keys())),
                            xaxis2=dict(
                            rangeslider=dict(visible=True), range = (-0.5, 50), type="linear",))

            st.plotly_chart(fig, use_container_width=True)

def kmer_general_stats(k, seqs, seq_type):
    with st.spinner('Loading...'):
        avgs_df = pd.DataFrame()
        kmers_df = pd.DataFrame()
        
        for seq_class in seqs:
            avg_df, kmer_df = seqs[seq_class].kmer_count(int(k))
            avgs_df = pd.concat([avgs_df, avg_df], axis = 1)
            kmers_df = pd.concat([kmers_df, kmer_df]).reset_index(drop=True)

        fig = px.bar(avgs_df, barmode='group', color_discrete_sequence = px.colors.qualitative.Dark2,
                        labels={
                            "index": "k-mer",
                            "value": "Average proportion"
                        })

        fig.update_layout(
            height=500,
            width=900,
            title_text= f"k-mer average proportion by {seq_type} class",
            legend_title_text= f"{seq_type} class",
            #legend = {"orientation":'h'}
        )

        st.plotly_chart(fig, use_container_width=True)

        return kmers_df
    
def kmer_individual_stats(kmers_df):
    st.markdown('######\n **k-mer individual proportion for sequences**')

    kmers_df = kmers_df.set_index("class")

    st.dataframe(kmers_df, use_container_width=True)

def nt_distribution(seqs, seq_type):
    st.markdown('### Nucleotide distribution insights')

    with st.spinner('Loading...'):
        df_plot = pd.DataFrame()

        for seq_class in seqs:
            new_df = pd.DataFrame(seqs[seq_class].df['seq'])
            new_df['class'] = seq_class
            
            df_plot = pd.concat([df_plot, new_df]).reset_index(drop=True)

        figures = {}

        for N in ['A', 'C', 'G', 'T']:
            df_plot[N] = df_plot['seq'].apply(lambda x : x.count(N) / len(x))
            figures[N] = px.violin(df_plot, x='class', y=N, color='class', color_discrete_sequence=px.colors.qualitative.Dark2)

        figures_traces = collections.defaultdict(list)

        for N in figures:
            for trace in range(len(figures[N]["data"])):
                figures_traces[N].append(figures[N]["data"][trace])

        fig = make_subplots(rows=2, cols=2,
            subplot_titles = ['Adenine', 'Cytosine', 'Guanine', 'Thymine'])

        for i, N in enumerate(['A', 'C', 'G', 'T']):
            for traces in figures_traces[N]:
                fig.add_trace(traces, row=(i//2) + 1, col=(i%2) + 1)

        fig.update_layout(height=575, width=575, showlegend=False, title_text=f"Violin plots for nucleotide proportion in sequences by {seq_type} class")

        st.plotly_chart(fig, use_container_width=True)

def summary_stats(seqs):
    df = pd.DataFrame()

    for seq_class in seqs:
        seq_desc = seqs[seq_class].desc()

        stats_df = pd.DataFrame({"class": seqs[seq_class].seq_class, 
                            "num_seqs": seq_desc['count'], 
                            "min_len (nt)": seq_desc['min'],  
                            "max_len (nt)": seq_desc['max'],
                            "avg_len (nt)": seq_desc['mean'],
                            "std_len (nt)": seq_desc['std'],
                            "sum_len (nt)": seqs[seq_class].seq_total_len(),
                            "Q1 (nt)": seq_desc['25%'],
                            "Q2 (nt)": seq_desc['50%'],
                            "Q3 (nt)": seq_desc['75%'],
                            "N50 (nt)": seqs[seq_class].calc_N50(),
                            "gc_content (%)": seqs[seq_class].avg_gc_content()}, 
                            index = [0])
        
        df = pd.concat([df, stats_df]).reset_index(drop=True)

    th_props = [
        ('text-align', 'center'),
        ('font-weight', 'bold')
        ]

    styles = [
        dict(selector="th", props=th_props)
    ]

    df = df.style.format({col: "{:.2f}" for col in df.columns if col != 'class'}).set_table_styles(styles)

    return df

def seq_stats(seqs):
    df = pd.DataFrame()

    for seq_class in seqs:
        
        stats_df = seqs[seq_class].gc_content()

        stats_df = pd.DataFrame({"class": seqs[seq_class].seq_class,
                            "nameseq": seqs[seq_class].df["name"],
                            "length (nt)": seqs[seq_class].seq_len(),     
                            "gc_content (%)": seqs[seq_class].gc_content()})

        df = pd.concat([df, stats_df]).reset_index(drop=True)

    th_props = [
        ('text-align', 'center'),
        ('font-weight', 'bold')
        ]

    styles = [
        dict(selector="th", props=th_props)
    ]

    cats = ['class', 'nameseq']

    df = df.set_index("class")
    df = df.style.format({col: "{:.2f}" for col in df.columns if col not in cats}).set_table_styles(styles)

    return df

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = "imgs/biotukey_icon.png", layout="wide")

    css_injection.inject_css()

    seqs = {}

    col1, col2 = st.columns([2, 6])

    with col1:
        #seq_type = st.selectbox("Select sequence type", ['', "DNA/RNA", "Protein"], 
        #                        help = "RNA sequences are back transcribed, being analyzed in the same way as DNA.")
        uploaded_files = st.file_uploader("Select your FASTA files by sequence class", 
                                        accept_multiple_files=True, type=["fasta", "fa", "faa"], 
                                        help="Each file must be named according to its class (e.g., sRNA.fasta). FASTA, FA, FAA files only.")

        study_example = st.selectbox("Or select study example", ['', "ncRNAs"])

        option = st.radio("Select option to load", ["Manual", "Example"], horizontal=True)

        match option:
            case "Manual":
                if not uploaded_files:
                    st.warning("Please select files.")
                else:
                    st.success("Files submitted with success.")
                    files, seq_type = process_files(uploaded_files)
            case "Example":
                if not study_example:
                    st.warning("Please select study example.")
                else:
                    st.success("Example submitted with success.")
                    files, seq_type = load_study(study_example)
    
    if (option == "Manual" and uploaded_files) or (option == "Example" and study_example):

        for seq_class in files:
            seq = seqdata.Seq(files[seq_class], seq_class, seq_type)
            seqs[seq_class] = seq

        df = summary_stats(seqs)

        with col2:
            if option == "Example":
                match study_example:
                    case "ncRNAs":
                        st.info("**Dataset from the following published paper:**\n \
                                Robson P Bonidia, Anderson P Avila Santos, Breno L S de Almeida, \
                                Peter F Stadler, Ulisses N da Rocha, Danilo S Sanches, \
                                André C P L F de Carvalho, BioAutoML: automated feature engineering \
                                and metalearning to predict noncoding RNAs in bacteria, \
                                Briefings in Bioinformatics, Volume 23, Issue 4, July 2022, \
                                bbac218, https://doi.org/10.1093/bib/bbac218")
                    
            st.markdown("Dataset provided has " + str(len(seqs)) + " " + seq_type + " class(es): " + ', '.join(seqs) + '. ' \
                        + 'Summary statistics for the sequences by class:')
            st.markdown("""  <div style="display: flex; justify-content: flex-end"><div class="tooltip"> 
                    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#66676e" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="icon">
                    <circle cx="12" cy="12" r="10"></circle><path d="M9.09 9a3 3 0 0 1 5.83 1c0 2-3 3-3 3"></path><line x1="12" y1="17" x2="12.01" y2="17"></line></svg>
                    <span class="tooltiptext">
                    <strong>num_seqs</strong>: Number of sequences;<br>
                    <strong>min_len</strong>: Minimum length of sequences;<br>
                    <strong>max_len</strong>: Maximum length of sequences;<br>
                    <strong>avg_len</strong>: Average length of sequences;<br>
                    <strong>std_len</strong>: Standard deviation for length of sequences;<br>
                    <strong>sum_len</strong>: Sum of length of all sequences;<br>
                    <strong>Q1</strong>: 25th percentile for length of sequences;<br>
                    <strong>Q2</strong>: 50th percentile for length of sequences;<br>
                    <strong>Q3</strong>: 75th percentile for length of sequences;<br>
                    <strong>N50</strong>: Length of the shortest read in the group of 
                    longest sequences that together represent (at least) 50% of the 
                    nucleotides in the set of sequences;<br>
                    <strong>gc_content</strong>: Average GC% content considering all sequences;</span>
                    </div></div> 
            """, unsafe_allow_html=True)

            st.table(df)

        tab1, tab2, tab3, tab4 = st.tabs(['Sequence statistics', 'Sequence alignment', 'k-mer distribution', 'Nucleotide distribution'])

        with tab1:
            df = seq_stats(seqs)

            st.dataframe(df, use_container_width=True)

        with tab2:
            seq_alignment(seqs)

        with tab3:
            st.markdown(f'### k-mer distribution for the {seq_type} class(es)')

            k = st.selectbox('Select size of k-mer:', ['1', '2', '3', '4', '5'])

            tab3_1, tab3_2 = st.tabs(['Average proportion', 'Individual proportion'])

            if k:
                with tab3_1:
                    kmers_df = kmer_general_stats(k, seqs, seq_type)
                with tab3_2:
                    kmer_individual_stats(kmers_df)

        with tab4:
            nt_distribution(seqs, seq_type)

        # st.markdown('---')


        # st.markdown('###### Hypothesis Testing')

        # st.markdown('Compare two sequence types to check for statistical significance in proportions related to GC% content.')

        # st.markdown('Null Hypothesis is:')

        # st.markdown('$$H_0: p_1 = p_2$$')

        # st.markdown('You can select the following Alternative Hypotheses:')

        # st.markdown('$$H_1: p_1 \\neq p_2; p_1 > p_2; p_1 < p_2$$')

        # with st.form("hypo_test"):
        #     st.markdown('Select the sequences types to be compared for GC% content:')

        #     type1 = st.selectbox('1st type ($p_1$):', [''] + list(seqs.keys()))

        #     type2 = st.selectbox('2nd type ($p_2$):', [''] + list(seqs.keys()))

        #     alternative = st.selectbox('Alternative Hypothesis:', ['two-sided', 'larger', 'smaller'])

        #     alpha = float(st.text_input('Significance level $\\alpha$:', 0.05))

        #     submitted = st.form_submit_button("Submit")
            
        #     st.markdown('---')

        #     if submitted and type1 and type2:    
        #         count1 = seqs[type1].nucleotide_count('G') + seqs[type1].nucleotide_count('C')
        #         nobs1 = seqs[type1].seq_total_len()
        #         count2 = seqs[type2].nucleotide_count('G') + seqs[type2].nucleotide_count('C')
        #         nobs2 = seqs[type2].seq_total_len()

        #         _, pvalue = prop.test_proportions_2indep(count1, nobs1, count2, nobs2, compare='diff', alternative=alternative)

        #         st.markdown('$p_1 = $ {:.4f}, $p_2 = $ {:.4f}'.format(seqs[type1].gc_content()/100, seqs[type2].gc_content()/100))

        #         st.markdown('***p*-value**: {:.4f}'.format(pvalue))

        #         if pvalue > alpha:
        #             st.markdown('Since *p*-value $> \\alpha$, we fail to reject the null hypothesis. At a {}% level of significance, there is not sufficient evidence to conclude the assumption of the alternative hypothesis.'.format(alpha*100))
        #         else:
        #             st.markdown('Since *p*-value $< \\alpha$, we reject the null hypothesis. At a {}% level of significance, the data support the assumption of the alternative hypothesis.'.format(alpha*100))

if __name__ == '__main__':
    runUI()