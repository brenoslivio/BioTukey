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
from pymsaviz import MsaViz, get_msa_testdata

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
        with open(save_path, mode='wb') as w:
            w.write(file.getvalue())

    dir_path += '/'

    files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}

    for seq_class in files:
        subprocess.run(['python', 'MathFeature/preprocessing/preprocessing.py', '-i', files[seq_class], '-o', dir_path + 'pre_' + seq_class + '.fasta'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        files[seq_class] = dir_path + 'pre_' + seq_class + '.fasta'

    return files

def load_study(study):
    files = {}
    seq_type = ''

    match study:
        case "ncRNAs":
            dir_path = "examples/example2/train/"
            files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}
            seq_type = "DNA/RNA"

    return files, seq_type

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = "imgs/biotukey_icon.png", layout="wide")

    css_injection.inject_css()

    seqs = {}

    col1, col2 = st.columns([2, 6])

    with col1:
        seq_type = st.selectbox("Select sequence type", ['', "DNA/RNA", "Protein"], 
                                help = "RNA sequences are back transcribed, being analyzed in the same way as DNA.")
        uploaded_files = st.file_uploader("Select your FASTA files by sequence class", 
                                        accept_multiple_files=True, type=["fasta", "fa", "faa"], 
                                        help="Each file must be named according to its class (e.g., sRNA.fasta). FASTA, FA, FAA files only.")

        study_example = st.selectbox("Or select study example", ['', "ncRNAs"])

        option = st.radio("Select option to load", ["Manual", "Example"], horizontal=True)

        match option:
            case "Manual":
                if not seq_type or not uploaded_files:
                    st.warning("Please select type and files.")
                else:
                    st.success("Files submitted with success.")
                    files = process_files(uploaded_files)
            case "Example":
                if not study_example:
                    st.warning("Please select study example.")
                else:
                    st.success("Example submitted with success.")
                    study = study_example.split(':')[0]
                    files, seq_type = load_study(study)

    if (option == "Manual" and seq_type and uploaded_files) or (option == "Example" and study_example):

        for seq_class in files:
            seq = seqdata.Seq(files[seq_class], seq_class, seq_type)
            seqs[seq_class] = seq

        df = pd.DataFrame()

        for seq_class in seqs:
            seq_desc = seqs[seq_class].desc()

            stats_df = pd.DataFrame({"class": seqs[seq_class].seq_class, 
                                "num_seqs": seq_desc['count'], 
                                "min_len (bp)": seq_desc['min'],  
                                "max_len (bp)": seq_desc['max'],
                                "avg_len (bp)": seq_desc['mean'],
                                "std_len (bp)": seq_desc['std'],
                                "Q1 (bp)": seq_desc['25%'],
                                "Q2 (bp)": seq_desc['50%'],
                                "Q3 (bp)": seq_desc['75%'],
                                "gc_content (%)": seqs[seq_class].gc_content()}, 
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

        with col2:
            if option == "Example":
                match study:
                    case "ncRNAs":
                        st.info("**Dataset from the following paper:**\n \
                                Robson P Bonidia, Anderson P Avila Santos, Breno L S de Almeida, \
                                Peter F Stadler, Ulisses N da Rocha, Danilo S Sanches, \
                                André C P L F de Carvalho, BioAutoML: automated feature engineering \
                                and metalearning to predict noncoding RNAs in bacteria, \
                                Briefings in Bioinformatics, Volume 23, Issue 4, July 2022, \
                                bbac218, https://doi.org/10.1093/bib/bbac218")
                    
            st.markdown("You have " + str(len(seqs)) + " " + seq_type + " class(es): " + ', '.join(seqs) + '. ' \
                        + 'We have the following summary statistics for the sequences:')
            st.markdown("""  <div style="display: flex; justify-content: flex-end"><div class="tooltip"> 
                    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#66676e" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="icon">
                    <circle cx="12" cy="12" r="10"></circle><path d="M9.09 9a3 3 0 0 1 5.83 1c0 2-3 3-3 3"></path><line x1="12" y1="17" x2="12.01" y2="17"></line></svg>
                    <span class="tooltiptext">
                    <strong>num_seqs</strong>: Number of sequences;<br>
                    <strong>min_len</strong>: Minimum length of sequences;<br>
                    <strong>max_len</strong>: Maximum length of sequences;<br>
                    <strong>avg_len</strong>: Average length of sequences;<br>
                    <strong>std_len</strong>: Standard deviation for length of sequences;<br>
                    <strong>Q1</strong>: 25th percentile for length of sequences;<br>
                    <strong>Q2</strong>: 50th percentile for length of sequences;<br>
                    <strong>Q3</strong>: 75th percentile for length of sequences;<br>
                    <strong>gc_content</strong>: Average GC% content considering all sequences;</span>
                    </div></div> 
            """, unsafe_allow_html=True)

            st.table(df)

        tab1, tab2, tab3 = st.tabs(['Sequences information', 'k-mer distribution', 'Nucleotide distribution'])

        with tab1:
            #a = st.selectbox("Select sequence class", seqs['sRNA'].df)
        
            sequences = [Sequence(name.encode(), 
                        seqs['sRNA'].df.loc[seqs['sRNA'].df['name'] == name]['seq'].item().encode()) 
                        for name in seqs['sRNA'].df['name'].sample(10)]

            aligner = Aligner(guide_tree="upgma")
            msa = aligner.align(sequences)
 
            with open("align.fasta", "w") as f:
                for sequence in msa:
                    f.write(">" + sequence.id.decode() + "\n" + sequence.sequence.decode() + "\n")

            #msa_file = get_msa_testdata("align.fasta")
            # mv = MsaViz("align.fasta", wrap_length=60, show_count=True)
            # mv.savefig("api_example01.png")

        with tab2:
            st.markdown(f'### k-mer distribution for the {seq_type} class(es)')

            k = st.selectbox('Select size of k-mer:', ['1', '2', '3', '4', '5'])

            tab2_1, tab2_2 = st.tabs(['Average proportion', 'Individual proportion'])

            if k:
                with tab2_1:
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
                            legend = {"orientation":'h'}
                        )

                        st.plotly_chart(fig, use_container_width=True)
                with tab2_2:
                    st.markdown('######\n **k-mer individual proportion for sequences**')
                    kmers_df = kmers_df.set_index("class")

                    st.dataframe(kmers_df, use_container_width=True)

        with tab3:
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
                    subplot_titles = ['Adenine', 'Cytosine', 'Guanine', 'Thymine/Uracil'])

                for i, N in enumerate(['A', 'C', 'G', 'T']):
                    for traces in figures_traces[N]:
                        fig.add_trace(traces, row=(i//2) + 1, col=(i%2) + 1)

                fig.update_layout(height=575, width=575, showlegend=False, title_text=f"Violin plots for nucleotide proportion in sequences by {seq_type} class")

                st.plotly_chart(fig, use_container_width=True)

            st.markdown('---')


            st.markdown('###### Hypothesis Testing')

            st.markdown('Compare two sequence types to check for statistical significance in proportions related to GC% content.')

            st.markdown('Null Hypothesis is:')

            st.markdown('$$H_0: p_1 = p_2$$')

            st.markdown('You can select the following Alternative Hypotheses:')

            st.markdown('$$H_1: p_1 \\neq p_2; p_1 > p_2; p_1 < p_2$$')

            with st.form("hypo_test"):
                st.markdown('Select the sequences types to be compared for GC% content:')

                type1 = st.selectbox('1st type ($p_1$):', [''] + list(seqs.keys()))

                type2 = st.selectbox('2nd type ($p_2$):', [''] + list(seqs.keys()))

                alternative = st.selectbox('Alternative Hypothesis:', ['two-sided', 'larger', 'smaller'])

                alpha = float(st.text_input('Significance level $\\alpha$:', 0.05))

                submitted = st.form_submit_button("Submit")
                
                st.markdown('---')

                if submitted and type1 and type2:    
                    count1 = seqs[type1].nucleotide_count('G') + seqs[type1].nucleotide_count('C')
                    nobs1 = seqs[type1].seq_total_len()
                    count2 = seqs[type2].nucleotide_count('G') + seqs[type2].nucleotide_count('C')
                    nobs2 = seqs[type2].seq_total_len()

                    _, pvalue = prop.test_proportions_2indep(count1, nobs1, count2, nobs2, compare='diff', alternative=alternative)

                    st.markdown('$p_1 = $ {:.4f}, $p_2 = $ {:.4f}'.format(seqs[type1].gc_content()/100, seqs[type2].gc_content()/100))

                    st.markdown('***p*-value**: {:.4f}'.format(pvalue))

                    if pvalue > alpha:
                        st.markdown('Since *p*-value $> \\alpha$, we fail to reject the null hypothesis. At a {}% level of significance, there is not sufficient evidence to conclude the assumption of the alternative hypothesis.'.format(alpha*100))
                    else:
                        st.markdown('Since *p*-value $< \\alpha$, we reject the null hypothesis. At a {}% level of significance, the data support the assumption of the alternative hypothesis.'.format(alpha*100))

if __name__ == '__main__':
    runUI()