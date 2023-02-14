import statsmodels.stats.proportion as prop
from plotly.subplots import make_subplots
import plotly.express as px
import collections
import seqdata
import pandas as pd
import streamlit as st
import subprocess

def general_stats(seqs):
    st.markdown("You have " + str(len(seqs)) + " RNA types: " + ', '.join(seqs) + '. ' \
                        + 'We have the following summary statistics for the sequences:')

    df = pd.DataFrame()

    for type in seqs:
        seq_desc = seqs[type].desc()

        stats_df = pd.DataFrame({"type": seqs[type].type, 
                            "num_seqs": seq_desc['count'], 
                            "min_len (bp)": seq_desc['min'],
                            "avg_len (bp)": seq_desc['mean'],  
                            "max_len (bp)": seq_desc['max'],
                            "std_len (bp)": seq_desc['std'],
                            "Q1 (bp)": seq_desc['25%'],
                            "Q2 (bp)": seq_desc['50%'],
                            "Q3 (bp)": seq_desc['75%'],
                            "gc_content (%)": seqs[type].gc_content()}, 
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

    with st.expander("Description on the variables"):
        st.write("""
            - **num_seqs**: Number of sequences;\n
            - **min_len**: Minimum length of sequences; \n
            - **avg_len**: Average length of sequences;\n
            - **max_len**: Maximum length of sequences; \n
            - **std_len**: Standard deviation for length of sequences; \n
            - **Q1**: 25th percentile for length of sequences; \n
            - **Q2**: 50th percentile for length of sequences; \n
            - **Q3**: 75th percentile for length of sequences;  \n
            - **gc_content**: GC% content considering all sequences. \n
        """)

    st.markdown('---')
        
    col1, col2 = st.columns([3, 4])

    with col1:
        st.markdown('### Nucleotide distribution insights')

        df_plot = pd.DataFrame()

        for type in seqs:
            new_df = pd.DataFrame(seqs[type].df['seq'])
            new_df['type'] = type
            
            df_plot = pd.concat([df_plot, new_df]).reset_index(drop=True)

        figures = {}

        for N in ['A', 'C', 'G', 'T']:
            df_plot[N] = df_plot['seq'].apply(lambda x : x.count(N) / len(x))
            figures[N] = px.violin(df_plot, x='type', y=N, color='type', color_discrete_sequence=px.colors.qualitative.Dark2)

        figures_traces = collections.defaultdict(list)

        for N in figures:
            for trace in range(len(figures[N]["data"])):
                figures_traces[N].append(figures[N]["data"][trace])

        fig = make_subplots(rows=2, cols=2,
            subplot_titles = ['Adenine', 'Cytosine', 'Guanine', 'Thymine/Uracil'])

        for i, N in enumerate(['A', 'C', 'G', 'T']):
            for traces in figures_traces[N]:
                fig.append_trace(traces, row=(i//2) + 1, col=(i%2) + 1)

        fig.update_layout(height=575, width=575, showlegend=False, title_text="Violin plots for nucleotide proportion in sequences by RNA type")

        st.plotly_chart(fig)

        st.markdown('---')

        st.markdown('###### Hypothesis Testing')

        st.markdown('Compare two sequence types to check for statistical significance in proportions related to GC% content.')

        st.markdown('Null Hypothesis is:')

        st.markdown('$$H_0: p_1 = p_2$$')

        st.markdown('You can select the following Alternative Hypotheses:')

        st.markdown('$$H_1: p_1 \\neq p_2;\; p_1 > p_2;\; p_1 < p_2$$')

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

                st.markdown('$p_1 = $ {:.4f}, $\;p_2 = $ {:.4f}'.format(seqs[type1].gc_content()/100, seqs[type2].gc_content()/100))

                st.markdown('***p*-value**: {:.4f}'.format(pvalue))

                if pvalue > alpha:
                    st.markdown('Since *p*-value $> \\alpha$, we fail to reject the null hypothesis. At a {}% level of significance, there is not sufficient evidence to conclude the assumption of the alternative hypothesis.'.format(alpha*100))
                else:
                    st.markdown('Since *p*-value $< \\alpha$, we reject the null hypothesis. At a {}% level of significance, the data support the assumption of the alternative hypothesis.'.format(alpha*100))

    with col2:
        st.markdown('### k-mer distribution for the RNA types')

        k = st.selectbox('Select size of k-mer:', ['', '1', '2', '3', '4', '5'])

        if k:
            avgs_df = pd.DataFrame()
            kmers_df = pd.DataFrame()
            
            for type in seqs:
                avg_df, kmer_df = seqs[type].kmer_count(int(k))
                avgs_df = pd.concat([avgs_df, avg_df], axis = 1)
                kmers_df = pd.concat([kmers_df, kmer_df]).reset_index(drop=True)

            fig = px.bar(avgs_df, barmode='group', orientation='h', color_discrete_sequence = px.colors.qualitative.Dark2,
                            labels={
                                "index": "k-mer",
                                "value": "Average proportion"
                            })

            fig.update_layout(
                height=1300,
                width=800,
                title_text="k-mer average proportion by RNA type",
                legend_title_text="RNA type"
            )

            fig.update_yaxes(autorange="reversed")
            fig.update_xaxes(showgrid=True)


            st.plotly_chart(fig)

            st.markdown('---')

            st.markdown('###### k-mer individual proportion for sequences')

            kmers_df = kmers_df.set_index("type")

            st.dataframe(kmers_df)