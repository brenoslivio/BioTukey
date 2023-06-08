import streamlit as st
import utils
import pandas as pd
import os, shutil
import numpy as np
import subprocess
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.manifold import TSNE
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif, RFE
from umap import UMAP
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.express as px
import utils

def feature_extraction(fasta_files, descriptors):
    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey/feat_engineering')

    try:
        shutil.rmtree(dir_path)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))
        print('Creating Directory...')

    os.makedirs(dir_path, exist_ok=True)

    features = pd.DataFrame()

    for seq_class in fasta_files:
        datasets = []

        if "Nucleotide acid composition (NAC)" in descriptors:
            dataset = dir_path + '/NAC.csv'

            subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques.py',
                            '-i', fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                            '-t', 'NAC', '-seq', '1'], 
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Dinucleotide composition (DNC)" in descriptors:
            dataset = dir_path + '/DNC.csv'

            subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques.py', '-i',
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                            '-t', 'DNC', '-seq', '1'], 
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Trinucleotide composition (TNC)" in descriptors:
            dataset = dir_path + '/TNC.csv'

            subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques.py', '-i',
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                            '-t', 'TNC', '-seq', '1'], 
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Xmer k-Spaced Ymer composition frequency (kGap)" in descriptors:
            dataset_di = dir_path + '/kGap_di.csv'
            dataset_tri = dir_path + '/kGap_tri.csv'

            subprocess.run(['python', 'MathFeature/methods/Kgap.py', '-i',
                            fasta_files[seq_class], '-o', dataset_di, '-l',
                            seq_class, '-k', '1', '-bef', '1',
                            '-aft', '2', '-seq', '1'],
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            subprocess.run(['python', 'MathFeature/methods/Kgap.py', '-i',
                            fasta_files[seq_class], '-o', dataset_tri, '-l',
                            seq_class, '-k', '1', '-bef', '1',
                            '-aft', '3', '-seq', '1'],
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset_di)
            datasets.append(dataset_tri)

        if "Open Reading Frame (ORF)" in descriptors:
            dataset = dir_path + '/ORF.csv'

            subprocess.run(['python', 'MathFeature/methods/CodingClass.py', '-i',
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class],
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Fickett Score" in descriptors:
            dataset = dir_path + '/Fickett.csv'

            subprocess.run(['python', 'MathFeature/methods/FickettScore.py', '-i',
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                            '-seq', '1'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Graphs" in descriptors:
            dataset = dir_path + '/ComplexNetworks.csv'

            subprocess.run(['python', 'MathFeature/methods/ComplexNetworksClass-v2.py', '-i', 
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class, 
                            '-k', '3'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Shannon entropy" in descriptors:
            dataset = dir_path + '/Shannon.csv'

            subprocess.run(['python', 'MathFeature/methods/EntropyClass.py', '-i',
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                            '-k', '5', '-e', 'Shannon'],
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)

        if "Tsallis entropy" in descriptors:
            dataset = dir_path + '/Tsallis.csv'

            subprocess.run(['python', 'other-methods/TsallisEntropy.py', '-i',
                            fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                            '-k', '5', '-q', '2.3'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            datasets.append(dataset)
        
    if datasets:
        dataframes = pd.concat([pd.read_csv(f) for f in datasets], axis=1)
        dataframes = dataframes.loc[:, ~dataframes.columns.duplicated()]
        dataframes = dataframes[~dataframes.nameseq.str.contains("nameseq")]

    features = dataframes.reset_index(drop=True)

    return features

def dimensionality_reduction(scaled_data, labels, nameseqs):
    dim_col1, dim_col2 = st.columns(2)

    with dim_col1:
        reduction = st.selectbox("Select dimensionality reduction technique", 
                                ["Principal Component Analysis (PCA)",
                                "t-Distributed Stochastic Neighbor Embedding (t-SNE)",
                                "Uniform Manifold Approximation and Projection (UMAP)"])
        
        if reduction == "t-Distributed Stochastic Neighbor Embedding (t-SNE)":
            perplexity = st.slider("Perplexity", min_value=5, max_value=50, value=30)
            learning_rate = st.slider("Learning rate", min_value=10, max_value=1000, value=200)
            n_iter = st.slider("Number of iterations", min_value=100, max_value=10000, value=1000)
            reducer = TSNE(n_components=3, perplexity=perplexity, learning_rate=learning_rate, n_iter=n_iter)
        elif reduction == "Uniform Manifold Approximation and Projection (UMAP)":
            n_neighbors = st.slider("Number of neighbors", min_value=2, max_value=100, value=15)
            min_dist = st.slider("Minimum distance", min_value=0.0, max_value=1.0, value=0.1)
            reducer = UMAP(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist)
        else:
            # No additional parameters for PCA
            reducer = PCA(n_components=3)

    names = {label: nameseqs[i] for i, label in enumerate(labels)}
            
    if reduction:
        with dim_col2:
            with st.spinner('Loading...'):

                reduced_data = reducer.fit_transform(scaled_data)

                fig = go.Figure()

                for i, label in enumerate(np.unique(labels)):
                    mask = labels == label
                    fig.add_trace(go.Scatter3d(
                        x=reduced_data[mask, 0],
                        y=reduced_data[mask, 1],
                        z=reduced_data[mask, 2],
                        mode='markers',
                        name=f'{label}',
                        marker=dict(
                            color=px.colors.qualitative.Dark2[i], size=3,
                        ),
                        hovertemplate=f"nameseq = {names[label]} <br> class = {label}",
                        textposition='top center',
                        hoverinfo='text'
                    ))

                fig.update_layout(
                    height=800,
                    scene=dict(
                        xaxis_title='Dimension 1',
                        yaxis_title='Dimension 2',
                        zaxis_title='Dimension 3'
                    ),
                    title=reduction
                )

                st.plotly_chart(fig, use_container_width=True)

def feature_correlation(features):
    correlation_method = st.selectbox('Select correlation method:', ['Pearson', 'Spearman'])

    if correlation_method == 'Pearson':
        correlation_matrix = features.corr(method='pearson')
    else:
        correlation_matrix = features.corr(method='spearman')
    
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Correlation between features sorted by the correlation coefficient**")
        correlation_df = pd.DataFrame(correlation_matrix.stack(), columns=['Correlation coefficient'])

        # Reset the index to get the feature pairs as separate columns
        correlation_df.reset_index(inplace=True)

        correlation_df.columns = ['Feature 1', 'Feature 2', 'Correlation coefficient']

        correlation_df = correlation_df.drop(correlation_df[correlation_df['Feature 1'] == correlation_df['Feature 2']].index)

        # Sort the DataFrame by correlation coefficient in descending order
        correlation_df = correlation_df.sort_values('Correlation coefficient', ascending=False).reset_index(drop=True)

        # Display the sorted dataframe
        st.dataframe(correlation_df, use_container_width=True)

    with col2:
        fig = px.imshow(correlation_matrix, x=correlation_matrix.columns, y=correlation_matrix.columns,
                        color_continuous_scale='RdBu', title='Correlation heatmap')
        
        fig.update_traces(hovertemplate='Feature 1 (x-axis): %{x}<br>Feature 2 (y-axis): %{y}<br>Correlation: %{z}<extra></extra>')

        fig.update_layout(height=800)
        st.plotly_chart(fig, use_container_width=True)

def feature_importance(features, labels):
    feature_selection_method = st.selectbox("Select feature selection method", ["Mutual Information", "Recursive Feature Elimination"])

    if feature_selection_method == "Mutual Information":
        scores = mutual_info_classif(features, labels)
    else:
        model = RandomForestClassifier() 
        rfe = RFE(model)
        rfe.fit(features, labels)
        scores = rfe.ranking_

    colnames = features.columns

    # Sort features by importance
    sorted_indices = np.argsort(scores)
    sorted_colnames = colnames[sorted_indices]
    sorted_scores = scores[sorted_indices]

    fig = go.Figure(data=go.Bar(
        x=sorted_colnames,
        y=sorted_scores,
        marker=dict(color=sorted_scores, colorscale='purples')
    ))
    fig.update_layout(
        title="Feature Importance",
        xaxis_title="Features",
        yaxis_title="Importance",
    )

    st.plotly_chart(fig, use_container_width=True)

def load(files, seq_type):
    col1, col2 = st.columns(2)

    with col1:
        with st.form("my_form"):
            descriptors = st.multiselect("Select descriptors for feature engineering of sequences provided:", 
                                        ["Nucleotide acid composition (NAC)", 
                                            "Dinucleotide composition (DNC)", 
                                            "Trinucleotide composition (TNC)", 
                                            "Xmer k-Spaced Ymer composition frequency (kGap)", 
                                            "Open Reading Frame (ORF)", "Fickett Score",
                                            "Graphs", "Shannon entropy", "Tsallis entropy"])
            scaling = st.selectbox("Select scaling method", ["StandardScaler", "MinMaxScaler"])

            submitted = st.form_submit_button("Submit")
        
    with col2:
        if submitted:
            with st.spinner('Loading...'):
                features = feature_extraction(files, descriptors)

                match scaling:
                    case "StandardScaler":
                        scaler = StandardScaler()
                    case "MinMaxScaler":
                        scaler = MinMaxScaler()

                nameseqs = features.pop("nameseq")
                labels = features.pop("label")
                colnames = features.columns
                scaled_data = scaler.fit_transform(features)

                st.session_state['data'] = [features, scaled_data, labels, nameseqs]
                
        if 'data' in st.session_state:
            data_shape = st.session_state['data'][0].shape
            st.markdown(f"{data_shape[0]} sequences with {data_shape[1]} features were generated from the descriptors. Preview of concatenated features from descriptors before scaling:")
            st.dataframe(st.session_state['data'][0])


    if 'data' in st.session_state:
        tab1, tab2, tab3 = st.tabs(["Dimensionality Reduction", "Feature Correlation", "Feature Importance"])

        with tab1:
            dimensionality_reduction(st.session_state['data'][1], st.session_state['data'][2], st.session_state['data'][3])

        with tab2:
            feature_correlation(st.session_state['data'][0])

        with tab3:
            feature_importance(st.session_state['data'][0], st.session_state['data'][2])