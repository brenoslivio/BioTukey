import streamlit as st
import pandas as pd
import os, shutil
import numpy as np
import subprocess
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.feature_selection import mutual_info_classif
from functools import partial
import plotly.figure_factory as ff
from umap import UMAP
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.express as px
import utils

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
            reducer = TSNE(n_components=3, perplexity=perplexity, learning_rate=learning_rate, n_iter=n_iter, random_state=0)
        elif reduction == "Uniform Manifold Approximation and Projection (UMAP)":
            n_neighbors = st.slider("Number of neighbors", min_value=2, max_value=100, value=15)
            min_dist = st.slider("Minimum distance", min_value=0.0, max_value=1.0, value=0.1)
            reducer = UMAP(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist, random_state=0)
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
                            color=utils.get_colors(len(np.unique(labels)))[i], size=3,
                        ),
                        hovertemplate=names[label],
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

def feature_importance(features, colnames, labels):

    with st.spinner('Loading...'):
        scores = partial(mutual_info_classif, random_state=0)(features, labels)

        # Sort features by importance
        sorted_indices = np.argsort(scores)
        sorted_colnames = colnames[sorted_indices]
        sorted_scores = scores[sorted_indices]

        fig = go.Figure(data=go.Bar(
            x=sorted_colnames,
            y=sorted_scores,
            marker=dict(color=sorted_scores, colorscale='purples'),
            hovertemplate='Feature: %{x}<br>Importance: %{y}<extra></extra>'
        ))
        fig.update_layout(
            title="Feature Importance using Mutual Information",
            xaxis_title="Features",
            yaxis_title="Importance",
        )

        st.plotly_chart(fig, use_container_width=True)

def feature_distribution(features, labels, nameseqs):
    col1, col2 = st.columns(2)

    # Select feature to plot
    with col1:
        selected_feature = st.selectbox("Select a feature", features.columns)

    # Get unique labels and assign colors
    unique_labels = labels.unique()
    color_map = utils.get_colors(len(unique_labels))[:len(unique_labels)]

    with col2:
        num_bins = st.slider("Number of bins", min_value=5, max_value=50, value=30)

    with st.spinner('Loading...'):
        fig_data = []
        fig_rug_text = []

        feature_data = features[selected_feature].values.astype(float)

        for label in unique_labels:
            group_indices = (labels == label).values
            group_data = feature_data[group_indices]
            fig_data.append(group_data)

            group_names = nameseqs[group_indices].tolist()
            fig_rug_text.append(group_names)

        bin_edges = np.histogram(fig_data[0], bins=num_bins)[1]

        fig = ff.create_distplot(
            fig_data,
            unique_labels,
            bin_size=bin_edges,
            colors=color_map,
            rug_text=fig_rug_text,
            histnorm="probability density"
        )

        fig.update_layout(
            title=f"Feature distribution for {selected_feature}",
            xaxis_title=selected_feature,
            yaxis_title="Density",
            height=800
        )

        st.plotly_chart(fig, use_container_width=True)

def load(files, seq_type):
    col1, col2 = st.columns(2)

    with col1:
        with st.form("feature_extraction"):
            if seq_type == "DNA/RNA":
                descriptors = st.multiselect("Select descriptors for feature engineering of DNA/RNA sequences provided:", 
                                            ["Nucleotide acid composition (NAC)", 
                                                "Dinucleotide composition (DNC)", 
                                                "Trinucleotide composition (TNC)", 
                                                "Xmer k-Spaced Ymer composition frequency (kGap)", 
                                                "Open Reading Frame (ORF)", "Fickett Score",
                                                "Graphs", "Shannon entropy", "Tsallis entropy"])
            else:
                descriptors = st.multiselect("Select descriptors for feature engineering of protein sequences provided:", 
                                            ["Amino acid composition (AAC)", 
                                                "Dipeptide composition (DPC)", 
                                                "Graphs", "Shannon entropy", "Tsallis entropy"])
            scaling = st.selectbox("Select scaling method", ["StandardScaler", "MinMaxScaler"])

            submitted = st.form_submit_button("Submit")
        
    with col2:
        if submitted:
            with st.spinner('Loading...'):
                features = utils.feature_extraction(files, descriptors, seq_type, False)

                match scaling:
                    case "StandardScaler":
                        scaler = StandardScaler()
                    case "MinMaxScaler":
                        scaler = MinMaxScaler()

                nameseqs = features.pop("nameseq")
                labels = features.pop("label")
                scaled_data = scaler.fit_transform(features)

                st.session_state['data'] = [features, scaled_data, labels, nameseqs, descriptors, scaling]
                
        if 'data' in st.session_state:
            data_shape = st.session_state['data'][0].shape
            st.markdown(f"{data_shape[0]} sequences with {data_shape[1]} features were generated from the descriptors. Preview of concatenated features from descriptors before scaling:")
            st.dataframe(st.session_state['data'][0])

    tab1, tab2, tab3, tab4 = st.tabs(["Feature Distribution", "Dimensionality Reduction", "Feature Correlation", "Feature Importance"])

    if 'data' in st.session_state:
        with tab1:
            feature_distribution(st.session_state['data'][0], st.session_state['data'][2], st.session_state['data'][3])

        with tab2:
            dimensionality_reduction(st.session_state['data'][1], st.session_state['data'][2], st.session_state['data'][3])

        with tab3:
            feature_correlation(st.session_state['data'][0])

        with tab4:
            feature_importance(st.session_state['data'][1], st.session_state['data'][0].columns, st.session_state['data'][2])