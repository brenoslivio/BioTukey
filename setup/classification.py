import streamlit as st
import utils
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import mutual_info_classif, SelectKBest
from sklearn.preprocessing import StandardScaler
from functools import partial
from umap import UMAP
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import xgboost as xgb
import pandas as pd

def manual_model_selection(model_selection, test):
    le = preprocessing.LabelEncoder()

    X_train = st.session_state['data'][0]
    y_train = le.fit_transform(st.session_state['data'][2])

    pipeline = make_pipeline()

    if st.session_state['data'][5] == 'StandardScaler':
        pipeline.steps.append(('StandardScaler', StandardScaler()))
    else:
        pipeline.steps.append(('MinMaxScaler', MinMaxScaler()))

    manual_col1, manual_col2 = st.columns(2)

    with manual_col1:
        featimportance = st.slider("Number of features for the model based on feature importance:", 1, X_train.shape[1], X_train.shape[1])
        
        if featimportance < X_train.shape[1]:
            pipeline.steps.append(('feat', SelectKBest(score_func = partial(mutual_info_classif, random_state = 0), k=featimportance)))

    with manual_col2:
        dimension = st.selectbox("Dimensionality reduction:", ['None', 'PCA', 'UMAP'])

        if dimension != 'None':
            components = st.slider("Number of components:", 1, min(featimportance, X_train.shape[1]), min(featimportance, X_train.shape[1]))

            if dimension == 'PCA':
                pipeline.steps.append(('PCA', PCA(n_components=components)))
            elif dimension == 'UMAP':
                pipeline.steps.append(('UMAP', UMAP(n_components=components)))

    if model_selection == 'Random Forest':
        st.markdown("**Random Forest hyperparameters**")

        hyper_col1, hyper_col2 = st.columns(2)

        with hyper_col1:
            n_estimators = st.slider("Number of estimators:", 1, 500, 100)
            max_depth = st.slider("Maximum depth:", 1, 20, 5)
        
        with hyper_col2:
            min_samples_split = st.slider("Minimum samples to split:", 2, 20, 2)
            min_samples_leaf = st.slider("Minimum samples per leaf:", 1, 10, 1)
        
        model = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, random_state=0, n_jobs=-1)
        pipeline.steps.append(('RandomForest', model))

    elif model_selection == 'XGBoost':
        st.markdown("**XGBoost hyperparameters**")

        hyper_col1, hyper_col2 = st.columns(2)

        with hyper_col1:
            learning_rate = st.slider("Learning rate:", 0.01, 1.0, 0.1)
            max_depth = st.slider("Maximum depth:", 1, 20, 5)

        with hyper_col2:
            subsample = st.slider("Subsample ratio:", 0.1, 1.0, 1.0)
            colsample_bytree = st.slider("Column subsample ratio per tree:", 0.1, 1.0, 1.0)
        
        model = xgb.XGBClassifier(learning_rate=learning_rate, max_depth=max_depth, subsample=subsample, colsample_bytree=colsample_bytree, random_state=0, n_jobs=-1)
        pipeline.steps.append(('XGBoost', model))

    with st.spinner('Loading...'):

        tab1, tab2 = st.tabs(['Performance Scoring', 'Confusion Matrix'])
        
        with tab1:
            st.markdown("**Cross-validation**")

            metric_col1, metric_col2 = st.columns(2)

            with metric_col1:
                # Perform cross-validation prediction
                y_pred_cv = cross_val_predict(pipeline, X_train, y_train, cv=10, n_jobs=-1)
                report = classification_report(le.inverse_transform(y_train), le.inverse_transform(y_pred_cv), output_dict=True)
                
                st.markdown("Cross-validated metrics:")
                st.dataframe(pd.DataFrame(report).transpose(), use_container_width=True)

            with metric_col2:
                st.markdown("**Accuracy**")
                scores = cross_val_score(pipeline, X_train, y_train, cv=10, scoring='accuracy', n_jobs=-1)
                st.write("Mean cross-validation accuracy score:", scores.mean(), "+-", scores.std())
                st.markdown("**Precision**")
                scores = cross_val_score(pipeline, X_train, y_train, cv=10, scoring='precision_weighted', n_jobs=-1)
                st.write("Mean cross-validation weighted precision score:", scores.mean(), "+-", scores.std())
                st.markdown("**F1-score**")
                scores = cross_val_score(pipeline, X_train, y_train, cv=10, scoring='recall_weighted', n_jobs=-1)
                st.write("Mean cross-validation weighted recall score:", scores.mean(), "+-", scores.std())
                st.markdown("**F1-score**")
                scores = cross_val_score(pipeline, X_train, y_train, cv=10, scoring='f1_weighted', n_jobs=-1)
                st.write("Mean cross-validation weighted F1-score score:", scores.mean(), "+-", scores.std())

            if test:
                st.markdown("---")
                if (st.session_state['data'][0].shape[1] == st.session_state['test'][0].shape[1]):
                    
                    st.markdown("**Test set**")

                    X_test = st.session_state['test'][0]
                    y_test = le.transform(st.session_state['test'][1])

                    pipeline.fit(X_train.values, y_train)
                    y_pred = pipeline.predict(X_test)

                    report = classification_report(le.inverse_transform(y_test), le.inverse_transform(y_pred), output_dict=True)
                    
                    st.markdown("Test set metrics:")
                    st.dataframe(pd.DataFrame(report).transpose(), use_container_width=True)
                else:
                    st.warning("Submit your test sets.")

        with tab2:

            st.markdown("**Cross-validation**")

            y_pred_cv = cross_val_predict(pipeline, X_train, y_train, cv=10, n_jobs=-1)

            fig = go.Figure(data=go.Heatmap(
                z=confusion_matrix(le.inverse_transform(y_train), le.inverse_transform(y_pred_cv)),
                x=le.classes_,
                y=le.classes_,
                colorscale='Blues',
                colorbar=dict(title='Counts')
            ))

            fig.update_layout(
                title='Confusion Matrix',
                xaxis=dict(title='Predicted label'),
                yaxis=dict(title='True label'),
                width=500,
                height=500,
            )

            st.plotly_chart(fig, use_container_width=True)

            if test:
                st.markdown("---")
                if (st.session_state['data'][0].shape[1] == st.session_state['test'][0].shape[1]):

                    st.markdown("**Test set**")

                    X_test = st.session_state['test'][0]
                    y_test = le.transform(st.session_state['test'][1])

                    pipeline.fit(X_train.values, y_train)
                    y_pred = pipeline.predict(X_test)

                    fig = go.Figure(data=go.Heatmap(
                        z=confusion_matrix(le.inverse_transform(y_test), le.inverse_transform(y_pred)),
                        x=le.classes_,
                        y=le.classes_,
                        colorscale='Blues',
                        colorbar=dict(title='Counts')
                    ))

                    fig.update_layout(
                        title='Confusion Matrix',
                        xaxis=dict(title='Predicted label'),
                        yaxis=dict(title='True label'),
                        width=500,
                        height=500,
                    )

                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("Submit your test sets.")

def load(seq_type, option, study_example):
    if 'data' not in st.session_state:
        st.warning("Please select and submit descriptors to use for classification within the Feature Engineering module.")
    else:
        col1, col2 = st.columns(2)

        with col1:
            model_selection = st.selectbox("Select a model:", ['Random Forest', 'XGBoost'])

        with col2:
            evaluation_selection = st.selectbox("Select an evaluation method:", ['10-fold cross-validation', '10-fold cross-validation and test set'])

        tab1, tab2 = st.tabs(['Manual Model', 'AutoML'])

        if evaluation_selection == '10-fold cross-validation':
            with tab1:
                manual_model_selection(model_selection, False)

        elif evaluation_selection == '10-fold cross-validation and test set':
            with st.spinner('Loading...'):
                
                if option == "Example":
                    st.success("Test set loaded from study example successfully.")
                    files, _ = utils.processing.load_study(study_example, False)
                    features = utils.feature_extraction(files, st.session_state['data'][4], seq_type, False)

                    features.pop("nameseq")
                    labels = features.pop("label")

                    st.session_state['test'] = [features.values.astype(float), labels]
                else:
                    with st.form("test_set"):
                        uploaded_files = st.file_uploader("Submit test set files:", accept_multiple_files=True)

                        submitted = st.form_submit_button("Submit")
                    
                    if submitted and uploaded_files:
                        files, _ = utils.processing.process_files(uploaded_files, seq_type, False)
                        features = utils.feature_extraction(files, st.session_state['data'][4], seq_type, False)

                        features.pop("nameseq")
                        labels = features.pop("label")

                        st.session_state['test'] = [features.values.astype(float), labels]

            if 'test' in st.session_state:
                manual_model_selection(model_selection, True)

        with tab2:
            st.write("Under construction")  # Placeholder for AutoML functionality



