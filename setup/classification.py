import streamlit as st
import utils
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import mutual_info_classif, SelectKBest
from sklearn.preprocessing import StandardScaler
from functools import partial
from umap import UMAP
from sklearn.decomposition import PCA
import xgboost as xgb
import pandas as pd

def manual_model_selection(model_selection):
    le = preprocessing.LabelEncoder()

    X_train = st.session_state['data'][1]
    y_train = le.fit_transform(st.session_state['data'][2])

    pipeline = make_pipeline()

    if st.session_state['data'][5] == 'StandardScaler':
        pipeline.steps.append(('StandardScaler', StandardScaler()))
    else:
        pipeline.steps.append(('MinMaxScaler', MinMaxScaler()))

    manual_col1, manual_col2 = st.columns(2)

    with manual_col1:
        hyperparameter = st.selectbox("Model hyperparameters:", ['Manual', 'Automatic (Bayesian Optimization)'])

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
        n_estimators = st.slider("Number of estimators:", 1, 100, 10)
        max_depth = st.slider("Maximum depth:", 1, 20, 5)
        model = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, random_state=0, n_jobs=-1)
        pipeline.steps.append(('RandomForest', model))

    elif model_selection == 'XGBoost':
        st.markdown("**XGBoost hyperparameters**")
        learning_rate = st.slider("Learning rate:", 0.01, 1.0, 0.1)
        max_depth = st.slider("Maximum depth:", 1, 20, 5)
        model = xgb.XGBClassifier(learning_rate=learning_rate, max_depth=max_depth, random_state=0, n_jobs=-1)
        pipeline.steps.append(('XGBoost', model))

    with st.spinner('Loading...'):

        tab1, tab2 = st.tabs(['Performance Metrics', 'Confusion Matrix'])
        
        with tab1:
            st.markdown("**Cross-validation**")

            metric_col1, metric_col2 = st.columns(2)

            with metric_col1:
                # Perform cross-validation prediction
                y_pred_cv = cross_val_predict(pipeline, X_train, y_train, cv=10, n_jobs=-1)
                report = classification_report(le.inverse_transform(y_train), le.inverse_transform(y_pred_cv), output_dict=True)
                
                st.markdown("Cross-validated estimates for each input data point:")
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

def load(seq_type):
    if 'data' not in st.session_state:
        st.warning("Please select and submit descriptors to use for classification within the Feature Engineering module.")
    else:
        col1, col2 = st.columns(2)

        with col1:
            model_selection = st.selectbox("Select a model:", ['Random Forest', 'XGBoost'])

        with col2:
            evaluation_selection = st.selectbox("Select an evaluation method:", ['10-fold cross-validation', '10-fold cross-validation and test set'])

    if evaluation_selection == '10-fold cross-validation and test set':
        uploaded_files = st.file_uploader("Upload test set files:", accept_multiple_files=True)
        if uploaded_files:
            for file in uploaded_files:
                # TODO: Perform actions with the uploaded test set files
                pass

    tab1, tab2 = st.tabs(['Manual Model', 'AutoML'])

    with tab1:
        manual_model_selection(model_selection)
    with tab2:
        st.write("Under construction")  # Placeholder for AutoML functionality


