<p align="center">
<img src="https://github.com/brenoslivio/BioTukey/blob/main/imgs/biotukey_logo.png?raw=true" alt="Logo BioTukey" width="200"/>
</p>

---

*BioTukey: a comprehensive toolkit for interactive data analysis, engineering, and classification of biological sequences*

---

## Summary 

- [Abstract](#1)

- [Project Description](#2)

- [Installation](#3)

- [License](#4)

---

## <a id="1" /> Abstract

This work proposes a toolkit named BioTukey, a graphical application that allows users to gather different kinds of insights about biological sequences such as DNA, RNA, and protein by providing a set of statistics, visualizations as multiple sequence alignment and sequence structure visualization, and machine learning approaches to proper show relevant characteristics and discriminate between types of these sequences. BioTukey provides different sequence descriptors, such as domain-specific biological and mathematical such as complex networks and entropy, to contribute to these insights. Machine learning models such as Random Forests and XGBoost can be selected along with choosing specific hyperparameters or using AutoML to recommend descriptors for classification. BioTukey can provide insightful perspectives to researchers working with multiple views on biological sequences, and even machine learning nonexperts could better apply ML pipelines considering the interactiveness of the application.

## <a id="2" /> Project Description

BioTukey was developed using Python programming, with the Streamlit framework as the foundation for its GUI. Streamlit is a web application framework that simplifies creating interactive data-driven applications. In this context, BioTukey leverages the capabilities of Streamlit to function as a web application, enabling users to interact with the application through its GUI. The application is designed as a responsive single-page application (SPA) to provide users with a user experience akin to a desktop application. It is organized into three main modules/pages: Overview, Feature Engineering, and Classification. Each module offers distinct functionalities to support analyzing and exploring biological sequences.

The Overview module provides comprehensive data visualization and summary statistics. It presents users with diverse statistics on the provided sequences, including sequence alignment visualization, $k$-mer distribution, and structure visualization for DNA/RNA and protein sequences. This module aims to provide users with a holistic understanding of the characteristics and properties of the input sequences.

The Feature Engineering module empowers users to extract different types of descriptors from the sequences, generating feature vectors that can be concatenated. These feature vectors serve as input for visualizations such as distribution plots, dimensionality reduction plots using techniques like PCA, t-SNE, or UMAP, correlation heatmaps, and feature importance bar charts. By allowing users to manipulate and explore the extracted features, this module facilitates a deeper understanding of the data and enables informed decision-making in subsequent analysis steps.

The Classification module allows users to perform classification tasks on the provided sequences. Users can select between two popular machine learning models, Random Forest, and XGBoost, as the algorithms for classification. They can further choose to apply 10-fold cross-validation, either with or without a separate test set. Additionally, users have the option to engage with the model manually by creating a custom pipeline. This includes selecting hyperparameters, filtering for the most important features, and applying dimensionality reduction techniques. The module also provides performance scores of the model, such as accuracy, precision, recall, and F1-score, along with a confusion matrix to evaluate the model's predictive capabilities.

In summary, BioTukey combines the capabilities of Python, Streamlit, and various data analysis techniques to create a user-friendly and comprehensive application for the analysis and classification of biological sequences. The application's modules, including Overview, Feature Engineering, and Classification, enable users to visualize, explore, and model biological data effectively, facilitating informed decision-making in bioinformatics research and applications.

## <a id="3" /> Installation

You need the packages in `environment.yml`. You can use Miniconda to easily configure the environment with the packages using:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

chmod +x Miniconda3-latest-Linux-x86_64.sh

./Miniconda3-latest-Linux-x86_64.sh

export PATH=~/miniconda3/bin:$PATH
```

After installing Miniconda, clone this repository using the following:

```bash
git clone https://github.com/brenoslivio/BioTukey.git BioTukey

cd BioTukey

git submodule init

git submodule update
```

This also guarantees that the MathFeature package is also appropriately configured. With Miniconda, create the environment with:

```bash
conda env create -f environment.yml -n biotukey
```

Activate the environment with:

```bash
conda activate biotukey
```

Finally, you can run the application with Streamlit using:

```bash
streamlit run main.py
```

You can access the application in a web browser in `http://localhost:8501`.

## <a id="4" /> License

Distributed under the MIT License. See `LICENSE` for more information.