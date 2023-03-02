import streamlit as st
import css_injection
import pandas as pd
from PIL import Image

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = 'imgs/biotukey_icon.png', layout="wide")

    css_injection.inject_css()

    st.title("BioTukey")
    st.markdown('''##### <span style="color:gray">Comprehensive toolkit for interactive data analysis, engineering, and classification of biological sequences.</span>
                ''', unsafe_allow_html=True)
                
if __name__ == '__main__':
    runUI()