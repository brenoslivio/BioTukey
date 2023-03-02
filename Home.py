import streamlit as st
import css_injection
import pandas as pd
from PIL import Image

def runUI():
    image = Image.open('imgs/biotukey_icon.png')

    st.set_page_config(page_title = "BioTukey", page_icon = image, layout="wide")

    #st.image('imgs/biotukey_logo.png')

    css_injection.inject_css()

    st.title("BioTukey")
    st.markdown('''##### <span style="color:gray">Comprehensive toolkit for interactive data analysis, engineering, and classification of biological sequences.</span>
                ''', unsafe_allow_html=True)
                
if __name__ == '__main__':
    runUI()