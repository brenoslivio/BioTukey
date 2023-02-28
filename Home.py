import streamlit as st
import css_injection
import os
import shutil
import pandas as pd
import subprocess

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = ":microscope:", layout="wide")

    css_injection.inject_css()

    st.sidebar.markdown("test")

    st.title("BioTukey")
    st.markdown('''##### <span style="color:gray">Discover unique information about RNA sequences with an overall interactive environment and Machine Learning approaches.</span>
                ''', unsafe_allow_html=True)
                
if __name__ == '__main__':
    runUI()