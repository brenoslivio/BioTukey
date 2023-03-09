import streamlit as st
import utils

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = 'imgs/biotukey_icon.png', layout="wide")
    
    utils.inject_css()
  
if __name__ == '__main__':
    runUI()