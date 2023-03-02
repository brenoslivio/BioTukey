import streamlit as st
import css_injection

def runUI():
    st.set_page_config(page_title = "BioTukey", page_icon = 'imgs/biotukey_icon.png', layout="wide")
    
    css_injection.inject_css()
  
if __name__ == '__main__':
    runUI()