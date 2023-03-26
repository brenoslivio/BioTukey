import streamlit as st
import base64

@st.cache_data
def get_base64(png_file):
    with open(png_file, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

def inject_css():
    with open("css/style.css", "r") as f:
        st.markdown(f'''<style>
                    [data-testid="stSidebar"] {{
                    background-image: url("data:image/png;base64,{get_base64('imgs/biotukey_logo.png')}");
                    padding-top: 150px;
                    background-repeat: no-repeat;
                    background-position: 50% 1%;
                    margin-top: 0%;
                    background-size: 80%;
                    }}
                    {f.read()}</style>''', unsafe_allow_html=True)