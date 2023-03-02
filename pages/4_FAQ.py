import streamlit as st
import css_injection

css_injection.inject_css()

st.markdown("**1. How my ZIP file has to look like for me to upload it?**")

st.markdown("Your ZIP file need to have FASTA files by RNA type with their respective names.")

st.image('imgs/zipfile.png')