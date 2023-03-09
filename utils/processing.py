import streamlit as st
import os
import re
from Bio import SeqIO
import shutil

def process_files(uploaded_files):
    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey')

    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

    os.makedirs(dir_path)

    for file in uploaded_files:
        save_path = os.path.join(dir_path, file.name)
        with open(save_path, mode='wb') as f:
            f.write(file.getvalue())

    dir_path += '/'

    files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}

    seq_type = None

    alphabets = {'nt': re.compile('^[acgtu]*$', re.I), 
                'aa': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}

    for seq_class in files:

        pre_file = dir_path + 'pre_' + seq_class + '.fasta'

        with open(pre_file, mode='a') as f:
            for record in SeqIO.parse(files[seq_class], 'fasta'):

                if alphabets['nt'].search(str(record.seq)) is not None:
                    if seq_type == 'Protein':
                        st.error("Error: Dataset contains both DNA/RNA and Protein sequences.")
                        return [], ""
                    else:
                        if seq_type is None:
                            seq_type = 'DNA/RNA'
                        
                        f.write(f">{record.id}\n")
                        f.write(f"{record.seq}\n")
                elif alphabets['aa'].search(str(record.seq)) is not None:
                    if seq_type == 'DNA/RNA':
                        st.error("Error: Dataset contains both DNA/RNA and Protein sequences.")
                        return [], ""
                    else:
                        if seq_type is None:
                            seq_type = 'Protein'
                        
                        f.write(f">{record.id}\n")
                        f.write(f"{record.seq}\n")

        files[seq_class] = pre_file
            
    return files, seq_type

def load_study(study):
    files = {}
    seq_type = ''

    match study:
        case "ncRNAs":
            dir_path = "examples/example2/train/"
            files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}
            seq_type = "DNA/RNA"

    return files, seq_type
