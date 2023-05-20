import streamlit as st
import os
import re
from Bio import SeqIO


def process_files(uploaded_files):

    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey')

    for file in uploaded_files:
        save_path = os.path.join(dir_path, file.name)
        with open(save_path, mode='wb') as f:
            f.write(file.getvalue())

    dir_path += '/'

    files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}

    seq_type = None

    alphabets = {'nt': re.compile('^[acgtu]*$', re.I), 
                'aa': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I),
                'aa_exclusive': re.compile('[defhiklmpqrsvwy]', re.I)} # aa exclusive characters

    for seq_class in files:

        pre_file = dir_path + 'processed_' + seq_class + '.fasta'

        with open(pre_file, mode='a') as f:
            protein = False
            for record in SeqIO.parse(files[seq_class], 'fasta'):

                if alphabets['aa_exclusive'].search(str(record.seq)) is not None:
                    protein = True
                    break
            
            for record in SeqIO.parse(files[seq_class], 'fasta'):

                if protein:
                    if alphabets['aa'].search(str(record.seq)) is not None:
                        f.write(f">{record.id}\n")
                        f.write(f"{record.seq}\n")
                else: 
                    if alphabets['nt'].search(str(record.seq)) is not None:
                        f.write(f">{record.id}\n")
                        f.write(f"{record.seq}\n")

            if protein:
                seq_type = 'Protein'
            else:
                seq_type = 'DNA/RNA'

        files[seq_class] = pre_file
            
    return files, seq_type

def load_study(study):
    files = {}
    seq_type = ''

    match study:
        case "ncRNAs":
            dir_path = "examples/example1/train/"
            files = {os.path.splitext(f)[0] : dir_path + f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))}
            seq_type = "DNA/RNA"

    return files, seq_type
