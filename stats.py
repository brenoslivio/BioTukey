#from Bio import SeqIO
import pandas as pd

class SeqData:

    def __init__(self, fasta):
        self.fasta = fasta

        with open(fasta) as file:
            lines = file.readlines()

        names = []
        seqs = []
        current_seq = ''
        
        for line in lines:
            if line.startswith(">"):
                if current_seq != '':
                    seqs.append(current_seq.replace(' ', ''))
                    current_seq = ''
                names.append(line[1:].strip())
            else:
                current_seq += line.strip()
        seqs.append(current_seq)

        self.df = pd.DataFrame({'name': names, 'seq': seqs})
    
    def desc(self):
        return self.df['seq'].apply(lambda x: len(x)).describe()










