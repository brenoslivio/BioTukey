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
                    seqs.append(current_seq.replace(' ', '').upper().replace('U', 'T'))
                    current_seq = ''
                names.append(line[1:].strip())
            else:
                current_seq += line.strip()
        seqs.append(current_seq)

        self.df = pd.DataFrame({'name': names, 'seq': seqs})
    
    def desc(self):
        return self.df['seq'].apply(lambda x: len(x)).describe()

    def gc_content(self):

        gc = 0
        total_len = 0

        for i in range(len(self.df['seq'])):
            gc += self.df['seq'][i].count('G') + self.df['seq'][i].count('C')
            total_len += len(self.df['seq'][i])

        pct_gc = (gc/total_len)*100

        return pct_gc, gc, total_len





