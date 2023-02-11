#from Bio import SeqIO
import pandas as pd

class Seq:
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

    def nucleotide_count(self, N):
        return sum(self.df['seq'].str.count(N))

    def seq_total_len(self):
        return sum(self.df['seq'].str.len())

    def gc_content(self):
        pct_gc = ((self.nucleotide_count('G') + self.nucleotide_count('C'))
                    /self.seq_total_len())*100

        return pct_gc




