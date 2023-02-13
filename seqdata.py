#from Bio import SeqIO
import pandas as pd
from itertools import product

class Seq:
    def __init__(self, fasta, type):
        self.fasta = fasta
        self.type = type

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
    
    def kmer_count(self, k: int):
        
        bases = ['A', 'C', 'G', 'T']
        
        def dict_kmer(seq):
            counts = {''.join(comb): 0 for comb in product(bases, repeat= k)}

            for i in range(len(seq) - k + 1):
                counts[seq[i:i+k]] += 1
                
            return counts
        
        kmer_df = pd.DataFrame.from_dict(dict(self.df['seq'].apply(lambda x: dict_kmer(x))), orient='index')

        kmer_df = kmer_df.div(self.df['seq'].str.len() - k + 1, axis=0)

        avg_df = kmer_df.mean()
        avg_df.name = self.type

        kmer_df.insert(0, 'nameseq', self.df['name'])
        kmer_df.insert(0, 'type', self.type)

        return avg_df, kmer_df
        
        