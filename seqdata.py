from Bio import SeqIO
import pandas as pd
from itertools import product
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

class Seq:
    def __init__(self, fasta, seq_class, seq_type):
        self.fasta = fasta
        self.seq_class = seq_class
        self.seq_type = seq_type

        names = []
        seqs = []

        if seq_type == "DNA/RNA":
            for record in SeqIO.parse(fasta, "fasta"):
                names.append(record.id)
                seqs.append(str(record.seq.back_transcribe()))
        elif seq_type == "Protein":
            for record in SeqIO.parse(fasta, "fasta"):
                names.append(record.id)
                seqs.append(str(record.seq))

        self.df = pd.DataFrame({'name': names, 'seq': seqs})
    
    def desc(self):
        return self.df['seq'].apply(lambda x: len(x)).describe()
    
    def avg_gc_content(self):
        pct_gc = (sum([gc_fraction(seq) for seq in self.df['seq']]) / len(self.df)) * 100

        return pct_gc
    
    def calc_N50(self):
        seq_lens = sorted([len(seq) for seq in self.df['seq']])
        half_length = sum(seq_lens) / 2

        cum_length = 0

        for seq_len in seq_lens:
            cum_length += seq_len
            if cum_length >= half_length:
                return seq_len


    def seq_total_len(self):
        return sum(self.df['seq'].str.len())
    
    def seq_len(self):
        seq_lens = [len(seq) for seq in self.df['seq']]

        return seq_lens

    def gc_content(self):
        gc = [gc_fraction(seq) * 100 for seq in self.df['seq']]

        return gc
    
    def kmer_count(self, k: int):
        
        bases = ['A', 'C', 'G', 'T']
        
        def dict_kmer(seq):
            counts = {''.join(comb): 0 for comb in product(bases, repeat= k)}
            L = len(seq) 
            for i in range(L - k + 1):
                counts[seq[i:i+k]] += 1
                
            return counts
        
        kmer_df = pd.DataFrame.from_dict(dict(self.df['seq'].apply(lambda x: dict_kmer(x))), orient='index')

        kmer_df = kmer_df.div(self.df['seq'].str.len() - k + 1, axis=0)

        avg_df = kmer_df.mean()
        avg_df.name = self.seq_class

        kmer_df.insert(0, 'nameseq', self.df['name'])
        kmer_df.insert(0, 'class', self.seq_class)

        return avg_df, kmer_df
        
        