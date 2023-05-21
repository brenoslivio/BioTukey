from Bio import SeqIO
import pandas as pd
from itertools import product
from Bio.SeqUtils import gc_fraction
import Bio.Seq
import numpy as np
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
    
    def orf_stats(self):
        table = 1

        num_orfs, min_len_orf, max_len_orf, avg_len_orf, std_len_orf  = [], [], [], [], []

        for seq in self.df['seq']:
            orfs_seq = 0
            orfs_len = []
            for strand, nuc in [(+1, Bio.Seq.Seq(seq)), (-1, Bio.Seq.Seq(seq).reverse_complement())]:
                for frame in range(3):
                    length = 3 * ((len(seq)-frame) // 3) #Multiple of three
                    for pro in nuc[frame:frame+length].translate(table).split("*"):
                        if len(pro) >= 30:
                            orfs_seq += 1
                            orfs_len.append(len(pro) * 3)

            num_orfs.append(orfs_seq)

            if len(orfs_len) > 0:
                min_len_orf.append(min(orfs_len))
                max_len_orf.append(max(orfs_len))
                avg_len_orf.append(np.mean(orfs_len))
                std_len_orf.append(np.std(orfs_len))
            else:
                min_len_orf.append(0)
                max_len_orf.append(0)
                avg_len_orf.append(0)
                std_len_orf.append(0)
                
        return num_orfs, min_len_orf, max_len_orf, avg_len_orf, std_len_orf
    
    def aa_stats(self):
        aliphatic = [((seq.count('G') + seq.count('A') + seq.count('V') + 
                      seq.count('L') + seq.count('I') + seq.count('P'))/len(seq))*100 for seq in self.df['seq']]
        
        aromatic = [((seq.count('F') + seq.count('W') + seq.count('Y'))/len(seq))*100 for seq in self.df['seq']]
        
        acidic = [((seq.count('D') + seq.count('E'))/len(seq))*100 for seq in self.df['seq']]

        basic = [((seq.count('K') + seq.count('R') + seq.count('H'))/len(seq))*100 for seq in self.df['seq']]

        hidroxylic = [((seq.count('S') + seq.count('T'))/len(seq))*100 for seq in self.df['seq']]

        sulphur = [((seq.count('C') + seq.count('M'))/len(seq))*100 for seq in self.df['seq']]

        amidic = [((seq.count('N') + seq.count('Q'))/len(seq))*100 for seq in self.df['seq']]
        
        return aliphatic, aromatic, acidic, basic, hidroxylic, sulphur, amidic

    def kmer_count(self, k):
        
        if self.seq_type == "DNA/RNA":
            chars = ['A', 'C', 'G', 'T']
        else:
            chars = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

        def dict_kmer(seq):
            counts = {''.join(comb): 0 for comb in product(chars, repeat= k)}
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
        
        