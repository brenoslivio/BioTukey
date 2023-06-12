
python utils/BioAutoML-feature.py --fasta_train examples/example1/train/sRNA.fasta examples/example1/train/tRNA.fasta examples/example1/train/rRNA.fasta --fasta_label_train sRNA tRNA rRNA --fasta_test examples/example1/test/sRNA.fasta examples/example1/test/tRNA.fasta examples/example1/test/rRNA.fasta --fasta_label_test sRNA tRNA rRNA --n_cpu -1 --output result

python utils/BioAutoML-feature-protein.py --fasta_train examples/example2/train/positive.fasta examples/example2/train/negative.fasta --fasta_label_train positive negative --fasta_test examples/example2/test/positive.fasta examples/example2/test/negative.fasta --fasta_label_test positive negative --n_cpu -1 --output result_protein
