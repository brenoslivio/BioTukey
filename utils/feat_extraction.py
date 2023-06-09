import subprocess, os, shutil
import pandas as pd

def feature_extraction(fasta_files, descriptors, seq_type):
    home_dir = os.path.expanduser('~')
    dir_path = os.path.join(home_dir, '.biotukey/feat_engineering/train')

    try:
        shutil.rmtree('.biotukey/feat_engineering')
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))
        print('Creating Directory...')

    os.makedirs(dir_path, exist_ok=True)

    features = pd.DataFrame()

    for seq_class in fasta_files:
        datasets = []

        if seq_type == "DNA/RNA":
            if "Nucleotide acid composition (NAC)" in descriptors:
                dataset = dir_path + '/NAC.csv'

                subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques.py',
                                '-i', fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-t', 'NAC', '-seq', '1'], 
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Dinucleotide composition (DNC)" in descriptors:
                dataset = dir_path + '/DNC.csv'

                subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-t', 'DNC', '-seq', '1'], 
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Trinucleotide composition (TNC)" in descriptors:
                dataset = dir_path + '/TNC.csv'

                subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-t', 'TNC', '-seq', '1'], 
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Xmer k-Spaced Ymer composition frequency (kGap)" in descriptors:
                dataset_di = dir_path + '/kGap_di.csv'
                dataset_tri = dir_path + '/kGap_tri.csv'

                subprocess.run(['python', 'MathFeature/methods/Kgap.py', '-i',
                                fasta_files[seq_class], '-o', dataset_di, '-l',
                                seq_class, '-k', '1', '-bef', '1',
                                '-aft', '2', '-seq', '1'],
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

                subprocess.run(['python', 'MathFeature/methods/Kgap.py', '-i',
                                fasta_files[seq_class], '-o', dataset_tri, '-l',
                                seq_class, '-k', '1', '-bef', '1',
                                '-aft', '3', '-seq', '1'],
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset_di)
                datasets.append(dataset_tri)

            if "Open Reading Frame (ORF)" in descriptors:
                dataset = dir_path + '/ORF.csv'

                subprocess.run(['python', 'MathFeature/methods/CodingClass.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class],
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Fickett Score" in descriptors:
                dataset = dir_path + '/Fickett.csv'

                subprocess.run(['python', 'MathFeature/methods/FickettScore.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-seq', '1'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Graphs" in descriptors:
                dataset = dir_path + '/ComplexNetworks.csv'

                subprocess.run(['python', 'MathFeature/methods/ComplexNetworksClass-v2.py', '-i', 
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class, 
                                '-k', '3'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Shannon entropy" in descriptors:
                dataset = dir_path + '/Shannon.csv'

                subprocess.run(['python', 'MathFeature/methods/EntropyClass.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-k', '5', '-e', 'Shannon'],
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Tsallis entropy" in descriptors:
                dataset = dir_path + '/Tsallis.csv'

                subprocess.run(['python', 'other-methods/TsallisEntropy.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-k', '5', '-q', '2.3'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)
        else:
            if "Amino acid composition (AAC)" in descriptors:
                dataset = dir_path + '/AAC.csv'
                subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques-Protein.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-t', 'AAC'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Dipeptide composition (DPC)" in descriptors:
                dataset = dir_path + '/DPC.csv'
                subprocess.run(['python', 'MathFeature/methods/ExtractionTechniques-Protein.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-t', 'DPC'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Graphs" in descriptors:
                dataset = dir_path + '/ComplexNetworks.csv'
                subprocess.run(['python', 'MathFeature/methods/ComplexNetworksClass-v2.py', '-i',
                                fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-k', '3'], stdout=subprocess.DEVNULL,
                                stderr=subprocess.STDOUT)
                datasets.append(dataset)
                                
            if "Shannon entropy" in descriptors:
                dataset = dir_path + '/Shannon.csv'
                subprocess.run(['python', 'MathFeature/methods/EntropyClass.py',
                                '-i', fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-k', '5', '-e', 'Shannon'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

            if "Tsallis entropy" in descriptors:
                dataset = dir_path + '/Tsallis.csv'
                subprocess.run(['python', 'other-methods/TsallisEntropy.py',
                                '-i', fasta_files[seq_class], '-o', dataset, '-l', seq_class,
                                '-k', '5', '-q', '2.3'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                datasets.append(dataset)

    if datasets:
        dataframes = pd.concat([pd.read_csv(f) for f in datasets], axis=1)
        dataframes = dataframes.loc[:, ~dataframes.columns.duplicated()]
        dataframes = dataframes[~dataframes.nameseq.str.contains("nameseq")]

    features = dataframes.reset_index(drop=True)

    return features