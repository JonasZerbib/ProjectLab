import pandas as pd
import numpy as np
from textwrap import wrap


# ###############TRANSLATE  FROM DNA TO AMINO ACID############

# This function translate codon to amino acide
# input: codon
# output: amino acid
def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

# Atchley's factor table
table_atchleys = {

    "A": (-0.591, -1.302, -0.733, 1.570, -0.146),
    "C": (-1.343, 0.465, -0.862, -1.020, -0.255),
    "D": (1.050, 0.302, -3.656, -0.259, -3.242),
    "E": (1.357, -1.453, 1.477, 0.113, -0.837),
    "F": (-1.006, -0.590, 1.891, -0.397, 0.412),
    "G": (-0.384, 1.652, 1.330, 1.045, 2.064),
    "H": (0.336, -0.417, -1.673, -1.474, -0.078),
    "I": (-1.239, -0.547, 2.131, 0.393, 0.816),
    "K": (1.831, -0.561, 0.533, -0.277, 1.648),
    "L": (-1.019, -0.987, -1.505, 1.266, -0.912),
    "M": (-0.663, -1.524, 2.219, -1.005, 1.212),
    "N": (0.945, 0.828, 1.299, -0.169, 0.933),
    "P": (0.189, 2.081, -1.628, 0.421, -1.392),
    "Q": (0.931, -0.179, -3.005, -0.503, -1.853),
    "R": (1.538, -0.055, 1.502, 0.440, 2.897),
    "S": (-0.228, 1.399, -4.760, 0.670, -2.647),
    "T": (-0.032, 0.326, 2.213, 0.908, 1.313),
    "V": (-1.337, -0.279, -0.544, 1.242, -1.262),
    "W": (-0.595, 0.009, 0.672, -2.128, -0.184),
    "Y": (0.260, 0.830, 3.097, -0.838, 1.512)
}


# ################# CDR3 FEATURE ONLY FROM DATASET FOR EACH PATIENT ###########################


# This function split the amino acid sequence into list of snippet
# imput: receive on amino acid sequence and snippet length (kmers)
# output: list of snippet of length k_mears
def snippet_generator(aa_seq, kmers):
    snippet_seq = ["" for x in range(len(aa_seq) - (kmers - 1))]
    if kmers == 4:
        for x in range(0, len(aa_seq) - (kmers - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3]

    if kmers == 5:
        for x in range(0, len(aa_seq) - (kmers - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3] + aa_seq[x + 4]

    if kmers == 6:
        for x in range(0, len(aa_seq) - (kmers - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3] + aa_seq[x + 4] + aa_seq[x + 5]

    if kmers == 7:
        for x in range(0, len(aa_seq) - (kmers - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3] + aa_seq[x + 4] + aa_seq[x + 5] + \
                             aa_seq[x + 6]

    return list(set(snippet_seq))

# This function transforms the snippets of an antibody into vectors of atschley's factor
# input: list of snippet for one antibody
# output: list of aschley factor , for each element of list a list of ashleys factor for one snippet
def Atchleys_of_snippet(snippet_seq, k_mears):
    Atchley_table = {

        "A": (-0.591, -1.302, -0.733, 1.570, -0.146),
        "C": (-1.343, 0.465, -0.862, -1.020, -0.255),
        "D": (1.050, 0.302, -3.656, -0.259, -3.242),
        "E": (1.357, -1.453, 1.477, 0.113, -0.837),
        "F": (-1.006, -0.590, 1.891, -0.397, 0.412),
        "G": (-0.384, 1.652, 1.330, 1.045, 2.064),
        "H": (0.336, -0.417, -1.673, -1.474, -0.078),
        "I": (-1.239, -0.547, 2.131, 0.393, 0.816),
        "K": (1.831, -0.561, 0.533, -0.277, 1.648),
        "L": (-1.019, -0.987, -1.505, 1.266, -0.912),
        "M": (-0.663, -1.524, 2.219, -1.005, 1.212),
        "N": (0.945, 0.828, 1.299, -0.169, 0.933),
        "P": (0.189, 2.081, -1.628, 0.421, -1.392),
        "Q": (0.931, -0.179, -3.005, -0.503, -1.853),
        "R": (1.538, -0.055, 1.502, 0.440, 2.897),
        "S": (-0.228, 1.399, -4.760, 0.670, -2.647),
        "T": (-0.032, 0.326, 2.213, 0.908, 1.313),
        "V": (-1.337, -0.279, -0.544, 1.242, -1.262),
        "W": (-0.595, 0.009, 0.672, -2.128, -0.184),
        "Y": (0.260, 0.830, 3.097, -0.838, 1.512)
    }

    k_mears = k_mears
    atchleys_seq = []
    for x in range(len(snippet_seq)):
        letter_factors = []
        for i in range(k_mears):
            letter = snippet_seq[x][i]
            letter_factors += Atchley_table[letter]
        atchleys_seq.append(letter_factors)
    return atchleys_seq


def atchley_factors(aa_seq, k_mears):
    return Atchleys_of_snippet(snippet_generator(aa_seq, k_mears), k_mears)

# This function calls on the function to translate DNA sequences into amino acid sequences, then transform these sequences into atchley's factor
# parameter : dataframe with all sequences CDR3_IMGT dna of patient , with k mears
# return: array of patient , for each row vector of ashcley facors of unique snippet + frequence of this snippet + 1
def get_patient_factors(patient, clean_dataset, kmers):
    size_taa = 0
    size_tag = 0
    size_tga = 0
    patient_factor_array = []
    for i in range(len(patient)):
        cdr3_dna = patient.iloc[i]['CDR3_IMGT']
        # check if there is a not codon stop before translate to amino acid
        check_codon_stop = wrap(cdr3_dna, 3)
        if 'TAA' in check_codon_stop:
            size_taa += 1
            continue
        elif 'TAG' in check_codon_stop:
            size_tag += 1
            continue
        elif 'TGA' in check_codon_stop:
            size_tga += 1
            continue
        else:
            seq_amino = translate(cdr3_dna)
            patient_factor_array += (atchley_factors(seq_amino, kmers))
    # count the number of each type of stop codon we have removed
    clean_dataset['codon_stop_taa'] = size_taa
    clean_dataset['codon_stop_tag'] = size_tag
    clean_dataset['codon_stop_tga'] = size_tga
    patient_factor_array = np.array(patient_factor_array, dtype=float)
    # receive uniques snippets
    patient_factor_unique, num_of_each_type = np.unique(patient_factor_array, axis=0, return_counts=True)
    frequence_of_each_snippet = num_of_each_type / np.sum(num_of_each_type) # calculate frequence of each snippet
    patient_factor_unique = np.c_[patient_factor_unique, frequence_of_each_snippet]
    patient_factor_unique = np.c_[patient_factor_unique, np.ones((len(patient_factor_unique), 1), float)] # add one for bias
    return patient_factor_unique, clean_dataset


# ################ LOAD DATA AND CLEAN HERE ##############################

# This function load CDR3_IMGT column of the repertoire of a patient
def get_table_data(i):
    path = "vdjbase_data/PROCESSED_V1/P1/_/P1_I" + str(i + 1) + "_S1/P1_I" + str(i + 1) + "_S1_genotyped_db-pass.tab"
    patient = pd.read_table(path, usecols=[37], dtype=str)
    return patient.applymap(str)


# This function clean the table of patient.
# input: This function receive CDR3_IMGT columns of a repertoire's patient and the list cleaned_data in which we must add the unwanted DNA sequences
# output: the function return the patient cleaned of undesirable DNA sequence, and the numbers of each sequence that we erased.
def clean_data(patient, clean_dataset):
    size_patient = len(patient)
    clean_dataset["size_patient"] = size_patient
    patient = patient[~patient.CDR3_IMGT.str.contains('-')]
    clean_dataset['num_of_-'] = size_patient - len(patient)
    patient = patient[~patient.CDR3_IMGT.str.contains('N')]
    clean_dataset['num_of_N'] = size_patient - clean_dataset['num_of_-'] - len(patient)
    patient = patient[~patient.CDR3_IMGT.str.contains('nan')]
    clean_dataset['num_of_nan'] = size_patient - clean_dataset['num_of_-'] - clean_dataset['num_of_N'] - len(patient)
    return patient, clean_dataset


# This function calls all others function to apply feature engineering.
# This is the principal function that will get each patient and will return it in a numeric form.
def get_dataset(start, end, dataset_x, clean_dataset, k_mears):
    for i in range(start, end):
        clean_patient = {}
        # get our patients from our dataset
        patient = get_table_data(i)
        # clean dataset
        patient, clean_patient = clean_data(patient, clean_patient)
        # get array of atchleys factors for each patient according to his snippets.
        patient_atchley_factor, clean_patient = get_patient_factors(patient, clean_patient, k_mears)
        dataset_x.append(patient_atchley_factor)
        clean_dataset.append(clean_patient)
    return dataset_x, clean_dataset


# input: length of snippet (kmers)
# output an array of 96 patients, and for each patient his array that each row correspond to atchleys factors of a snippet.
def get_dataset_celiac_x(k_mears):
    dataset_x = list()
    clean_dataset = []
    for num in range(5):
        if num == 0:
            dataset_x, clean_dataset = get_dataset(0, 31, dataset_x, clean_dataset, k_mears)
        if num == 1:
            dataset_x, clean_dataset = get_dataset(32, 35, dataset_x, clean_dataset, k_mears)
        if num == 2:
            dataset_x, clean_dataset = get_dataset(36, 60, dataset_x, clean_dataset, k_mears)
        if num == 3:
            dataset_x, clean_dataset = get_dataset(61, 95, dataset_x, clean_dataset, k_mears)
        if num == 4:
            dataset_x, clean_dataset = get_dataset(96, 100, dataset_x, clean_dataset, k_mears)
    dataset_x = np.array(dataset_x)
    return dataset_x, clean_dataset

# input: labels contains dataframe with status of patients: healthy or sick
# output: y numeric label
def set_numeric_value_label(labels):
    labels.applymap(str)
    label = list()
    # one hot encoding methode:
    for i in range(len(labels)):
        # For healthy patient y=0
        if labels.iloc[i]['Health Status'] == 'Healthy':
            label.append(0)
        # For sick patient y=1
        else:
            label.append(1)
    return np.array(label).T

# This function load labels y
# output: dataset y for label
def get_dataset_celiac_y():
    # get the csv file of label
    csv_path = r'P1_CELIAC_METADATA.csv'
    labels = pd.read_csv(csv_path, encoding='utf-8')
    # transform this label into numerics value.
    labels = set_numeric_value_label(labels)
    return labels


# return input x and y label in a numeric form
# Return clean data with x_data in an another variable and labels y
def get_dataset_celiac(k_mears):
    x_data, clean_dataset = get_dataset_celiac_x(k_mears)
    y_data = get_dataset_celiac_y()
    return x_data, y_data
