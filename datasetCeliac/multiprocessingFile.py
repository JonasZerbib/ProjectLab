import multiprocessing
import numpy as np
from multiprocessing import Process
from sklearn.model_selection import KFold
from textwrap import wrap
import pandas as pd
import pickle
from scipy.special import expit as sigmoid
import os.path
import sys


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


# ############### FIRST PART: IMPORT AND TREAT DATA #######################

# ############ TRANSLATE FROM AMINO ACID TO ASCHLEYS FACTOR

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


# parameter: receive on amino acid sequence
# return: list of snippet of length k_mears
# This function split the amino acid sequence into list of snippet
def snippet_generator(aa_seq, k_mears):
    snippet_seq = ["" for x in range(len(aa_seq) - (k_mears - 1))]
    if k_mears == 4:
        for x in range(0, len(aa_seq) - (k_mears - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3]

    if k_mears == 5:
        for x in range(0, len(aa_seq) - (k_mears - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3] + aa_seq[x + 4]

    if k_mears == 6:
        for x in range(0, len(aa_seq) - (k_mears - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3] + aa_seq[x + 4] + aa_seq[x + 5]

    if k_mears == 7:
        for x in range(0, len(aa_seq) - (k_mears - 1), 1):
            snippet_seq[x] = aa_seq[x] + aa_seq[x + 1] + aa_seq[x + 2] + aa_seq[x + 3] + aa_seq[x + 4] + aa_seq[x + 5] + \
                             aa_seq[x + 6]

    return list(set(snippet_seq))


# parameters: list of snippets for one antibody
# return: list of aschley factor , for each element of list a list of ashleys factor for for one snippet
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


# parameter : dataframe with all sequences dna of patient , with k mears
# return: array of patient , for each row vector of ashcley facors of unique snippet + frequence of this snippet + 1
def get_patient_factors(patient, k_mears):
    patient_factor_array = []
    for i in range(len(patient)):
        cdr3_dna = patient.iloc[i]['CDR3_IMGT']
        # check if there is a not codon stop before translate to amino acid
        check_codon_stop = wrap(cdr3_dna, 3)
        if 'TAA' in check_codon_stop:
            continue
        elif 'TAG' in check_codon_stop:
            continue
        elif 'TGA' in check_codon_stop:
            continue
        else:
            seq_amino = translate(cdr3_dna)
            patient_factor_array += (atchley_factors(seq_amino, k_mears))
    patient_factor_array = np.array(patient_factor_array, dtype=float)
    # receive uniques snippets
    patient_factor_unique, num_of_each_type = np.unique(patient_factor_array, axis=0, return_counts=True)
    # calculate frequence of each snippet ln(q)
    frequence_of_each_snippet = num_of_each_type / np.sum(num_of_each_type)
    patient_factor_unique = np.c_[patient_factor_unique, frequence_of_each_snippet]
    patient_factor_unique = np.c_[patient_factor_unique, np.ones((len(patient_factor_unique), 1), float)]
    return patient_factor_unique


# ################ LOAD DATA AND CLEAN HERE ##############################

# get table of patient
def get_table_data(i):
    path = "vdjbase_data/PROCESSED_V1/P1/P1_I" + str(i + 1) + "_S1/P1_I" + str(i + 1) + "_S1_genotyped_db-pass.tab"
    patient = pd.read_table(path, usecols=[37], dtype=str)
    print("patient " + str(i))
    return patient.applymap(str)


# get table of patient
# def get_table_data(i):
#     path = "../vdjbase_data/P1_I" + str(i + 1) + "_S1/P1_I" + str(i + 1) + "_S1_genotyped_db-pass.tab"
#     patient = pd.read_table(path, usecols=[37], dtype=str)
#     return patient.applymap(str)


# clean the table of patient
def clean_data(patient):
    patient = patient[~patient.CDR3_IMGT.str.contains('-')]
    patient = patient[~patient.CDR3_IMGT.str.contains('N')]
    return patient[~patient.CDR3_IMGT.str.contains('nan')]


# principal function that will get each patient and will return it in a numeric form.
def get_dataset(start, end, dataset_x, k_mears):
    for i in range(start, end):
        # get our patients from our dataset
        patient = get_table_data(i)
        # clean dataset
        patient = clean_data(patient)
        # get array of atchleys factors for each patient according to his snippets.
        patient_atchley_factor = get_patient_factors(patient, k_mears)
        dataset_x.append(patient_atchley_factor)
    return dataset_x


# return an array of 96 patients, and for each patient his array that each row correspond to
# atchleys factors of a snippet.
def get_dataset_celiac_x(k_mears):
    dataset_x = list()
    for num in range(5):
        if num == 0:
            dataset_x = get_dataset(0, 31, dataset_x, k_mears)
        if num == 1:
            dataset_x = get_dataset(32, 35, dataset_x, k_mears)
        if num == 2:
            dataset_x = get_dataset(36, 60, dataset_x, k_mears)
        if num == 3:
            dataset_x = get_dataset(61, 95, dataset_x, k_mears)
        if num == 4:
            dataset_x = get_dataset(96, 100, dataset_x, k_mears)
    dataset_x = np.array(dataset_x)
    return dataset_x


# recover y numeric label
def set_numeric_value_label(labels):
    labels.applymap(str)
    label = list()
    for i in range(len(labels)):
        if labels.iloc[i]['Health Status'] == 'Healthy':
            label.append(0)
        else:
            label.append(1)
    return np.array(label).T


# return dataset y for label
def get_dataset_celiac_y():
    # get the csv file of label
    csv_path = r'P1_CELIAC_METADATA.csv'
    labels = pd.read_csv(csv_path, encoding='utf-8')
    # transform this label into numerics value.
    labels = set_numeric_value_label(labels)
    return labels


# return input x and y label in a numeric form
def get_dataset_celiac(k_mears):
    x_data = get_dataset_celiac_x(k_mears)
    y_data = get_dataset_celiac_y()
    return x_data, y_data


# ############################# MODEL LOGISTIC RERGRESSION ##################

# before each epoch, we want to shuffle example in order to avoid overfitting
def shuffle(trains_x, trains_y):
    s = np.arange(trains_x.shape[0])
    np.random.shuffle(s)
    trains_x = trains_x[s]
    trains_y = trains_y[s]
    return [trains_x, trains_y]


# sigmoid
# def sigmoid(x):
#     return 1 / (1 + np.exp(-x))


# forward part of our algorithm:
# arguments: w, b0, b16 and patient x
# return: dictonnary with the high score(y_hat), snippet that give the best score, and our parameters
def forward(parameters, x):
    logit = np.matmul(x, np.transpose(parameters))
    score = sigmoid(logit)  # sigmoid
    # find the snippet that will give the best score.
    y_hat = np.max(score)
    # index of snippet that give the best score
    indice_max = np.where(score == y_hat)
    indice_max = indice_max[0]
    # factors of snippets that give the best score
    x_input_max = x[indice_max]
    propo_ret = {"y_hat": y_hat, "parameters": parameters, "x_input_max": x_input_max, "score": score}
    return propo_ret


#  loss function
def loss_function(y, y_hat):
    loss = ((-y) * np.log(y_hat)) - ((1 - y) * np.log(1 - y_hat))
    return loss


# calculate gradient for each parameter:
def gradients_for_each_parameters(param):
    y, y_hat, w, x = [param[key] for key in ('y', 'y_hat', 'parameters', 'x_input_max')]
    lnq = x[0][-2]
    x_factor = x[0][:-2]
    # derivative of w (weight), bias (b0) and bq (parameter for frequence of the snippet in this patient)
    dw = -(y - y_hat) * x_factor  # gradient of w
    db0 = -(y - y_hat)  # gradient of bias
    dbq = -(y - y_hat) * lnq  # gradient of lnq (q: frequence of snippet by patient)
    gradient = np.zeros((1, len(np.transpose(x))), float)
    gradient[0][:-2] = dw
    gradient[0][-1] = db0
    gradient[0][-2] = dbq
    return gradient


# GRADIENTS DESCENT OPTIMIZER
# update rules for logistic regression. in w there is a vector weight, bias, and bq
def update_parameters(gradient, lr, param):
    w = param['parameters']
    w = w - (lr * gradient)
    return w


# trainning our model
def train(parameters, epochs, lr, threshold, trains_x, trains_y, test_x, test_y, index, send_end):
    results = {}
    # loss_per_epoch = list()
    accuracy_per_Fold = list()
    acc = 0
    for i in range(epochs):
        correct = 0
        example = 0
        sumLoss = 0.0
        print("Epoch no. {0}".format(i + 1))
        trains_x, trains_y = shuffle(trains_x, trains_y)
        for x, y in zip(trains_x, trains_y):
            param = forward(parameters, x)
            param['y'] = y
            y_hat = param["y_hat"]
            # check if the score > thresold
            if y_hat > threshold:
                if 1 == y:
                    correct += 1
                example += 1
                loss = loss_function(y, y_hat)
                sumLoss += loss
                gradient = gradients_for_each_parameters(param)
                parameters = update_parameters(gradient, lr, param)
            # otherwise, check if the score < thresold
            else:
                if 0 == y:
                    correct += 1
                example += 1
                loss = loss_function(y, y_hat)
                sumLoss += loss
                gradient = gradients_for_each_parameters(param)
                parameters = update_parameters(gradient, lr, param)
        if (i % 100) == 0:
            print("{}/{} : {}".format(correct, example, (correct / example)))
            accuracy_per_Fold.append("{}/{} : {}".format(correct, example, round(correct / example, 2)))
        if i == 0:
            weightVector = parameters
        if acc <= correct / example:
            weightVector = parameters
            acc = correct / example

    # testing our model
    correct_test = 0
    example_test = 0
    param = forward(parameters, test_x)
    param['y'] = test_y
    y_hat = param["y_hat"]
    # check if the score > thresold
    if y_hat > threshold:
        if 1 == test_y:
            correct_test += 1
        example_test += 1
    else:
        if 0 == test_y:
            correct_test += 1
        example_test += 1
    results["acc_test"] = "{}/{} : {}".format(correct_test, example_test, round(correct_test / example_test, 2))
    results["index"] = index
    results["best_acc"] = acc
    results["best_weight"] = weightVector
    results["correct_test"] = correct_test

    send_end.send(results)


# ######## RECEIVE DATA, TRAIN AND TEST OUR LOGISTIC REGRESSION ###########

# in our parameters, we have the k_mear*5 first parameters for weight , and the 2 lasts parameters for the frequence
# of snippet and for the bias
def initialization_of_parameters(k_mears, num_seed):
    np.random.seed(num_seed)
    w = np.random.uniform(-0.5, 0.5, [1, k_mears * 5 + 2])
    w[0][-1] = 0
    return w


if __name__ == "__main__":
    k_mers = int(sys.argv[1])
    numOfFile = int(sys.argv[2])
    epochs = int(sys.argv[3])
    numOfFold = int(sys.argv[4])
    num_seed = int(sys.argv[5])
    part = int(sys.argv[6])

    x, y = get_dataset_celiac(k_mers)
    x_train, x_test, y_train, y_test, processes, pipe_list = [], [], [], [], [], []

    # Split data for leave one out
    if part == 1:
        for i in range(0, 24):
            x_1 = x
            y_1 = y
            x_test.append(x[i])
            y_test.append(y[i])
            x_train.append(np.delete(x_1, i))
            y_train.append(np.delete(y_1, i))

    elif part == 2:
        for i in range(24, 48):
            x_1 = x
            y_1 = y
            x_test.append(x_1[i])
            y_test.append(y_1[i])
            x_train.append(np.delete(x_1, i))
            y_train.append(np.delete(y_1, i))

    elif part == 3:
        for i in range(48, 72):
            x_1 = x
            y_1 = y
            x_test.append(x_1[i])
            y_test.append(y_1[i])
            x_train.append(np.delete(x_1, i))
            y_train.append(np.delete(y_1, i))

    else:
        for i in range(72, 96):
            x_1 = x
            y_1 = y
            x_test.append(x_1[i])
            y_test.append(y_1[i])
            x_train.append(np.delete(x_1, i))
            y_train.append(np.delete(y_1, i))

    # create processes
    for i in range(0, 24):
        w = initialization_of_parameters(k_mers, num_seed)
        recv_end, send_end = multiprocessing.Pipe(False)
        processes.append(Process(target=train, args=(
            w, epochs, 0.01, 0.5, x_train[i], y_train[i], x_test[i], y_test[i], i, send_end)))
        pipe_list.append(recv_end)
    # start processes in parallele
    for i in range(0, 24):
        processes[i].start()
    for i in range(0, 24):
        processes[i].join()
    result_list = [x.recv() for x in pipe_list]

    fileName = "Results" + str(k_mers) + "Part" + str(part) + "File" + str(numOfFile)
    savefile = open(fileName, "wb")
    pickle.dump(result_list, savefile)
    savefile.close()
