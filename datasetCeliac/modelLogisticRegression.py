import numpy as np
from scipy.special import expit as sigmoid
from sklearn.model_selection import KFold


# before each epoch, we want to shuffle example in order to avoid overfitting
def shuffle(trains_x, trains_y):
    s = np.arange(trains_x.shape[0])
    np.random.shuffle(s)
    trains_x = trains_x[s]
    trains_y = trains_y[s]
    return [trains_x, trains_y]


# forward part of our algorithm:
# input: w, b0, b16 and patient x
# output: dictonnary with the high score(y_hat), snippet that give the best score, and our parameters
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
    gradient[0][:-2] = dw  # gradient of weight vector
    gradient[0][-1] = db0  # gradient of bias
    gradient[0][-2] = dbq  # gradient of bfreq
    return gradient


# GRADIENTS DESCENT OPTIMIZER
# update rules for logistic regression. in w there is a vector weight, bias, and bq
def update_parameters(gradient, lr, param):
    w = param['parameters']
    w = w - (lr * gradient)
    return w


# trainning our model
def train(parameters, epochs, lr, threshold, trains_x, trains_y):
    loss_per_epoch = list()
    accuracy_per_epoch = list()
    last_accuracy = 0.0
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

            else:
                if 0 == y:
                    correct += 1
                example += 1
                loss = loss_function(y, y_hat)
                sumLoss += loss
                gradient = gradients_for_each_parameters(param)
                parameters = update_parameters(gradient, lr, param)
        if i == (epochs-1):
            last_accuracy = (correct/example)
        loss_per_epoch.append(sumLoss / example)
        accuracy_per_epoch.append((correct / example) * 100)
    results = {"loss_per_epoch": loss_per_epoch, "accuracy_per_epoch": accuracy_per_epoch,
               "best_parameters": parameters, "last_acc": last_accuracy}
    return results


# testing our model
def test(validation_x, validation_y, threshold, data):
    true_positive = 0
    true_negative = 0
    false_positive = 0
    false_negative = 0
    sick = list()
    not_sick = list()
    parameters = data['best_parameters']
    for x, y in zip(validation_x, validation_y):
        result = forward(parameters, x)
        if y == 1:
            sick.append(result['score'])
        else:
            not_sick.append(result['score'])
        y_predict = result["y_hat"]
        if y_predict > threshold:
            if y == 1:
                true_positive += 1
            else:
                false_positive += 1
        else:
            if y == 0:
                true_negative += 1
            else:
                false_negative += 1
    results = {"threshold": threshold, "true_positive": true_positive, "false_positive": false_positive,
               "true_negative": true_negative,
               "false_negative": false_negative, 'sick_list': sick, 'healthy_list': not_sick}
    results.update(data)
    return results


# ######## RECEIVE DATA, TRAIN AND TEST OUR LOGISTIC REGRESSION ###########

def split_data(len_trainning_set, dataset_x, dataset_y):
    len_split = round((len_trainning_set / 100) * len(dataset_x))
    dataset_x, dataset_y = shuffle(dataset_x, dataset_y)
    train_x = dataset_x[:len_split]
    test_x = dataset_x[len_split:]
    train_y = dataset_y[:len_split]
    test_y = dataset_y[len_split:]
    return train_x, train_y, test_x, test_y


# in our parameters, we have the k_mear*5 first parameters for weight , and the 2 lasts parameters for the frequence
# of snippet and for the bias
def initialization_of_parameters(k_mears, sid):
    np.random.seed(sid)
    w = np.random.uniform(-0.5, 0.5, [1, k_mears * 5 + 2])
    w[0][-1] = 0
    return w


# This function test lasts parameters for the last patient in the method of leave on out cross validation
def test_crossVal(x_test, y_test, threshold, train_result):
    parameters = train_result['best_parameters']
    good_pred = 0
    for x, y in zip(x_test, y_test):
        result = forward(parameters, x)
        y_predict = result["y_hat"]
        if y_predict > threshold:
            if y == 1:
                good_pred = 1
        else:
            if y == 0:
                good_pred = 1
    return good_pred


# These function allows check model. The leave one out cross validation is used.
def crossVal(parameters, epochs, learning_rate, threshold, x_dataset, y_dataset):
    x_train, x_test, y_train, y_test, last_acc = [], [], [], [], []
    good_prediction = 0
    # Function from sklearn that allow to separate our dataset in 96 folders. For each index among 96, we keep only
    # one patient for test, and the 95 others for trainning
    kf = KFold(n_splits=96)
    for train_index, test_index in kf.split(x_dataset):
        x_train.append(x_dataset[train_index])
        x_test.append(x_dataset[test_index])
        y_train.append(y_dataset[train_index])
        y_test.append(y_dataset[test_index])
    # train the model and test it.
    for i in range(0, 96):
        train_result = train(parameters, epochs, learning_rate, threshold, x_train[i], y_train[i])
        last_acc.append(train_result['last_acc'])
        good_prediction = good_prediction + test_crossVal(x_test[i], y_test[i], threshold, train_result)

    return good_prediction, last_acc
