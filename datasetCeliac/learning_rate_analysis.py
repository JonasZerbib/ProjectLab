import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def learn_rate_analysis(curve):
    loss_per_epoch = curve['loss_per_epoch']
    accuracy_per_epoch = curve['accuracy_per_epoch']
    epoch_num = len(loss_per_epoch)
    k_mears = int((len(curve['best_parameters'][0]) - 2) / 5)

    title1 = "Epochs based learning accuracy analysis " + '\n' + '(k-mers = ' + str(
        k_mears) + ' ' + '| epochs number = ' + str(epoch_num) + ')'

    title2 = "Epochs based learning loss function analysis " + '\n' + '(k-mers = ' + str(
        k_mears) + ' ' + '| epochs number = ' + str(epoch_num) + ')'

    if epoch_num > 100:
        label = round(epoch_num / 100)
        loss_per_epoch = loss_per_epoch[::label]
        accuracy_per_epoch = accuracy_per_epoch[::label]
        new_len = list(range(epoch_num))
        xticks = new_len[::label]

        plt.figure()
        plt.plot(xticks, loss_per_epoch, label='loss', color='blue')
        plt.title(title2)
        plt.xlabel('Epoch number')
        plt.legend(loc='best')
        fig11 = plt.gcf()

        plt.figure()
        plt.plot(xticks, accuracy_per_epoch, label='accuracy (%)', color='orangered')
        plt.title(title1)
        plt.xlabel('Epoch number')
        plt.legend(loc='best')
        fig12 = plt.gcf()

        return fig11, fig12


    else:

        plt.figure()
        plt.plot(loss_per_epoch, label='loss', color='blue')
        plt.title(title2)
        plt.xlabel('Epoch number')
        plt.legend(loc='best')
        fig11 = plt.gcf()

        plt.figure()
        plt.plot(accuracy_per_epoch, label='accuracy (%)', color='orangered')
        plt.title(title1)
        plt.xlabel('Epoch number')
        plt.legend(loc='best')
        fig12 = plt.gcf()

        return fig11, fig12