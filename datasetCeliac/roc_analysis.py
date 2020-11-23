import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def roc_curve(curve):  # receive a dictionnary that contains 4 keys: treshold, true_positif, false_positif,
    # true_negatif, false_negatif.

    treshold = curve['threshold']
    true_positif = curve['true_positive']
    false_positif = curve['false_positive']
    true_negatif = curve['true_negative']
    false_negatif = curve['false_negative']

    title = 'Receiver operating characteristic for tresh = ' + str(treshold)

    if true_positif + false_negatif == 0 or false_positif + true_negatif == 0 or true_positif + false_positif == 0:

        plt.figure()
        lw = 2
        plt.plot([0, 1], [0, 1], color='orangered', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.05])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.grid()
        plt.legend(loc='best', bbox_to_anchor=(1, 0.5),
                   title="INPUT DATA TO TEST NOT RELEVENT PLEASE CHANGE INPUT DATA")
        # plt.figure(figsize=(20, 20))

        fig4 = plt.gcf()
        return fig4




    else:

        true_pos_rate = true_positif / (true_positif + false_negatif)
        false_pos_rate = false_positif / (false_positif + true_negatif)
        precision = true_positif / (true_positif + false_positif)

        legend = 'Precision = ' + str(precision) + '\n' + 'T.P = ' + str(true_positif) + '\n' + 'T.N = ' + str(
            true_negatif) + '\n' + 'F.P = ' + str(false_positif) + '\n' + 'F.N = ' + str(false_negatif)

        print(legend)

        plt.figure()
        lw = 2
        plt.scatter([false_pos_rate], [true_pos_rate], color='m', s=50,
                    )
        plt.plot([0, 1], [0, 1], color='orangered', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.grid()
        plt.legend(loc='best', bbox_to_anchor=(1, 0.5), title=legend)
        # plt.figure(figsize=(20, 20))

        fig4 = plt.gcf()
        return fig4
