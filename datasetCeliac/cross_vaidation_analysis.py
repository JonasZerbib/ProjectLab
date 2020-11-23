import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def cross_validation_accuracy(folds_accuracies, tests_results):
    height = folds_accuracies
    bars_int = folds_accuracies
    bars = [str(n) for n in bars_int]
    y_pos = np.arange(len(bars))

    # Create bars
    plt.bar(y_pos, height)

    # Create names on the x-axis

    legend = 'Model precision = ' + str(tests_results) + '/96'

    plt.xlabel('Fold_number')
    plt.ylabel('Accuracy of last epoch')
    plt.legend(loc='best', title=legend)
    plt.title('Chosen model analysis')
    # Show graphic
    fig7 = plt.gcf()
    return fig7

