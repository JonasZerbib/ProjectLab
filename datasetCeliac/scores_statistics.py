import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ####### FUNCTION THAT PLOT STATISTICS ABOUT SCORES ARRAYS FROM ALL HEALTHY AND SICK PATIENTS OBTAINED FROM TEST SET RESULTS
def fig_plot(curve):
    healthy_list_of_score_arrays = curve['healthy_list']

    healthy_scores = np.vstack(healthy_list_of_score_arrays)

    healthy_scores_mean = np.mean(healthy_scores)
    healthy_scores_mean = '%.2E' % healthy_scores_mean

    healthy_min_score = np.min(healthy_scores)
    healthy_min_score = '%.2E' % healthy_min_score

    healthy_max_score = np.max(healthy_scores)
    healthy_max_score = '%.2E' % healthy_max_score

    healthy_scores_std = np.std(healthy_scores)
    healthy_scores_std = '%.2E' % healthy_scores_std

    sick_list_of_score_arrays = curve['sick_list']

    sick_scores = np.vstack(sick_list_of_score_arrays)

    sick_scores_mean = np.mean(sick_scores)
    sick_scores_mean = '%.2E' % sick_scores_mean

    sick_min_score = np.min(sick_scores)
    sick_min_score = '%.2E' % sick_min_score

    sick_max_score = np.max(sick_scores)
    sick_max_score = '%.2E' % sick_max_score

    sick_scores_std = np.min(sick_scores)
    sick_scores_std = '%.2E' % sick_scores_std

    legend_h = 'Summary H: \n' + 'mean = ' + str(healthy_scores_mean) + '\n' + 'min = ' + str(
        healthy_min_score) + '\n' + 'max = ' + str(healthy_max_score) + '\n' + 'std = ' + str(healthy_scores_std)

    legend_s = 'Summary S: \n' + 'mean = ' + str(sick_scores_mean) + '\n' + 'min = ' + str(
        sick_min_score) + '\n' + 'max = ' + str(sick_max_score) + '\n' + 'std = ' + str(sick_scores_std)

    plt.hist(healthy_scores, bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], histtype='bar',
             edgecolor='black', color='cyan', alpha=0.5, label=legend_h)
    plt.hist(sick_scores, bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], histtype='bar', edgecolor='black',
             color='crimson', alpha=0.5, label=legend_s)

    plt.yscale('log')
    plt.title("Histogram of healthy and sick patients snippets scores")
    plt.legend(loc='best')
    fig6 = plt.gcf()
    return fig6
