from datasetCeliac.scores_statistics import fig_plot
from datasetCeliac.roc_analysis import roc_curve
from datasetCeliac.learning_rate_analysis import learn_rate_analysis
from datasetCeliac.biochemical_weights_analysis import biochemical_analysis
import matplotlib.pyplot as plt
import base64
import io
from urllib import parse
from datasetCeliac.cross_vaidation_analysis import cross_validation_accuracy


def presentation_of_results(results):
    plt.close('all')
    fig1 = biochemical_analysis(results)
    img1 = save_figures(fig1)
    plt.close('all')
    fig2, fig6 = learn_rate_analysis(results)
    img2 = save_figures(fig2)
    plt.close('all')
    img6 = save_figures(fig6)
    plt.close('all')
    fig3 = roc_curve(results)
    img3 = save_figures(fig3)
    plt.close('all')
    fig4 = fig_plot(results)
    img4 = save_figures(fig4)
    plt.close('all')
    return img1, img2, img3, img4, img6


def cross_val_results(last_acc, prediction):
    plt.close('all')
    fig1 = cross_validation_accuracy(last_acc, prediction)
    img1 = save_figures(fig1)
    plt.close('all')
    return img1


def save_figures(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='PNG', bbox_inches='tight')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    img_ready = parse.quote(string)
    buf.close()
    return img_ready
