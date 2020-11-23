# ########### BIOCHEMICAL ANALYSIS OF THE DIAGNOSTIC WEIGHTS VECTOR #########
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ############ FUNCTION THAT RECEIVE THE BEST SNIPPET A.FACTORS AND RETURNS ITS A.A SEQUENCE AS A STRING ##########
def biochemical_analysis(curve):  # le vecteur poids recu en parametre n'a pas le bon type
    w_vec = curve['best_parameters']
    w_vec = w_vec[0]
    w_vec = np.array(w_vec)
    w = w_vec[:-2]
    b = w_vec[-2:]
    k_mears = len(w) / 5

    factors_1 = w[::5]
    factors_2 = w[1:][::5]
    factors_3 = w[2:][::5]
    factors_4 = w[3:][::5]
    factors_5 = w[4:][::5]

    b0 = b[1]
    b1 = b[0]

    if k_mears == 4:
        df = pd.DataFrame([factors_1, factors_2, factors_3, factors_4, factors_5, b],
                          index=['Polarity', 'Secondary\nstructure', 'Molecular\nvolume', 'Codon\ndiversity',
                                 'Electrostatic\ncharge', 'Bq | B0'],
                          columns=pd.Index(['1st A.A', '2nd A.A', '3rd A.A', '4th A.A'],
                                           name='Genus')).round(2)

        df.plot(kind='bar', figsize=(10, 4))

        ax = plt.gca()
        pos = []
        for bar in ax.patches:
            pos.append(bar.get_x() + bar.get_width() / 2.)

        ax.set_xticks(pos, minor=False)
        lab = []
        for i in range(len(pos)):
            l = df.columns.values[i // len(df.index.values)]
            lab.append(l)

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.legend(loc='best')
        ax.set_xticklabels(lab, minor=True)
        ax.tick_params(axis='x', which='major', pad=10, size=0)
        plt.setp(ax.get_xticklabels(), rotation=0)
        plt.title("Biochemical Weights Analysis")
        # plt.figure(figsize=(20, 20))

        fig5 = plt.gcf()
        return fig5

    if k_mears == 5:
        df = pd.DataFrame([factors_1, factors_2, factors_3, factors_4, factors_5, b],
                          index=['Polarity', 'Secondary\nstructure', 'Molecular\nvolume', 'Codon\ndiversity',
                                 'Electrostatic\ncharge', 'Bq | B0'],
                          columns=pd.Index(['1st A.A', '2nd A.A', '3rd A.A', '4th A.A', '5th A.A'],
                                           name='Genus')).round(2)

        df.plot(kind='bar', figsize=(10, 4))

        ax = plt.gca()
        pos = []
        for bar in ax.patches:
            pos.append(bar.get_x() + bar.get_width() / 2.)

        ax.set_xticks(pos, minor=False)
        lab = []
        for i in range(len(pos)):
            l = df.columns.values[i // len(df.index.values)]
            lab.append(l)

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_xticklabels(lab, minor=True)
        ax.tick_params(axis='x', which='major', pad=10, size=0)
        plt.setp(ax.get_xticklabels(), rotation=0)

        # plt.figure(figsize=(20, 20))

        fig5 = plt.gcf()
        return fig5

    if k_mears == 6:
        df = pd.DataFrame([factors_1, factors_2, factors_3, factors_4, factors_5, b],
                          index=['Polarity', 'Secondary\nstructure', 'Molecular\nvolume', 'Codon\ndiversity',
                                 'Electrostatic\ncharge', 'Bq | B0'],
                          columns=pd.Index(
                              ['1st A.A', '2nd A.A', '3rd A.A', '4th A.A', '5th A.A', '6th A.A'],
                              name='Genus')).round(2)

        df.plot(kind='bar', figsize=(10, 4))

        ax = plt.gca()
        pos = []
        for bar in ax.patches:
            pos.append(bar.get_x() + bar.get_width() / 2.)

        ax.set_xticks(pos, minor=False)
        lab = []
        for i in range(len(pos)):
            l = df.columns.values[i // len(df.index.values)]
            lab.append(l)

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Snippet residue positions")
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_xticklabels(lab, minor=True)
        ax.tick_params(axis='x', which='major', pad=10, size=0)
        plt.setp(ax.get_xticklabels(), rotation=0)
        plt.title("Biochemical Weights Analysis")
        # plt.figure(figsize=(20, 20))

        fig5 = plt.gcf()
        return fig5

    if k_mears == 7:
        df = pd.DataFrame([factors_1, factors_2, factors_3, factors_4, factors_5, b],
                          index=['Polarity', 'Secondary\nstructure', 'Molecular\nvolume', 'Codon\ndiversity',
                                 'Electrostatic\ncharge', 'Bq | B0'],
                          columns=pd.Index(
                              ['1st A.A', '2nd A.A', '3rd A.A', '4th A.A', '5th A.A', '6th A.A',
                               '7th A.A'],
                              name='Snippet at:')).round(2)

        df.plot(kind='bar', figsize=(10, 4))

        ax = plt.gca()
        pos = []
        for bar in ax.patches:
            pos.append(bar.get_x() + bar.get_width() / 2.)

        ax.set_xticks(pos, minor=False)
        lab = []
        for i in range(len(pos)):
            l = df.columns.values[i // len(df.index.values)]
            lab.append(l)

        ax.legend(loc='uper left', bbox_to_anchor=(1, 0.5), title="Snippet residue positions")
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_xticklabels(lab, minor=True)
        ax.tick_params(axis='x', which='major', pad=10, size=0)
        plt.setp(ax.get_xticklabels(), rotation=0)
        plt.title("Biochemical Weights Analysis")
        # plt.figure(figsize=(20, 20))

        fig5 = plt.gcf()
        return fig5







