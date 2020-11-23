import pickle
from statistics import mean

def recover_variable():
    # Ouverture du fichier si il existe et récupération de la list
    fichierSauvegarde = open("C:/Users/Jonas/PycharmProjects/ProjectLab/datasetCeliac/cleanData", "rb")
    curve = pickle.load(fichierSauvegarde)
    fichierSauvegarde.close()
    return curve


def display_metaData():
    # recover_result
    repertoires = recover_variable()
    repertoire_len = []
    for repertoire in repertoires:
        repertoire_len.append(int(repertoire["size_patient"]))
    repertoire_len.sort(reverse=True)
    length = len(repertoire_len)
    max = repertoire_len[0]
    min = repertoire_len[len(repertoire_len)-1]
    meann = int(mean(repertoire_len))
    sick = 46
    unsick = 40
    numOfMEn = 18
    numOfWomen = 45
    return {"length": length, "max": max, "min": min, "mean": meann, "sick": sick, "unsick": unsick, "numOfMen": numOfMEn, "numOfWomen": numOfWomen}
