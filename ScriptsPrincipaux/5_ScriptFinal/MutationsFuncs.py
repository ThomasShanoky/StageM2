import pandas as pd



def getTypesOfMutationsAndInd(Genes, gene, data_mutation):
    """Retourne un dictionnaire avec les types de mutations et les individus associés"""
    
    geneMutations = data_mutation[Genes.index(gene)]
    dico_mut = {}
    dico_IndAndMut = {}
    for _, row in geneMutations.iterrows():
        Id = row["sampleID"]
        local = row["localisation"]
        alt = row["alt"]
        ref = row["ref"]

        if [local, alt, ref] not in dico_mut.values():
            type_mut = "mut" + str(len(dico_mut)+1)
            dico_mut[type_mut] = [local, alt, ref]
            dico_IndAndMut[type_mut] = []

        dico_IndAndMut[type_mut].append(Id[:6])

    return dico_IndAndMut, dico_mut


def format_mutations(dico_mut, dico_IndAndMut, minIndThreshold=4):
    formatted_mutations = ["Toutes les mutations"]
    for mut, details in dico_mut.items():
        count = len(dico_IndAndMut[mut])
        if count >= minIndThreshold:
            formatted_mutations.append(f"{mut}: {details[2]}>{details[1]} ({count} patients)")
    return formatted_mutations


def checkIfPatientsAreInExpressions(data_expressions, dico_IndAndMut):
    """Retourne la liste d'individus à supprimer (qui ne sont pas dans les données d'expressions)"""
    Ind = dico_IndAndMut.values()
    Ind = [item for sublist in Ind for item in sublist]

    IndToRemove = []
    for ind in Ind:
        if ind+"R" not in data_expressions.columns:
            IndToRemove.append(ind)

    return IndToRemove


def filterDicoIndAndMut(dico_IndAndMut, minIndThreshold=4):
    """Supprime les mutations qui n'ont pas assez d'individus associés"""

    mutToDel = []
    for mut in dico_IndAndMut:
        if len(dico_IndAndMut[mut]) < minIndThreshold:
            mutToDel.append(mut)

    dico_IndAndMut = {key: value for key, value in dico_IndAndMut.items() if key not in mutToDel}

    return dico_IndAndMut


def getTable(data_expressions, gene, dico_IndAndMut):
    """Retourne un tableau avec 3 colonnes : l'ID sample, le type de mutation et l'expression du gène"""
    Tableau = pd.DataFrame(columns=["ID_sample", "TypeMutation", "ExpressionGene"])

    for mut in dico_IndAndMut.keys():
        for Id in dico_IndAndMut[mut]:
            Expression = data_expressions.loc[data_expressions['display_label'] == gene, Id + "R"].values[0]
            Tableau = Tableau._append({"ID_sample": Id, "TypeMutation": mut, "ExpressionGene": Expression}, ignore_index=True)

    return Tableau


def getNonMutatedInd(data_expressions, AllIndMutated):
    NonMutatedInd = []
    AllInd = list(data_expressions.columns[4:]) #Tous les individus présents dans le fichier d'expressions géniques
    AllInd = [ind[:6] for ind in AllInd]
    
    for ind in AllInd:
        if ind not in AllIndMutated:
            NonMutatedInd.append(ind)
    return NonMutatedInd


def checkIfIndExpressionAreInIndex(IndIndex, NonMutatedInd):
    IndToDel = []

    for ind in NonMutatedInd:
        if ind not in IndIndex:
            IndToDel.append(ind)

    new_NonMutatedInd = [ind for ind in NonMutatedInd if ind not in IndToDel]

    return new_NonMutatedInd


def addNonMutatedIndExpression(data_expressions, Tableau, gene, nonMutatedInd):

    for ind in nonMutatedInd:
        Expression = data_expressions.loc[data_expressions['display_label'] == gene, ind].values[0]
        Tableau = Tableau._append({"ID_sample": ind, "TypeMutation": "NonMut", "ExpressionGene": Expression}, ignore_index=True)

    return Tableau