import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
import numpy as np
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold
from chembl_webresource_client.new_client import new_client
import pandas as pd
#from flask_table import Table, Col
# from IPython.display import display
# from tabulate import tabulate


# class ItemTable(Table):
#     cross_references = Col('cross_references')
#     organism = Col('organism')
#     pref_name = Col('pref_name')
#     score = Col('score')
#     species_group_flag = Col('species_group_flag')
#     target_chembl_id = Col('target_chembl_id')
#     target_components = Col('target_componentss')
#     target_type = Col('target_type')
#     tax_id = Col('tax_id')


def searchQuery(input_value):
    target = new_client.target
    search_query = target.search(input_value)

    moleculeTargets = pd.DataFrame.from_dict(search_query)
    return moleculeTargets

    # table = ItemTable(search_query)
    # return table


def selectTarget(n__o):
    activity = new_client.activity
    activities = activity.filter(
        target_chembl_id=n__o).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(activities)
    df2 = df[df.standard_value.notna()]
    df2.to_csv('pic50_data_preparation.csv', index=False)

    return df2


def phase_i():
    df = pd.read_csv('pic50_data_preparation.csv')

    molecule_bioactivity = []
    for i in df.standard_value:
        if float(i) >= 10000:
            molecule_bioactivity.append("inactive")
        elif float(i) <= 1000:
            molecule_bioactivity.append("active")
        else:
            molecule_bioactivity.append("intermediate")

    mol_cid = []
    for i in df.molecule_chembl_id:
        mol_cid.append(i)

    canonical_smiles = []
    for i in df.canonical_smiles:
        canonical_smiles.append(i)

    standard_value = []
    for i in df.standard_value:
        standard_value.append(i)

    data_tuples = list(zip(mol_cid, canonical_smiles,
                       molecule_bioactivity, standard_value))
    df3 = pd.DataFrame(data_tuples, columns=[
                       'mol_cid', 'canonical_smiles', 'molecule_bioactivity', 'standard_value'])

    # Droping the empty rows
    df4 = df3.dropna()

    # Saving it to the csv file
    df4.to_csv('ssri_bioactivity.csv', index=False)

    return df4


# def selectTarget(ev, n__o):
#     target = new_client.target
#     search_query = target.search(ev)
#     moleculeTargets = pd.DataFrame.from_dict(search_query)

#     selectedTarget = moleculeTargets.target_chembl_id[n__o]
#     return selectedTarget

# moleculeTargets = pd.DataFrame.from_dict(search_query).to_html()
# # write html to file
# moleculeTargets_html_file = open("moleculeTargets.html", "w")
# moleculeTargets_html_file.write(moleculeTargets)
# moleculeTargets_html_file.close()

# # open html file
# file = codecs.open("moleculeTargets.html", "r", "utf-8")

# return file.read()

# return (tabulate(moleculeTargets, headers='keys', tablefmt='psql'))
# moleculeTargets = pd.DataFrame.from_dict(search_query)
# return HTML(moleculeTargets.to_html(classes='table table-striped'))
# mt_display = display(moleculeTargets)
# return display(moleculeTargets)

# selectedTarget = moleculeTargets.target_chembl_id[0]
# selectedTarget

#     activity = new_client.activity
#     activities = activity.filter(target_chembl_id=selectedTarget).filter(standard_type="IC50")
#     df = pd.DataFrame.from_dict(activities)

#     dataframe2 = df[df.standard_value.notna()]
#     dataframe2

# dataframe2.to_csv('data_preparation.csv', index=False)
