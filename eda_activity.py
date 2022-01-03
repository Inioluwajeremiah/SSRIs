import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import seed
from numpy.random import randn
import seaborn as sns
from scipy.stats import mannwhitneyu

sns.set(style='ticks')


def phaseii_Lipinski():
    df = pd.read_csv('ssri_bioactivity.csv')

    df_lipinski = generateDescriptors(df.canonical_smiles)
    dfCombined = pd.concat([df, df_lipinski], axis=1)
    df_norm = norm_value(dfCombined)
    df_final = convert_to_pIC50(df_norm)
    df_rem_int_class = df_final[df_final.molecule_bioactivity !=
                                'intermediate']
    df_final.to_csv('lipinski_bioactivity.csv', index=False)

    return df_final


# generate lipinski descriptors from the smiles notation

def generateDescriptors(smiles, verbose=False):

    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHAcceptors = Descriptors.NumHAcceptors(mol)
        desc_NumHDonors = Descriptors.NumHDonors(mol)

        row = np.array([desc_MolWt, desc_MolLogP,
                       desc_NumHAcceptors, desc_NumHDonors])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["MolWt", "MolLogP", "NumHAcceptors", "NumHDonors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

# print("Total descriptors generated: 17 x "+str(len(smiles)))
    return descriptors


def convert_to_pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        #         convert from nM to M
        molar = i*(10**-6)
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
    return x


# normalize funtion.

def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 20000:
            i = 20000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)
    return x


def eda_bioactivity_plot():
    bioactivity_df = pd.read_csv('lipinski_bioactivity.csv')

    plt.figure(figsize=(5.5, 5.5))

    sns.countplot(x='molecule_bioactivity',
                  data=bioactivity_df, edgecolor='black')

    plt.xlabel('Bioactivity', fontsize=14, fontWeight='bold')
    plt.ylabel('Frequency', fontSize=14, fontWeight='bold')
    plt.savefig('plot_bioactivity.csv')
    plt.savefig('plot_bioactivity.png')
    bioactivity_plot = plt.savefig('plot_bioactivity.png')
    return bioactivity_plot


def eda_sp_mw_vs_logp():
    bioactivity_df = pd.read_csv('lipinski_bioactivity.csv')
    plt.figure(figsize=(5.5, 5.5))
    sns.scatterplot(x='MolWt', y='MolLogP', data=bioactivity_df,
                    hue='molecule_bioactivity', size='pIC50', edgeColor='black', alpha=0.7)
    plt.xlabel('MW', fontsize=14, fontWeight='bold')
    plt.ylabel('LogP', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.savefig('plot_W_vs_LogP.pdf')
    mv_vs_log_img = plt.savefig('plot_W_vs_LogP.png')
    return mv_vs_log_img


def pIC50_boxPlot():
    df_final = pd.read_csv('lipinski_bioactivity.csv')
    plt.figure(figsize=(5.5, 5.5))
    sns.boxplot(x='molecule_bioactivity', y='pIC50', data=df_final)
    plt.xlabel('Bioactivity value', fontsize=14, fontWeight='bold')
    plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')

    plt.savefig('plot_IC50.pdf')
    pic50_box_img = plt.savefig('plot_IC50.png')
    man_PIC50 = mannWhitney('pIC50')

    pIC50_data = man_PIC50
    return pIC50_data


def molwt_boxplot():
    df_final = pd.read_csv('lipinski_bioactivity.csv')
    plt.figure(figsize=(5.5, 5.5))

    sns.boxplot(x='molecule_bioactivity', y='MolWt', data=df_final)
    plt.xlabel('Bioactivity value', fontsize=14, fontWeight='bold')
    plt.ylabel('MW', fontsize=14, fontweight='bold')

    plt.savefig('plot_MW.pdf')

    molwt_box_img = plt.savefig('plot_molwt.png')
    man_molwt = mannWhitney('MolWt')

    molwt_data = [molwt_box_img, man_molwt]
    return molwt_data


def logp_boxplot():
    df_final = pd.read_csv('lipinski_bioactivity.csv')
    plt.figure(figsize=(5.5, 5.5))

    sns.boxplot(x='molecule_bioactivity', y='MolLogP', data=df_final)
    plt.xlabel('Bioactivity value', fontsize=14, fontWeight='bold')
    plt.ylabel('LogP', fontsize=14, fontweight='bold')

    plt.savefig('plot_LogP.pdf')
    LogP_box_img = plt.savefig('plot_LogP.png')
    man_LogP = mannWhitney('MolLogP')

    LogP_data = [LogP_box_img, man_LogP]
    return LogP_data


def hdonor_boxplot():
    df_final = pd.read_csv('lipinski_bioactivity.csv')
    plt.figure(figsize=(5.5, 5.5))

    sns.boxplot(x='molecule_bioactivity', y='NumHDonors', data=df_final)
    plt.xlabel('Bioactivity value', fontsize=14, fontWeight='bold')
    plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')

    plt.savefig('plot_num_hydrogen_donors.pdf')
    hdonor_box_img = plt.savefig('plot_hdonor.png')
    man_hdonor = mannWhitney('NumHDonors')

    hdonor_data = [hdonor_box_img, man_hdonor]
    return hdonor_data


def hacceptor_boxplot():
    df_final = pd.read_csv('lipinski_bioactivity.csv')
    plt.figure(figsize=(5.5, 5.5))

    sns.boxplot(x='molecule_bioactivity', y='NumHAcceptors', data=df_final)
    plt.xlabel('Bioactivity value', fontsize=14, fontWeight='bold')
    plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')

    plt.savefig('plot_num_hydrogen_acceptors.pdf')
    hacceptor_box_img = plt.savefig('plot_hacceptor.png')
    man_hacceptor = mannWhitney('NumHAcceptors')

    hacceptor_data = [hacceptor_box_img, man_hacceptor]
    return hacceptor_data


def mannWhitney(descriptor, verbose=False):

    # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python

    # seed the random number generator
    seed(1)

    # actives and inactives
    df_final = pd.read_csv('lipinski_bioactivity.csv')
    df_rem_int_class = df_final[df_final.molecule_bioactivity !=
                                'intermediate']

    selection = [descriptor, 'molecule_bioactivity']
    df = df_rem_int_class[selection]
    active = df[df.molecule_bioactivity == 'active']
    active = active[descriptor]

    selection = [descriptor, 'molecule_bioactivity']
    df = df_rem_int_class[selection]
    inactive = df[df.molecule_bioactivity == 'inactive']
    inactive = inactive[descriptor]

    #     compare samples
    stat, p = mannwhitneyu(active, inactive)
    # print('statistics%.3f, p=%.3f' %(stat, p))
    # interpret
    # print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
    results = pd.DataFrame({'Descriptor': descriptor, 'Statistics': stat, 'P': p,
                            'alpha': alpha, 'Interpretation': interpretation}, index=[0])
    filename = 'manwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)
    return results
