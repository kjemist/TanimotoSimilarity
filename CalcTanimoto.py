# CalcTanimoto.py
# Author: Illimar Rekand
# email: illimar.rekand@uib.no
# Date: 2024-10-15

'''
usage: CalcTanimoto.py 
'''

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
import csv
from operator import itemgetter

def CreateIterableLists(input_df):
    smiles_list = input_df["SMILES"].to_list()
    mol_list = [Chem.MolFromSmiles(x) for x in smiles_list]
    fingerprint_list = [FingerprintMols.FingerprintMol(x) for x in mol_list] # Convert MolFiles to fingerprints in order to calculate Tanimoto score
    return fingerprint_list, smiles_list

smiles_template_df =  pd.read_csv('./template_w_scaffold.csv', delimiter=";")
print(smiles_template_df)
smiles_comparison_df =  pd.read_csv('./comparison_set_w_scaffold.csv', delimiter=";")

fp_template, smiles_template = CreateIterableLists(smiles_template_df)
fp_compset, smiles_comparison  = CreateIterableLists(smiles_comparison_df)


for tempmol, tempsmiles in zip(fp_template, smiles_template):
    rows = []
    for compmol, compsmiles in zip(fp_compset, smiles_comparison):
        row = []
        tanimoto_score = DataStructs.TanimotoSimilarity(tempmol, compmol)
        print(f'{tempsmiles}-{compsmiles}', tanimoto_score)
        row.append(tempsmiles + "-" + compsmiles)
        row.append(tanimoto_score)
        rows.append(row)
    
    file_name = "tanimoto_scores_" + str(fp_template.index(tempmol)+1) + "_w_scaffold.csv"
    rows_sorted = sorted(rows, key = itemgetter(1), reverse=True) #sort rows by highest score on top
    rows_sorted.insert(0,["SMILES pair","Tanimoto Score"])
    with open(file_name, "w", newline="\n") as file:
        writer = csv.writer(file, delimiter= ",")
        writer.writerows(rows_sorted)