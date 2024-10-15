# CalcTanimoto.py
# Author: Illimar Rekand
# email: illimar.rekand@uib.no
# Date: 2024-10-15 

import RDKit
from RDKit import Chem, DataStructs
import csv


mol1 = Chem.MolFromSmiles("ccn")
mol2 = Chem.MolFromSmiles("ccc")

with open('template.csv')