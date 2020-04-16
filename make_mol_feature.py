#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import pandas as pd
import os 
import re 
import glob
import argparse as arg
import sys
import source
import numpy as np

from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.AtomPairs.Sheridan import GetBPFingerprint
from rdkit.Chem.EState.Fingerprinter import FingerprintMol
#from rdkit.Avalon.pyAvalonTools import GetAvalonFP #GetAvalonCountFP  #int vector version
from rdkit.Chem.AllChem import  GetMorganFingerprintAsBitVect, GetErGFingerprint
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import rdkit.DataStructs.cDataStructs

class molutil():
    def s2m(m):
        try:
            mol = Chem.MolFromSmiles(m)
            if mol:
                return mol
            else:
                return "error"
        except:
            return "error"

    def ExplicitBitVect_to_NumpyArray(bitvector):
        bitstring = bitvector.ToBitString()
        intmap = map(int, bitstring)
        return np.array(list(intmap))

class fingerprint():
    def __init__(self, fp_fun, name):
        self.fp_fun = fp_fun
        self.name = name
        self.x = []

    def apply_fp(self, mols):
        for mol in mols:
            #print ("****************** ", mol , " ***************************")
            fp = self.fp_fun(mol)
            if isinstance(fp, tuple):
                fp = np.array(list(fp[0]))
            if isinstance(fp, rdkit.DataStructs.cDataStructs.ExplicitBitVect):
                fp = molutil.ExplicitBitVect_to_NumpyArray(fp)
            if isinstance(fp,rdkit.DataStructs.cDataStructs.IntSparseIntVect):
                fp = np.array(list(fp))

            self.x += [fp]

            if (str(type(self.x[0])) != "<class 'numpy.ndarray'>"):
                print("WARNING: type for ", self.name, "is ", type(self.x[0]))




if __name__ == '__main__':
    opt=arg.ArgumentParser(description='Create dose response upload file from csv')
    opt.add_argument('-i','--input_file', help='input file')
    opt.add_argument('-o','--output_file', type=str, help='output file')
    opt.add_argument('--length', type=int, default=1024, help='nBit length when builting fingerprint')
    opt.add_argument('--fp', default='RDKit', type=str, choices=['RDKit', 'rdkit', 'Morgan', 'morgan'], help='Fingerprint name')
    args=opt.parse_args()

    #read in input file
    try:
        t= pd.read_csv(args.input_file, dtype = str)
    except:
        t= pd.read_csv(args.input_file, dtype = str, encoding='ISO-8859-1')
    print ('Original file length: {0}'.format(len(t))) 

    t.drop_duplicates(subset=['Ligand_ID'], inplace=True)
    print ('after drop duplicates: {0}'.format(len(t))) 
    t['Mol'] = t['SMILES'].apply(molutil.s2m)
    t= t[t["Mol"]!=None]
    print ('after convert mol length: {0}'.format(len(t))) 

    if args.fp == 'RDKit' or args.fp == 'rdkit':
        fp = fingerprint(lambda x : RDKFingerprint(x, fpSize=args.length), "RDKit fingerprint")
    elif args.fp == 'Morgan' or args.fp == 'morgan':
        fp = fingerprint(lambda x : GetMorganFingerprintAsBitVect(x, 2, nBits = args.length), "Morgan circular")

    mol_id = t['Ligand_ID'].values
    fp.apply_fp(list(t['Mol']))
    print ('Size of the fingerprint {0}'.format(len(fp.x), len(fp.x[0])))

    col_names=[]
    for i in range(len(fp.x[0])):
        col_names.append('col'+str(i)) 
    ft = pd.DataFrame(fp.x, columns=col_names)
    ft.insert(loc=0, column='Ligand_ID', value=mol_id)
    print (len(ft))
    ft.to_csv(args.output_file,index=False)


