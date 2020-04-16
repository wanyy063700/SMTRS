#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import pandas as pd
#import util
import os 
import re 
#import db 
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


#using one RNA ID and a RNA feature file and a ligand feature file to make feature file for the pair 


if __name__ == '__main__':
    opt=arg.ArgumentParser(description='Create dose response upload file from csv')
    opt.add_argument('-i','--input_file', help='input file')
    opt.add_argument('--rnaid',help='rna id file')
    opt.add_argument('--rf',help='rna feature file')
    opt.add_argument('--lf',help='ligand feature file')
    opt.add_argument('-o','--output_file', type=str, help='output file')
    args=opt.parse_args()

    
    t_rna= pd.read_csv(args.rf, dtype = str)
    t_ligd = pd.read_csv(args.lf, dtype = str)
    ligandIDs = list(set(t_ligd['Ligand_ID']))
    pairs = [[args.rnaid, lid, args.rnaid+":"+lid, ''] for lid in ligandIDs]
    t_data= pd.DataFrame(pairs, columns=['RNA_ID', 'Ligand_ID', 'Pair','status'] )
    print ('Positive pair length: {0}'.format(len(t_data)) )

    t_data=pd.merge(t_data, t_rna, left_on=['RNA_ID'], right_on=['RNA_ID'], how = 'left')
    print ('all data after rna feature merge length: {0}'.format(len(t_data)) )

    t_data=pd.merge(t_data, t_ligd, left_on=['Ligand_ID'], right_on=['Ligand_ID'], how = 'left')
    t_data.drop(['RNA_ID', 'Ligand_ID'], axis=1, inplace=True)
    print ('all data after ligand feature merge length: {0}'.format(len(t_data)) )
    #print (t_data)
    t_data.to_csv(args.output_file,index=False)
