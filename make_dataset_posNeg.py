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
from rdkit.Chem.AllChem import  GetMorganFingerprintAsBitVect, GetErGFingerprint
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import rdkit.DataStructs.cDataStructs


#create positive and negative datasets, matched rnaid and ligandid are positive, other ranmomized are negative sets

if __name__ == '__main__':
    opt=arg.ArgumentParser(description='Create dose response upload file from csv')
    opt.add_argument('-i','--input_file', help='input file')
    opt.add_argument('--rf',help='rna feature file')
    opt.add_argument('--lf',help='ligand feature file')
    opt.add_argument('-o','--output_file', type=str, help='output file')
    args=opt.parse_args()

    #read in input file
    try:
        t= pd.read_csv(args.input_file, dtype = str)
    except:
        t= pd.read_csv(args.input_file, dtype = str, encoding='ISO-8859-1')
    print ('Original file length: {0}'.format(len(t))) 
    t.drop_duplicates(subset=['Ligand_ID','RNA_ID'], inplace=True)
    print ('Original file after drop duplicates length: {0}'.format(len(t))) 
    ligandIDs = list(set(t['Ligand_ID'])) 
    rnaIDs = list(set(t['RNA_ID']))
    print ('Unique Ligand IDs: {0}, RNA IDs: {1}'.format(len(ligandIDs), len(rnaIDs)))
    t['Pair'] = t['RNA_ID']+":"+t['Ligand_ID']

    #get all positive pairs
    t_pos = t[['Pair']].copy()
    t_pos['status'] = 1
    print ('Positive pair length: {0}'.format(len(t_pos)) )

    #get all negative pairs
    rls = [[ r+":"+l, 0 ] for r in rnaIDs for l in ligandIDs]
    t_all = pd.DataFrame(rls, columns=['Pair','status'] )
    print ('all possible pair length: {0}'.format(len(t_all)) )
    t_all = t_all[~t_all['Pair'].isin(t_pos['Pair'])] 
    t_all.sort_values(by=['Pair'], inplace=True)
    t_all.index = range(0, len(t_all))
    print ('all possible pair length: {0}'.format(len(t_all)) )
    t_neg_sample = t_all.sample(n=3*len(t_pos), random_state=1000)
    t.drop_duplicates(inplace=True)
    print ('all possible pair length after sample: {0}'.format(len(t_neg_sample)) )
    #print (t_neg_sample[0:10]) 

    #join positive and neg_sample for data sets
    t_data= pd.concat([t_pos, t_neg_sample], ignore_index=True)
    t_data.index = range(0, len(t_data))

    t_data[['RNA_ID','Ligand_ID']] = t_data['Pair'].str.split(":",expand=True)
    
    #merge with features
    try:
        t_rna= pd.read_csv(args.rf, dtype = str)
    except:
        t_rna= pd.read_csv(args.rf, dtype = str, encoding='ISO-8859-1')
    print ('rna feature lngth: {0}'.format(len(t_rna)) )
    #t=pd.merge(t, t_hcs[cols], left_on=['abarcode', 'row', 'col'], right_on=[pIDHeader, 'row','col'], how = 'left', suffixes=('',merge_suffix))
    t_data=pd.merge(t_data, t_rna, left_on=['RNA_ID'], right_on=['RNA_ID'], how = 'left')
    print ('all data after rna feature merge length: {0}'.format(len(t_data)) )

    try:
        t_ligd= pd.read_csv(args.lf, dtype = str)
    except:
        t_ligd= pd.read_csv(args.lf, dtype = str, encoding='ISO-8859-1')
    print ('ligand feature lngth: {0}'.format(len(t_ligd)) )
    t_data=pd.merge(t_data, t_ligd, left_on=['Ligand_ID'], right_on=['Ligand_ID'], how = 'left')
    t_data.drop(['RNA_ID', 'Ligand_ID'], axis=1, inplace=True)
    print ('all data after ligand feature merge length: {0}'.format(len(t_data)) )

    #print (t_data[0:10])
    #print (t_data[len(t_data)-10:])

    #sample and split to train and test
    t_train = t_data.sample(frac=0.9, random_state=1000)
    t_test = t_data[~t_data['Pair'].isin(t_train['Pair'])] 
    print (len(t_train), len(t_test), len(t_data))
    output_train = args.output_file +'.train.csv'
    output_test = args.output_file +'.test.csv'
    t_train.to_csv(output_train,index=False)
    t_test.to_csv(output_test,index=False)

