#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import pandas as pd
import os 
import re 
import glob
import argparse as arg
import sys
import source


if __name__ == '__main__':
    opt=arg.ArgumentParser(description='Create dose response upload file from csv')
    opt.add_argument('-i','--input_file', help='input file')
    opt.add_argument('-o','--output_file', type=str, help='output file')
    args=opt.parse_args()

    #read in input file
    try:
        t= pd.read_csv(args.input_file, dtype = str)
    except:
        t= pd.read_csv(args.input_file, dtype = str, encoding='ISO-8859-1')
    print ('Original file length: {0}'.format(len(t))) 

    if 'SMILES' in t.columns:
        t['Mol'] = t['SMILES'].apply(source.s2m)
        t= t[t["Mol"]!='error']
        print ('file after molecule check length: {0}'.format(len(t))) 
        t=t.drop(['Mol'], axis=1)

    if 'Target Sequence refined' in t.columns:
        t['seq_check'] = t['Target Sequence refined'].apply(source.checkseq)
        t= t[t["seq_check"]!='error']
        print ('file after sequence check length: {0}'.format(len(t))) 

    if 'Target Sequence refined' in t.columns:
        #fitted['Half Life (min)'] = fitted.apply(lambda x: Prodrug.half_life(x['Conversion_120'], x['Conversion-Media_T120'], x['T0.5']), axis=1)
        t['secondary_Struc'] = t.apply(lambda x: source.foldRNASeq(x['RNA_ID'],x['Target Sequence refined']), axis=1)
        print ('file after secondary structure prediction, length: {0}'.format(len(t))) 

        t= t[t["secondary_Struc"]!='error']

    t.to_csv(args.output_file,index=False)



