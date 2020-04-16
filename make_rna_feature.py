#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import pandas as pd
import os 
import re 
import glob
import argparse as arg
import sys
import source

class Seqfeature(object):
    ##get miRNA-mRNA (RNA) length
    def getMiRNALen (miRNASeq):
        return len(miRNASeq)

    ##get miRNA-mRNA GC content
    def getGCContent (miRNASeq):
        Len = Seqfeature.getMiRNALen(miRNASeq)
        gc = miRNASeq.count('G')
        gc = gc + miRNASeq.count('g')
        gc = gc + miRNASeq.count('C')
        gc = gc + miRNASeq.count('c')
        GCContent = 1.0*gc/Len
        return "{0:.4f}".format(GCContent)

    ##get mono nucleotide
    def getMonoNucleotideContent (miRNASeq, ch):
        m = miRNASeq.count(ch)
        return 1.0*m/len(miRNASeq)
    
    ##get all mono nucleotide 
    def getAllMonoNucleotideContent( miRNASeq ):
        ContA = Seqfeature.getMonoNucleotideContent( miRNASeq, "A")
        ContC = Seqfeature.getMonoNucleotideContent( miRNASeq, "C")
        ContG = Seqfeature.getMonoNucleotideContent( miRNASeq, "G")
        ContU = Seqfeature.getMonoNucleotideContent( miRNASeq, "U")
        return ["{0:.4f}".format(i) for i in [ContA, ContC, ContG, ContU]]

    ##get dinucleotide content
    def getDiNucleotideContent( miRNASeq, ch1, ch2 ):
        subStr = ch1 + ch2
        m = miRNASeq.count( subStr )
        res = 1.0*m/len(miRNASeq)
        return "{0:.4f}".format(res)

    def getAllDiNucleotideContent( miRNASeq ):
        chs = ["A", "C", "G", "U"]
        res = []
        for i in range(len(chs)):
            ch1 = chs[i]
            for j in range(len(chs)):
                ch2 = chs[j]
                res.append( Seqfeature.getDiNucleotideContent( miRNASeq, ch1, ch2 ) )
        return res

    #Maximum length on the miRNA-mRNA without bulges and its percent
    def getMaximumLengthWithoutBulges(mirnaStruct):
        Len = Max = 0
        for i in range(len(mirnaStruct)):
            ch = mirnaStruct[i]
            if( ch == '(' or ch == ')' ):
                Len = Len + 1
            else:
                Len = 0
            if( Max < Len ):
                Max = Len
        return [Max, "{0:.4f}".format(1.0*Max/len(mirnaStruct))]

    #Base pairs in duplexe miRNA-miRNA*
    def getBasePairsInDuplex( mirnaStruct ):
        p = mirnaStruct.count( '(' ) + mirnaStruct.count( ')' )
        return p

    def isParenthese(ch):
        if( ch=='(' or ch==')' ):
            return True
        else:
            return False

    #number of bulges in miRNA-mRNA
    def getNumberOfBulges( miRNAStruct):
        num = 0
        for i in range(len(miRNAStruct)-2):
            ch1 = miRNAStruct[i]
            ch2 = miRNAStruct[i+1]
            if( ch1 == '.' and Seqfeature.isParenthese(ch2) ):
                num = num + 1 
        return num

    ##20mer base paired
    def getPresenceOfPerfect20MerBasePair(mirnaStruct):
        if( mirnaStruct.find( "((((((((((((((((((((" ) >= 0 or mirnaStruct.find( ")))))))))))))))))))))" ) >= 0 ):
            return 1
        else:
            return 0

    ##start of 20mer base paired
    def getStartOfPerfect20MerBasePair(mirnaStruct):
        if( Seqfeature.getPresenceOfPerfect20MerBasePair(mirnaStruct) == 1 ):
            s = mirnaStruct.find("((((((((((((((((((((")
            if( s == -1 ):
                s = mirnaStruct.find(")))))))))))))))))))))")
            return s+1
        else:
            return -1

    ##10mer base paired
    def getPresenceOfPerfect10MerBasePair(mirnaStruct):
        if( mirnaStruct.find( "((((((((((" ) >= 0 or mirnaStruct.find( "))))))))))" ) >= 0 ):
            return 1
        else:
            return 0

    ##start of 10mer base paired
    def getStartOfPerfect10MerBasePair(mirnaStruct):
        if( Seqfeature.getPresenceOfPerfect10MerBasePair(mirnaStruct) == 1 ):
            s = mirnaStruct.find("((((((((((")
            if( s == -1 ):
                s = mirnaStruct.find("))))))))))")
            return s+1
        else:
            return -1

    ##5mer base paired
    def getPresenceOfPerfect5MerBasePair(mirnaStruct):
        if( mirnaStruct.find( "(((((" ) >= 0 or mirnaStruct.find( ")))))" ) >= 0 ):
            return 1
        else:
            return 0


    ##start of 5mer base paired
    def getStartOfPerfect5MerBasePair(mirnaStruct):
        if( Seqfeature.getPresenceOfPerfect5MerBasePair(mirnaStruct) == 1 ):
            s = mirnaStruct.find("(((((")
            if( s == -1 ):
                s = mirnaStruct.find(")))))")
            return s+1
        else:
            return -1

    ##get average number base pairs in window w
    def getAverageNumberOfBasePairInWindow( mirnaStruct, w ):
        numList = []
        total = 0
        for i in range(len(mirnaStruct) - w):
            subStruct = mirnaStruct[i:i+w+1]
            m1 = subStruct.count( "(")
            m2 = subStruct.count( ")")
            numList.append( m1 + m2 )
            total = total + m1
            total = total + m2
        return "{0:.4f}".format(1.0*total/len(numList))

    ##get all mono-Nucleotide&Structure
    def getMonoSeqStructure( miRNASeq,  miRNAStruct):
        miRNALen = Seqfeature.getMiRNALen( miRNASeq )
        resDic = {"A.":0, "A(":0, "A)":0,"C.":0, "C(":0, "C)":0, "G.":0, "G(":0, "G)":0, "U.":0, "U(":0, "U)":0, "others":0}
        for i in range(miRNALen):
            curKey = miRNASeq[i] + miRNAStruct[i]
            if( curKey in resDic.keys() ):
                resDic[curKey] = resDic[curKey] + 1
            else:
                resDic["others"] = resDic["others"] + 1

        #sort keys
        allKeys = sorted( resDic.keys() )
        res = []
        for curKey in allKeys:
            res.append( "{0:.4f}".format(1.0*resDic[curKey]/miRNALen) )
        return res


    def getDiSeqStructure(miRNASeq, miRNAStruct):
        chars = ["A", "C", "G", "U"]
        dots = ["(", ".", ")"]
        resDic = {}
        for i in range(len(chars)):
            for j in range(len(chars)):
                for x in range(len(dots)):
                    for y in range(len(dots)):
                        curKey = chars[i] + chars[j] + dots[x] + dots[y]
                        resDic[curKey] = 0
        resDic["others_Di"] = 0
    
        for i in range(1, len(miRNASeq)):
            ch1 = miRNASeq[i-1]
            ch2 = miRNASeq[i]
            dot1 = miRNAStruct[i-1]
            dot2 = miRNAStruct[i]
            curKey = ch1 + ch2 + dot1 + dot2
            if( curKey in resDic.keys() ):
                resDic[curKey] = resDic[curKey] + 1
            else:
                resDic["others_Di"] = resDic["others_Di"] + 1

        #sort keys
        allKeys = sorted( resDic.keys() )
        res = []
        for curKey in allKeys:
            res.append( "{0:.4f}".format(1.0*resDic[curKey]/len(miRNASeq)) )
        return res


    def getAllTriplets(miRNASeq, miRNAStruct):
        resDic = {"A...":0,"C...":0,"G...":0,"U...":0,"A(..":0,"C(..":0,"G(..":0,"U(..":0,"A((.":0,"C((.":0,"G((.":0,"U((.":0,"A.((":0,"C.((":0,"G.((":0,"U.((":0,"A(((":0,"C(((":0,"G(((":0,"U(((":0,"A.(.":0,"C.(.":0,"G.(.":0,"U.(.":0,"A..(":0,"C..(":0,"G..(":0,"U..(":0,"A(.(":0,"C(.(":0,"G(.(":0,"U(.(":0,"other_tri":0}
        for i in range(1,len(miRNASeq)-2):
            curKey = miRNASeq[i] + miRNAStruct[i-1:i+2]
            if curKey in resDic.keys():
                resDic[curKey] = resDic[curKey] + 1
            else:
                resDic["other_tri"] = resDic["other_tri"] + 1
    
        #sort keys
        allKeys = sorted( resDic.keys() )
        res = []
        for curKey in allKeys:
            res.append( 1.0*resDic[curKey]/len(miRNASeq) )

        return res



class Features(object):
    def __init__(self):
        self.features=[]
        self.rna_names=[]
        self.feature_names=['miRNALen', 'miRNAGC']
        self.feature_names+=['ContA', 'ContC', 'ContG', 'ContU']
        self.feature_names+= Features._DiNucleotideName()
        self.feature_names+=['mlBulge', 'mlBulge_Perct']
        self.feature_names+=['bpNum']
        self.feature_names+=['numBulges']
        self.feature_names+=['start_20Mer', 'persent_20Mer', 'start_10Mer', 'persent_10Mer', 'start_5Mer', 'persent_5Mer', 'NumBP_Win7', 'NumBP_Win5', 'NumBP_Win3']

        resDic = {"A.":0, "A(":0, "A)":0,"C.":0, "C(":0, "C)":0, "G.":0, "G(":0, "G)":0, "U.":0, "U(":0, "U)":0, "others":0}
        allKeys = sorted( resDic.keys() )
        self.feature_names+=allKeys

        chars = ["A", "C", "G", "U"]
        dots = ["(", ".", ")"]
        resDic = {}
        for i in range(len(chars)):
            for j in range(len(chars)):
                for x in range(len(dots)):
                    for y in range(len(dots)):
                        curKey = chars[i] + chars[j] + dots[x] + dots[y]
                        resDic[curKey] = 0
        resDic["others_Di"] = 0
        allKeys = sorted( resDic.keys() )
        self.feature_names+=allKeys


        resDic = {"A...":0,"C...":0,"G...":0,"U...":0,"A(..":0,"C(..":0,"G(..":0,"U(..":0,"A((.":0,"C((.":0,"G((.":0,"U((.":0,"A.((":0,"C.((":0,"G.((":0,"U.((":0,"A(((":0,"C(((":0,"G(((":0,"U(((":0,"A.(.":0,"C.(.":0,"G.(.":0,"U.(.":0,"A..(":0,"C..(":0,"G..(":0,"U..(":0,"A(.(":0,"C(.(":0,"G(.(":0,"U(.(":0,"other_tri":0}
        allKeys = sorted( resDic.keys() )
        self.feature_names+=allKeys


        self.ft = pd.DataFrame()
        print ("builting features")

    @staticmethod
    def _DiNucleotideName():
        chs = ["A", "C", "G", "U"]
        res=[]
        for i in range(len(chs)):
            ch1 = chs[i]
            for j in range(len(chs)):
                ch2 = chs[j]
                res.append(ch1+ch2)
        return res
    def tf2file(self, fn="featuretable.csv"):
        ft2 = self.ft.copy()
        ft2.insert(loc=0, column='RNA_ID', value=self.rna_names)
        ft2.to_csv(fn,index=False)

    def t2f(self, t=pd.DataFrame()):
        print ("t2f")
        
        for i in t.index:
            curSeq = t.ix[i, 'Target Sequence refined']
            curStruc = t.ix[i, 'secondary_Struc']
            curRNAname = t.ix[i, 'RNA_ID']
            fs = Features.deCodingMiRNASequence(curSeq, curStruc)
            self.features.append(fs)
            self.rna_names.append(curRNAname)
            
        self.ft = pd.DataFrame(self.features, columns=self.feature_names)

    @staticmethod
    def deCodingMiRNASequence(seq,struc):
        #1st feature:miRNA-mRNA length
        miRNALen = Seqfeature.getMiRNALen(seq)
        #2nd feature: miRNA-mRNA GCContent
        miRNAGC = Seqfeature.getGCContent(seq)

        #3-6 feature: content of mono nucleotide 
        MonoNucleotideContent = Seqfeature.getAllMonoNucleotideContent( seq)

        #7-22 features: dinucleotide content 
        DiNucleotideContent = Seqfeature.getAllDiNucleotideContent( seq)

        #feature: maximal length of miRNA-mRNA without bulge, mlBulge percent, 
        [mlBulge, mlBulge_Perct] = Seqfeature.getMaximumLengthWithoutBulges( struc)

        #feature: number of base pairs in miRNA-mRNA duplex
        bpNum = Seqfeature.getBasePairsInDuplex( struc)

        #skip dist2loop, dist2Helix here, since it is cmp between miRNA and preMiRNA
        #1 feature: number of bulges in miRNA sequence
        numBulges = Seqfeature.getNumberOfBulges(  struc)

        start_20Mer = Seqfeature.getStartOfPerfect20MerBasePair(struc)
        persent_20Mer = Seqfeature.getPresenceOfPerfect20MerBasePair(struc)
        start_10Mer = Seqfeature.getStartOfPerfect10MerBasePair(struc)
        persent_10Mer = Seqfeature.getPresenceOfPerfect10MerBasePair(struc)
        start_5Mer = Seqfeature.getStartOfPerfect5MerBasePair(struc)
        persent_5Mer = Seqfeature.getPresenceOfPerfect5MerBasePair(struc)
        NumBP_Win7 = Seqfeature.getAverageNumberOfBasePairInWindow( struc, 7 )
        NumBP_Win5 = Seqfeature.getAverageNumberOfBasePairInWindow( struc, 5 )
        NumBP_Win3 = Seqfeature.getAverageNumberOfBasePairInWindow( struc, 3 )

        #12 features: mono-nucleotide & strucutre features
        MonoSeqStructFeatures = Seqfeature.getMonoSeqStructure( seq, struc)

        # features: di-nucleotide and structure features
        DiSeqStructFeatures = Seqfeature.getDiSeqStructure( seq, struc)

        #33 features: triplets sequence & strucutre features
        tripleSeqStructFeatures = Seqfeature.getAllTriplets(seq, struc)






        return([miRNALen, miRNAGC]+MonoNucleotideContent + DiNucleotideContent + [mlBulge, mlBulge_Perct] + [bpNum] +[numBulges] + [start_20Mer, persent_20Mer, start_10Mer, persent_10Mer,start_5Mer, persent_5Mer, NumBP_Win7, NumBP_Win5, NumBP_Win3]+MonoSeqStructFeatures + DiSeqStructFeatures + tripleSeqStructFeatures)



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
    t.drop_duplicates(subset=['RNA_ID'], inplace=True)
    print ('After drop duplicate RNA_ID length: {0}'.format(len(t))) 
    t.index=range(0, len(t))

    f = Features()
    f.t2f(t=t)
    f.tf2file(fn=args.output_file)

