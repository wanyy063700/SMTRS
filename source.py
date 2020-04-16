from rdkit import Chem
import re
import os

RNAfoldDic = "/depts/ChemInfo/yzhai/scripts/ViennaRNA/bin/"

tmpFileDir = "tmp_seq.txt"
RNAfoldResDir = "tmp_struct.txt"

def dummyfun():
    pass

def s2m(m):
    try:
        mol = Chem.MolFromSmiles(m)
        if mol:
            return mol
        else:
          return "error"
    except:
        return "error"


def checkseq(s):
    if not s or type(s) != str:
        return "error" 

    if re.search(r"\/", s):
        return "error" 
    elif re.search(r"-", s):
        return "error" 
    elif re.search(r"^\s*$", s):
        return "error" 
    else:
        return "OK" 

def readLinesFromFile (fileDir):
    fpr = open(fileDir, "r")
    lines = fpr.readlines()
    fpr.close()
    return lines

##parse RNAfold result
def parseRNAfoldResult( RNAfoldResDir ):
    lines = readLinesFromFile(RNAfoldResDir)
    curStruct = lines[2].split(" ")[0]
    return curStruct

'''
  for ii in range(seqNum):
    if( (6*ii + 5) > len(lines) ):
      break
    curID = lines[6*ii]
    curID = curID[1:]
    curID = curID.strip()
    curStruct = lines[6*ii+2]
    elements = curStruct.split(" ")
    resDic[curID] = elements[0]
  #end for ii
  return resDic
'''

def foldRNASeq(rid, s):
    fpw = open(tmpFileDir, "w")
    fpw.write(">" + rid + "\n" +  s + "\n" )
    fpw.close()

    foldCommand = RNAfoldDic + "RNAfold -p -d2 --noLP < " + tmpFileDir + " > " + RNAfoldResDir
    res = os.system(foldCommand)    

    if( res != 0 ):
        print ("Error to run RNAfold using the command:", rid, foldCommand)
        return "error" 
        #sys.exit()
    else:
        os.remove(rid+'_dp.ps') #RSMDR_5_dp.ps
        os.remove(rid+'_ss.ps') #RSMDR_6_ss.ps
        curStruct = parseRNAfoldResult(RNAfoldResDir)
        if curStruct and re.search(r"[\.\(\)]", curStruct):
            return curStruct
        else: 
            return "error" 
        #print (curStruct)








