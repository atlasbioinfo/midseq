'''
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
                Copyright 2019 Haopeng Yu atlasbioin4@gmail.com
        Licensed under the Apache License, Version 2.0 (the "License");
        you may not use this file except in compliance with the License.
'''     
import sys,re
import argparse

baseT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "U": "A",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
    "u": "a",
    "N": "N",
    "n": "n",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "W": "W",
    "S": "S",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "r": "y",
    "y": "r",
    "m": "k",
    "k": "m",
    "w": "w",
    "s": "s",
    "h": "d",
    "d": "h",
    "b": "v",
    "v": "b"
}

def tranStrand(seq):
    rSeq = seq[::-1]
    newSeq = ""
    for i in range(len(rSeq)):
        newSeq += baseT[rSeq[i]]
    return newSeq

parser = argparse.ArgumentParser(description='Get the internal sequence between 5\' and 3\' boundary sequence from fastq file. E.g., get string "ATLAS" from "GGGATLASAAA" with -f "GGG" and -t "AAA".')
parser.add_argument("Fastq",nargs='*',metavar='Fastq files',help="HTS reads in fastq format (not in .gz format). Support muti-fastq-files.",type=str)
parser.add_argument("-f",help="The 5' sequence. Default is: CCGGCGACGTTGGGTCAACT", default="CCGGCGACGTTGGGTCAACT",type=str)
parser.add_argument("-t",help="The 3' sequence. Default is: TGTCCTCTTCCTCTTTAGCG", default="TGTCCTCTTCCTCTTTAGCG",type=str)
parser.add_argument("-o","--output",help="The output file name.",default="output_stat.tsv", type=str)
parser.add_argument("-of","--outputfasta",help="The output fasta name. The default is not to output this file unless assigned its file name.", type=str)

args = parser.parse_args()

five=args.f
three=args.t
outputStat=open(args.output,"w")
outputStat.write("#Seq\tlength\tSeq_num\n")
fastaTrigger=0
if (args.outputfasta):
    outputFa=open(args.outputfasta,"w")
    fastaTrigger=1

regFive=re.compile(five,re.I)
regThree=re.compile(three,re.I)
# 负链搜索参数 通常用不到
# regMFive=re.compile(tranStrand(three),re.I)
# regMThree=re.compile(tranStrand(five),re.I)

header=""
lab=0
seqDict={}
for f in range(len(args.Fastq)):
    print("Checking file: "+args.Fastq[f])
    count=1
    fileLines=0
    print("Initial statistics...")
    with open(args.Fastq[f],"r") as f1:
        for line1 in f1:
            fileLines+=1
    if (fileLines % 4 !=0):
        print("The file: "+args.Fastq[f]+" may not a fastq file...")
    allReads=int(fileLines/4)
    print("Processing...")
    with open(args.Fastq[f],"r") as f2:
        for line2 in f2:
            if (re.match(r'^@',line2)):
                header=line2.strip()
                lab=1
            else:
                if(lab):
                    count+=1
                    tmp=line2.strip()
                    fiveMat=regFive.search(tmp)
                    threeMat=regThree.search(tmp)
                    if (fiveMat and threeMat):
                        if (fiveMat.span()[1]<threeMat.span()[0]):
                            tseq=tmp[fiveMat.span()[1]:threeMat.span()[0]]
                            if( fastaTrigger):
                                outputFa.write(">"+header+"\n"+tseq+"\n")
                            if tseq in seqDict:
                                seqDict[tseq]+=1
                            else:
                                seqDict[tseq]=1
                    # 负链搜索，一般用不到
                    # mfiveMat=regMFive.search(tmp)
                    # mthreeMat=regMThree.search(tmp)
                    # if (mfiveMat and mthreeMat):
                    #     if (mfiveMat.span()[1]<mthreeMat.span()[0]):
                    #         tseq=tranStrand(tmp[mfiveMat.span()[1]:mthreeMat.span()[0]])
                    #         if(fastaTrigger):
                    #             outputFa.write(">"+header+"\n"+tseq+"\n")
                    #         if tseq in seqDict:
                    #             seqDict[tseq]+=1
                    #         else:
                    #             seqDict[tseq]=1
                    lab=0
                    header=""
                    if (count % 10000 == 0):
                        print("Checked "+ str(count)+ " reads ( "+ str(int(count/allReads*10000)/100)+"% )...\r",end="")
    print("The file, "+args.Fastq[f]+", analysis completed...")
                        
seqDict = sorted( seqDict.items(),key = lambda x:x[1],reverse = True)
for i in seqDict:
    outputStat.write(i[0]+"\t"+str(len(i[0]))+"\t"+str(i[1])+"\n")
    
    