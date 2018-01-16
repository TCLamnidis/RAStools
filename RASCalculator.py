#!/usr/bin/env python3

import sys, argparse, re
from math import sqrt
from time import strftime

class FreqSumParser:
    def __init__(self, filehandle):
        self.handle = filehandle
        self.__readHeader()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        line = next(self.handle)
        fields=line.strip().split()
        Chrom=int(fields[0])-1 #Convert Chromosome numbering to 0-based)
        Pos=int(fields[1])
        Ref=fields[2]
        Alt=fields[3]
        AlleleFreqs = [int(x) for x in fields[4:]]
        afDict = dict(zip(self.popNames, AlleleFreqs))
        return (Chrom, Pos, Ref, Alt, afDict)
        
    def __readHeader(self):
        line = next(self.handle)
        fields=line.strip().split()
        self.sizes = {}
        self.popNames = []
        Pops=fields[4:]
        for p in Pops:
            splitPopName = re.split('[(|)]', p)
            popName = splitPopName[0]
            self.popNames.append(popName)
            popSize = int(splitPopName[1])
            self.sizes[popName] = popSize
        print("#Available populations in Input File and their respective sizes: ", self.sizes, file=args.Output)
    
########## MAIN ############

parser = argparse.ArgumentParser(description="Extract the frequency of shared rare variants between each left population and all right populations from a freqsum file. Also preforms error estimation using jackknifing, using the number of observed sites for normalisation.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input freqsum file. Omit to read from stdin.", required=False)
parser.add_argument("-O", "--Output", metavar="<OUTPUT FILE>", type=argparse.FileType('w'), help="The output file. Omit to print in stdout.")
parser.add_argument("-M", "--maxAF", metavar="<MAX ALLELE COUNT>", type=int, default=10, help="The maximum number of alleles (total) in the reference populations. The default maximum allele value is 10.", required=False)
parser.add_argument("-m", "--minAF", metavar="<MIN ALLELE COUNT>", type=int, default=2, help="The minimum number of alleles (total) in the reference populations. The default minimum allele count is 2.", required=False)
parser.add_argument("-L", "--LeftPops", type=str, metavar="POP1,POP2,...", required=True, help="Set the Test populations/individuals. RAS will be calculated between the Test and all Right populations.")
parser.add_argument("-R", "--RightPops", type=str, metavar="POP1,POP2,...", required=False, help="A list of comma-separated population names that should be considered when computing the allele frequency. Consider all populations if not provided.")
parser.add_argument("-NT", "--NoTransitions", action='store_true', help="When present, No Transitions are included in the output. Useful for ancient samples with damaged DNA.")
parser.add_argument("-P", "--Private", action='store_true', required=False, help="Restrict the RAS calculation to privately shared rare variants only.")
parser.add_argument("-C", "--NrChroms", type=int, metavar="<INT>", default=22, required=False, help="The number of chromosomes in the dataset. [22]")
parser.add_argument("-d", "--details", action='store_true', help="Print RAS calculations for each allele frequency, in addition to the total.")
args = parser.parse_args()

if args.minAF<1:
    parser.error("--minAF cannot be lower than 1.")
print ("Program began running at:", strftime("%D %H:%M:%S"), file=sys.stderr)
#If no input file given, read from stdin
if args.Input == None:
    args.Input = sys.stdin
#If no output is given, use stdout
if args.Output == None:
    args.Output = sys.stdout

NumBins=args.NrChroms
minAF=args.minAF
maxAF=args.maxAF
Transitions = {"A":"G", "G":"A","C":"T","T":"C"}
# RightIndex={} #Holds the INDICES of the Right pops
# LeftsIndex={} #Holds the INDICES of the Left pops

# #Read -S argument into sample list
# if args.LeftPops!=None:
#     Lefts.append(LeftPops)

freqSumParser = FreqSumParser(args.Input)
LeftPops = args.LeftPops.split(",") #Holds the NAMES of the Left pops
RightPops = args.RightPops.split(",") if args.RightPops != None else [n for n in freqSumParser.popNames if n not in LeftPops] #Holds the NAMES of the Right pops

for x in LeftPops:
    assert (x in freqSumParser.popNames), "{} not found in FreqSum".format(x)
for x in RightPops:
    assert (x in freqSumParser.popNames), "{} not found in FreqSum".format(x)

RAS = [[[[0 for i in range(NumBins)] for j in range(maxAF+1)] for k in range(len(RightPops))] for x in range(len(LeftPops))]
mj = [[[[0 for i in range(NumBins)] for j in range(maxAF+1)] for k in range(len(RightPops))] for x in range(len(LeftPops))]

for (Chrom, Pos, Ref, Alt, afDict) in freqSumParser:
    #Exclude transitions if the option is given.
    if args.NoTransitions:
        if Ref in Transitions and Alt == Transitions[Ref]:
            continue 
    AfSum=0
    #Calculate AfSum for each position
    for pop in afDict:
        if pop in RightPops:
            if afDict[pop] > 0:
                AfSum+=afDict[pop]
    #Only consider sites with ascertained AF between the provided ranges.
    if AfSum > maxAF or AfSum < minAF:
        continue
    for Lftidx, leftPop in enumerate(LeftPops):
        for Rgtidx, rightPop in enumerate(RightPops):
                        
            #Only consider Privately shared sites when the --Private option is provided.
            if args.Private:
                if AfSum != afDict[rightPop]:
                    continue
            #Handle missing data by reducing pop size for that position (take out two chromosomes). If left is missing, count site. no RAS.
            leftSize=freqSumParser.sizes[leftPop]
            rightSize=freqSumParser.sizes[rightPop]
             #For now we will assume that Rights cannot have missing data, since they are grouped. RAS calculation with some missingness in rights will come once we start grouping individuals in populations internally.
            # if afDict[rightPop] == -1:
            #     rightSize-=2
            
            #Case where right and left pops are different
            if afDict[leftPop] >= 0 and afDict[rightPop] >= 0 and leftPop != rightPop:
                RAS[Lftidx][Rgtidx][AfSum][Chrom]+=(afDict[leftPop]*afDict[rightPop])/(leftSize*rightSize)
                mj[Lftidx][Rgtidx][AfSum][Chrom]+=1
                RAS[Lftidx][Rgtidx][minAF-1][Chrom]+=(afDict[leftPop]*afDict[rightPop])/(leftSize*rightSize) #within "minAF-1" we store total RAS and observed sites, for Jackknife estimations on the totals.
                mj[Lftidx][Rgtidx][minAF-1][Chrom]+=1 #within "minAF-1" we store total RAS and observed sites, for Jackknife estimations on the totals.
                
            #Case where left pop is also right pop (within population)
            elif afDict[leftPop] >= 0 and afDict[rightPop] >= 0 and leftPop == rightPop:
                RAS[Lftidx][Rgtidx][AfSum][Chrom]+=(afDict[leftPop]*(afDict[leftPop]-1))/(leftSize*(leftSize-1))
                mj[Lftidx][Rgtidx][AfSum][Chrom]+=1
                RAS[Lftidx][Rgtidx][minAF-1][Chrom]+=(afDict[leftPop]*(afDict[leftPop]-1))/(leftSize*(leftSize-1)) #within "minAF-1" we store total RAS and observed sites, for Jackknife estimations on the totals.
                mj[Lftidx][Rgtidx][minAF-1][Chrom]+=1 #within "minAF-1" we store total RAS and observed sites, for Jackknife estimations on the totals.

#Jackknifing
Thetahat = [[[0 for j in range(maxAF+1)] for k in range(len(RightPops))] for x in range(len(LeftPops))]
Thetaminus=[[[[0 for c in range(NumBins)] for j in range(maxAF+1)] for k in range(len(RightPops))] for x in range(len(LeftPops))]
for i in range(minAF-1,maxAF+1): # M+1 to pick up all chromosomes (0-based to 1-based). 
    for x in range(len(LeftPops)):
        for j in range(len(RightPops)):
            if sum(mj[x][j][i])==0:
                Thetahat[x][j][i]=0
            else:
                Thetahat[x][j][i]=(sum(RAS[x][j][i])/sum(mj[x][j][i]))
            for c in range(NumBins):
                if mj[x][j][i][c] == sum(mj[x][j][i]):
                    Thetaminus[x][j][i][c]=0
                else:
                    # Thetaminus[j][i][c]=(sum(RAS[j][i]) - RAS[j][i][c]) / (sum(lengths) - lengths[c])
                    Thetaminus[x][j][i][c]=(sum(RAS[x][j][i]) - RAS[x][j][i][c]) / (sum(mj[x][j][i]) - mj[x][j][i][c])

ThetaJ=[[[0 for j in range(maxAF+1)] for k in range(len(RightPops))] for x in range(len(LeftPops))]
for i in range(minAF-1,maxAF+1): # M+1 to pick up all chromosomes (0-based to 1-based). mAF-1 to get the Thetas for the Sum of AFs too.
    for x in range(len(LeftPops)): 
        for j in range(len(RightPops)):
            Sum1=0
            Sum2=0
            for c in range(NumBins):
                Sum1+=Thetahat[x][j][i]-Thetaminus[x][j][i][c]
                if sum(mj[x][j][i]) >0:
                    Sum2+=((mj[x][j][i][c] * Thetaminus[x][j][i][c])/sum(mj[x][j][i]))
            ThetaJ[x][j][i]=Sum1+Sum2

Sigma2=[[[0 for j in range(maxAF+1)] for k in range(len(RightPops))] for x in range(len(LeftPops))]
for i in range(minAF-1,maxAF+1): # M+1 to pick up all chromosomes (0-based to 1-based). mAF-1 to get the Thetas for the Sum of AFs too.
    for x in range(len(LeftPops)):
        for j in range(len(RightPops)):
            for c in range(NumBins):
                if mj[x][j][i][c]==0:
                    pseudovalue=Thetahat[x][j][i]
                    Sigma2[x][j][i]+=0
                else:
                    hj=sum(mj[x][j][i])/mj[x][j][i][c]
                    pseudovalue=(hj*Thetahat[x][j][i])-((hj-1) * Thetaminus[x][j][i][c])
                    Sigma2[x][j][i]+=(((pseudovalue-ThetaJ[x][j][i])**2)/(hj-1))/NumBins

# print ("#FREQSUM POPULATIONS & SIZES:",*Pops, file=args.Output, sep=" ", end="\n")
print ("#Left Populations: ", *LeftPops, sep=" ", file=args.Output, end="\n")
print ("#Populations considered for allele frequency calculation (Rights):", *RightPops, file=args.Output, sep="\t", end="\n")
# RAS, number of sites, RAS /Site, stderr of (RAS/site), Allele Freq
print("RightPop","LeftPop","RAS","Number of sites","RAS/site", "Jackknife Error", "Allele Frequency", sep="\t", file=args.Output)
for leftidx, leftPop in enumerate(LeftPops):
    for rightidx, rightPop in enumerate(RightPops):
        if args.details:
            for m in range(minAF,maxAF+1):
                print (rightPop, leftPop, "{:.5}".format(float(sum(RAS[leftidx][rightidx][m]))), "{:.15e}".format(sum(mj[leftidx][rightidx][m])), "{:.15e}".format(Thetahat[leftidx][rightidx][m]), "{:.15e}".format(sqrt(Sigma2[leftidx][rightidx][m])),m, sep="\t", file=args.Output)
        m=minAF-1
        print (rightPop, leftPop, "{:.5}".format(float(sum(RAS[leftidx][rightidx][m]))), "{:.15e}".format(sum(mj[leftidx][rightidx][m])), "{:.15e}".format(Thetahat[leftidx][rightidx][m]), "{:.15e}".format(sqrt(Sigma2[leftidx][rightidx][m])),"Total [{},{}]".format(minAF,maxAF), sep="\t", file=args.Output)
        #print ("", file=args.Output)

print ("Program finished running at:", strftime("%D %H:%M:%S"), file=sys.stderr)

