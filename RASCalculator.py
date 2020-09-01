#!/usr/bin/env python3

import sys, argparse, re
from math import sqrt
from time import strftime
import RASUtils as ras

########## MAIN ############

parser = argparse.ArgumentParser(description="Compute rare allele sharing statistics between two populations with respect to an outgroup, as well as outgroup F3 statistics. Also preforms error estimation using jackknifing, using the number of observed sites for normalisation.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input freqsum file. Omit to read from stdin.", required=False)
parser.add_argument("-O", "--Output", metavar="<OUTPUT FILE>", type=argparse.FileType('w'), help="The output file. Omit to print in stdout.")
parser.add_argument("-o", "--outgroup", metavar="POP", type=str, help="The outgroup to polarize all alleles with. As a useful standard choice, use Chimp as outgroup", required=True)
parser.add_argument("-M", "--maxAF", metavar="<MAX ALLELE COUNT>", type=int, default=10, help="The maximum number of alleles (total) in the reference populations. The default maximum allele value is 10.", required=False)
parser.add_argument("-m", "--minAF", metavar="<MIN ALLELE COUNT>", type=int, default=2, help="The minimum number of alleles (total) in the reference populations. The default minimum allele count is 2.", required=False)

parser.add_argument("-a", "--ascertainIn", type=str, metavar="POP1,POP2,...", required=False, help="The populations to ascertain in (i.e. Rights). Rare alleles will be defined as alleles with minor allele freuency between <minAF> and <maxAF> in this set of populations. By default, all non-Left populations will be used for ascertainment.")
parser.add_argument("-c", "--calculateFor", type=str, metavar="POP1,POP2,...", required=True, help="The populations to be used in RAS calculations. RAS will be calculated between each Left population and each of these populations. By default, RAS will be calculated for all non-Left populations.")


parser.add_argument("-L", "--LeftPops", type=str, metavar="POP1,POP2,...", required=True, help="Set the Left populations/individuals. RAS will be calculated between each Left and all populations in the calculation list.")
# parser.add_argument("-R", "--RightPops", type=str, metavar="POP1,POP2,...", required=False, help="A list of comma-separated population names that should be considered when computing the allele frequency. Consider all populations if not provided.")

parser.add_argument("-x", "--MissingnessCutoff", type=float, metavar="<CUTOFF>", default=0.0, help="Missingness cutoff proportion for Right populations. E.g. 0.1: If more than 10%% of individuals in Right populations show missing data, the variant will be ignored. [default=0]")
parser.add_argument("-NT", "--NoTransitions", action='store_true', help="When present, No Transitions are included in the output. Useful for ancient samples with damaged DNA.")
parser.add_argument("-P", "--Private", action='store_true', required=False, help="Restrict the RAS calculation to privately shared rare variants only.")
parser.add_argument("-C", "--NrChroms", type=int, metavar="<INT>", default=22, required=False, help="The number of chromosomes in the dataset. [22]")
parser.add_argument("-d", "--details", action='store_true', help="Print RAS calculations for each allele frequency, in addition to the total.")
parser.add_argument("--f3FreqCutoff", type=float, metavar="<FREQ>", default=0.05, help="minimum minor allele frequency for f3 values. Defaults to 0.05. This is recommended to skip rare mutations and reduce the sensitivity of the F3 statistics on sequencing errors and DNA damage.")
parser.add_argument("--skipJackknife", action='store_true', default=False, help="When provided, the jackknife error is not estimated, but instead set to 0. [False]")
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
skipJK=args.skipJackknife

# RightIndex={} #Holds the INDICES of the Right pops
# LeftsIndex={} #Holds the INDICES of the Left pops

# #Read -S argument into sample list
# if args.LeftPops!=None:
#     Lefts.append(LeftPops)

freqSumParser = ras.FreqSumParser(args.Input, args.Output)
LeftPops = args.LeftPops.split(",") #Holds the NAMES of the Left pops
RightPops = args.ascertainIn.split(",") if args.ascertainIn != None else [n for n in freqSumParser.popNames if n not in LeftPops] #Holds the NAMES of the Right pops
TestPops = args.calculateFor.split(",") if args.calculateFor != None else [n for n in freqSumParser.popNames if n not in LeftPops]

# if Focal != None:
#     assert (Focal in freqSumParser.popNames), "Focal population '{}' not found in FreqSum".format(Focal)
for x in LeftPops:
    assert (x in freqSumParser.popNames), "Population {} not found in FreqSum".format(x)
for x in RightPops:
    assert (x in freqSumParser.popNames), "Population {} not found in FreqSum".format(x)
for x in TestPops:
    assert (x in freqSumParser.popNames), "Population {} not found in FreqSum".format(x)
assert (args.outgroup in freqSumParser.popNames), "Population {} not found in FreqSum".format(args.outgroup)

def getMissingness(afDict):
    missing = 0
    total = 0
    for x in RightPops:
        n = freqSumParser.sizes[x]
        if afDict[x] == -1:
            missing += n
        total += n
    return float(missing) / float(total)

def isTransition(ref, alt):
    Transitions = {"A":"G", "G":"A", "C":"T", "T":"C"}
    return (ref in Transitions and alt == Transitions[ref])

def getTotalDerivedAC(afDict):
    #Calculate AfSum for each position
    nonRefAC = 0
    totalCount = 0
    for pop in RightPops:
        if afDict[pop] >= 0:
            nonRefAC += afDict[pop]
            totalCount += freqSumParser.sizes[pop]
    xOutgroup = afDict[args.outgroup] / freqSumParser.sizes[args.outgroup]
    if xOutgroup < 0.5:
        return nonRefAC
    else:
        return totalCount - nonRefAC

# Bin minAF - 1: Total rare allele sharing
# Bins minAF -> maxAF: Rare Allele sharing per allele count
# Bin maxAF + 1: Outgroup F3 stats
RAS = [[[[0 for i in range(NumBins)] for j in range(maxAF+2)] for k in range(len(TestPops))] for x in range(len(LeftPops))]
# The normalization only records a total for all allele frequencies.
blockSizes =  [0 for i in range(NumBins)] ## Bin sizes are now stable across all RAS calculations.
mj =  [[ [0 for i in range(NumBins)]                          for k in range(len(TestPops))] for x in range(len(LeftPops))]


totalRightPopSize = sum(freqSumParser.sizes[p] for p in RightPops)
lineCount = 0
for (Chrom, Pos, Ref, Alt, afDict) in freqSumParser:

    #Skip sites where missingness in rightpops is above specified cutoff.
    lineCount += 1
    if lineCount % 10000 == 0:
        print("processing position {}:{}".format(Chrom + 1, Pos), file=sys.stderr)
    
    #Exclude transitions if the option is given.
    if args.NoTransitions and isTransition(Ref, Alt):
        continue     
    
    blockSizes[Chrom] += 1
    
    missingness = getMissingness(afDict)
    if missingness > args.MissingnessCutoff:
        continue
    
    if afDict[args.outgroup] < 0:
        continue
    
    
    derivedAC = getTotalDerivedAC(afDict)
    derivedAlleleFreq = derivedAC / totalRightPopSize
    minorAlleleFreq = derivedAlleleFreq if derivedAlleleFreq < 0.5 else 1.0 - derivedAlleleFreq

    for Lftidx, leftPop in enumerate(LeftPops):
        for Tstidx, testPop in enumerate(TestPops):
                        
            #Only consider Privately shared sites when the --Private option is provided.
            isPrivate = (derivedAC == afDict[testPop])
            leftSize = freqSumParser.sizes[leftPop]
            rightSize = freqSumParser.sizes[testPop]
            
            if afDict[leftPop] >= 0 and afDict[testPop] >= 0:
                mj[Lftidx][Tstidx][Chrom] += 1
                xLeft = afDict[leftPop] / leftSize
                xRight = afDict[testPop] / rightSize
                xOutgroup = afDict[args.outgroup] / freqSumParser.sizes[args.outgroup]
                add = (xLeft - xOutgroup) * (xRight - xOutgroup)
                if minorAlleleFreq >= args.f3FreqCutoff:
                    RAS[Lftidx][Tstidx][maxAF+1][Chrom] += add # For Outgroup F3 Stats
                if derivedAC >= minAF and derivedAC <= maxAF and (not args.Private or isPrivate):
                    
                    # For rare allele sharing, we are rounding the outgroup to 1 or 0, since we noticed that sequencing errors and spurious DNA damage can cause negative ras statistics when the same formula is used as for F3 stats. To see this, consider the following case: xOutgroup > 0, xLeft = 1 (sequencing error or DNA damage), rRight = 0. Then add < 0 and in some cases the entire statistics will be zero. For F3 stats, this effect will average out, but for low frequency RAS the total effect can systematically shift the statistics towards negative numbers.
                    
                    add = xLeft * xRight if xOutgroup < 0.5 else (1.0 - xLeft) * (1.0 - xRight)
                    #Only consider sites with ascertained minor AF between the provided ranges.
                    RAS[Lftidx][Tstidx][derivedAC][Chrom] += add
                    #within "minAF-1" we store total Rare allele sharing.
                    RAS[Lftidx][Tstidx][minAF-1][Chrom] += add

#Jackknifing
ThetaJ=[[[0 for j in range(maxAF+2)] for k in range(len(TestPops))] for x in range(len(LeftPops))]
Sigma2=[[[0 for j in range(maxAF+2)] for k in range(len(TestPops))] for x in range(len(LeftPops))]
for x in range(len(LeftPops)):
    for j in range(len(TestPops)):
        for i in range(minAF-1,maxAF+2):
            thetaJ,sigma2 = ras.getJackknife(RAS[x][j][i], mj[x][j], blockSizes, skipJK)
            ThetaJ[x][j][i] = thetaJ
            Sigma2[x][j][i] = sigma2

# print ("#FREQSUM POPULATIONS & SIZES:",*Pops, file=args.Output, sep=" ", end="\n")
print ("#Left Populations: ", *LeftPops, sep=" ", file=args.Output, end="\n")
print ("#Tested Populations: ", *TestPops, sep=" ", file=args.Output, end="\n")
print ("#Populations considered for allele frequency calculation (Rights):", *RightPops, file=args.Output, sep="\t", end="\n")
print ("#Outgroup: ", args.outgroup, file=args.Output, sep="\t", end="\n")
# RAS, number of sites, RAS /Site, stderr of (RAS/site), Allele Freq
print("TestPop","LeftPop","RAS","NumSites","RAS/site", "Error", "AlleleCount", sep="\t", file=args.Output)
for leftidx, leftPop in enumerate(LeftPops):
    for tstidx, testPop in enumerate(TestPops):
        if args.details:
            for m in range(minAF,maxAF+1):
                print (testPop, leftPop, "{:.5}".format(float(sum(RAS[leftidx][tstidx][m]))), "{:.15e}".format(sum(mj[leftidx][tstidx])), "{:.15e}".format(ThetaJ[leftidx][tstidx][m]), "{:.15e}".format(sqrt(Sigma2[leftidx][tstidx][m])),m, sep="\t", file=args.Output)
        m=minAF-1
        print (testPop, leftPop, "{:.5}".format(float(sum(RAS[leftidx][tstidx][m]))), "{:.15e}".format(sum(mj[leftidx][tstidx])), "{:.15e}".format(ThetaJ[leftidx][tstidx][m]), "{:.15e}".format(sqrt(Sigma2[leftidx][tstidx][m])),"Total [{},{}]".format(minAF,maxAF), sep="\t", file=args.Output)
        m=maxAF+1
        print (testPop, leftPop, "{:.5}".format(float(sum(RAS[leftidx][tstidx][m]))), "{:.15e}".format(sum(mj[leftidx][tstidx])), "{:.15e}".format(ThetaJ[leftidx][tstidx][m]), "{:.15e}".format(sqrt(Sigma2[leftidx][tstidx][m])),"Outgroup F3", sep="\t", file=args.Output)
        #print ("", file=args.Output)

print ("Program finished running at:", strftime("%D %H:%M:%S"), file=sys.stderr)

