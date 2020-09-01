import sys, argparse, re

class FreqSumParser:
    def __init__(self, filehandle, outputfilehandle):
        self.handle = filehandle
        self.output = outputfilehandle
        self.__readHeader()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        # while(True):
        line = next(self.handle)
        fields=line.strip().split()
        if fields[0][0:3] == "chr":
            Chrom = int(fields[0][3:])-1 #ignore "chr" if in start of chromosome name
        else:
            Chrom=int(fields[0])-1 #Convert Chromosome numbering to 0-based)
        Pos=int(fields[1])
        Ref=fields[2]
        Alt=fields[3]
        AlleleFreqs = [int(x) for x in fields[4:]]
        afDict = dict(zip(self.popNames, AlleleFreqs))
        # if self.noTransitions and ...isTransition...:
        #     continue
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
        print("#Available populations in Input File and their respective sizes: ", self.sizes, file=self.output)

def getJackknife(blockValues, totalObservations, blockSizes, skipJackknife):
    thetaminus=[0 for x in range(len(blockSizes))]
    sum1=0
    sum2=0
    jackknifeStdErr=0
    if sum(totalObservations)==0:
        ## If no observations were made, the rare allele sharing rate is 0.
        thetahat=0
    else:
        ## thetahat is normalised by number of observations
        thetahat=sum(blockValues)/sum(totalObservations)

    for c in range(len(blockSizes)):
        if totalObservations[c] == sum(totalObservations):
            ## If all rare alleles are on a single chunk, the ras/site without that chunk is 0.
            thetaminus[c]=0
        else:
            ## thetaminus is normalised by number of observations
            thetaminus[c]=( (sum(blockValues)-blockValues[c]) / (sum(totalObservations)-totalObservations[c]) )
        
        ## Jackknife estimator calculation
        sum1+=thetahat-thetaminus[c]
        sum2+=(blockSizes[c]*thetaminus[c])/sum(blockSizes)
    jackknifeEstimator=sum1+sum2
    
    ## Standard error calculation
    if not skipJackknife:
        for c in range(len(blockSizes)):
            hj=sum(blockSizes)/blockSizes[c]
            pseudoval = (hj*thetahat)-((hj-1)*thetaminus[c])
            jackknifeStdErr+=(1/len(blockSizes))*(((pseudoval-jackknifeEstimator)**2)/(hj-1))
    return (jackknifeEstimator,jackknifeStdErr)
