#WE WANT TO FIND PEAKS
#PEAKS ARE AREAS WHERE COVERAGE IN SEQUENCING DATA IS MUCH HIGHER, HIGHER READ COUNTS
#we need to approach this as a interval scheduling problem ,keeping track of number of reads per base pair
#we need a metric to define how many counts and how an increasing count trend can act as a peak.


#OUTLINE:

# IMPORTS
import numpy as np
import scipy.stats as stats
import os
import pwd
import re

# Define the tags, control tags, and output as command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('tags_tsv', type=str, help='input tags tsv file path')
parser.add_argument('control_tags', type=str, help='input control tags tsv file path')
parser.add_argument('-O', '--output', type=str, help='output file path')
args = parser.parse_args()

#We can call args.tags_tsv, args.control_tags, or args.output


#Function to read TSV Files from each chromosome, and then calling the peak finding on each
def processTagTSVs(dir_name):
    allCounts = {}
    files = os.listdir(dir_name)
    for file_name in files:
        if file_name.endswith(".tags.tsv"):
            processEachFile(file_name)
    return

#function that processes each tag.tsv file
#returns a dictionary of {position: counts} given each bp in the chromosome of choice
def processEachFile(file_name):
    fileCounts = {}
    with open(file_name) as tsv:
        for line in tsv:
            line = line.strip("\n")
            line = line.split("\t")
            #go from range position to length of read
            #start with forward strand ones
            if (int(line[3]) == 0):
                for i in range(int(line[2]), int(line[2]) + int(line[-1])):
                    key = i
                    #if position isn't there, create a key with it
                    if key not in fileCounts.keys():
                        fileCounts[key] = int(float(line[4]))
                    else:
                        #update the position already there
                        fileCounts.update({key : fileCounts[key]+int(float(line[4]))})
            #now do reverse strand
            if (int(line[3])==1):
                for i in range(int(line[2]), int(line[2]) - int(line[-1]),-1):
                    key = i
                    #if position isn't there, create a key with it
                    if key not in fileCounts.keys():
                        fileCounts[key] = int(float(line[4]))
                    else:
                        #update the position already there
                        fileCounts.update({key : fileCounts[key]+int(float(line[4]))})

                        

def calcMaxFoldChange(windowStart, windowEnd, cr_dict, input_cr_dict):
    #maxVal = 0
    maxPosition = (windowStart + windowEnd)/2
    TFcounter = 0
    ControlCounter = 0
    #in between startPos and startPos+windowSize, find the position with MAX value
    for i in range(windowStart, windowEnd):
        #check if it is the new max
        if i in cr_dict.keys():
            TFcounter += cr_dict[i]
        if i in input_cr_dict.keys():
            ControlCounter += input_cr_dict[i]
    #print ("TF: "+ str(TFcounter))
    #print ("CC: " + str(ControlCounter))
    #if control or TF isn't 0, which it could be which would make value undefined
    if (ControlCounter != 0 and TFcounter != 0):
        return TFcounter/ControlCounter, maxPosition
    return 0, 0
def findPeaks(cr_dict,input_cr_dict, fcThreshold, windowSize):
    possiblePeaks = []
    count = 0
    windowStart = next(iter(cr_dict))
    windowEnd = windowStart + windowSize
    #for key, value in input_cr_dict.items():  
      #  print ("key: " + str(key) + " value: " + str(value))
    #hunter code outline below:
    #do it per chromosome
    #given a window size, that is how much you check. incrementn i by window size
    #temporarily work only on chr17
    #start window is first position
    count = 0
    startofWindow = True
    for pos in cr_dict:
        #if position is still in the window, since we already ran it, continue
        if pos < windowEnd:
            continue
        else:
            #update window otherwise
            windowStart = pos
            windowEnd = pos + windowSize
        maxFoldChange, maxPosition = calcMaxFoldChange(windowStart, windowEnd, cr_dict, input_cr_dict)
        #print (maxFoldChange)
        #check if this maxFoldChange is above our threshold
        if (maxFoldChange >= fcThreshold):
            count+=1
            #peak is the middle of the window +- half peak size on each end
            peak = (maxPosition+(estimatedPeakSize/2),maxPosition-(estimatedPeakSize/2))
            possiblePeaks.append((windowStart,windowEnd, maxFoldChange))

    
    return possiblePeaks
       
    
#function to process output BED file, write using python Output stream
# test -> paired_list = [(3,78),(99,174),(201,276)]
def makeBED(paired_list):
    file = open("myfile.bed", "w")
    with file as f:
            for pair in paired_list:
                f.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(1) + "\n")

makeBED(paired_list)

  
