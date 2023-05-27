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


#Function to read TSV Files from each chromosome, and select the chromosome of either user choice or by default select one with most data
def processTagTSVs(dir_name):
    allCounts = {}
    files = os.listdir(dir_name)
    for file_name in files:
        if file_name.endswith(".tags.tsv"):
            fileCounts = {}
            with open(os.path.join(dir_name, file_name)) as tsv:
                for line in tsv:
                    line = line.strip("\n")
                    line = line.split("\t")
                    for i in range(int(line[2]), int(line[2]) + int(line[-1])):
                        key = i
                        if key not in fileCounts.keys():
                            fileCounts[key] = int(float(line[4]))
                        else:
                            fileCounts.update({key : fileCounts[key]+int(float(line[4]))})
            fileCounts = {key : value for key, value in fileCounts.items() if value > 1}
            sortKeys = list(fileCounts.keys())
            sortKeys.sort()
            fileCounts = {i : fileCounts[i] for i in sortKeys}
            allCounts[file_name.split(".")[0]] = fileCounts
    
    sortKeys = list(allCounts.keys())
    sortKeys.sort()
    allCounts = {i : allCounts[i] for i in sortKeys}

    return allCounts
  
  
 return
#function to process trends in data and idefntiy peaks given thresholds of significance
#need a function to process foldChange
def calcMaxFoldChange(windowStart, windowEnd, cr_dict, input_cr_dict):
    maxVal = 0
    #in between startPos and startPos+windowSize, find the position with MAX value
    for i in range(windowStart, windowEnd):
        #check if it is the new max
        if i in cr_dict.keys():
            if cr_dict[i] > maxVal:
                maxVal = cr_dict[i]
        
    return maxVal

def findPeaks(cr_dict,input_cr_dict, fcThreshold, windowSize):
    possiblePeaks = []
    #hunter code outline below:
    #do it per chromosome
    #given a window size, that is how much you check. incrementn i by window size
    #temporarily work only on chr17
    #start window is first position
    windowStart = next(iter(cr_dict))
    windowEnd = windowStart + windowSize
    for pos in cr_dict:
        #if position is still in the window, since we already ran it, continue
        if pos < windowEnd:
            print ("COnt")
            continue
        else:
            #update window otherwise
            windowStart = pos
            windowEnd = pos + windowSize
        print ("start: " + str(windowStart))
        print ("end: " + str(windowEnd))
 
        print (pos)
        if (pos > 3000408):
            break
        #run maxFoldChange

        maxFoldChange = calcMaxFoldChange(windowStart, windowEnd, cr_dict, input_cr_dict)
        print ("MAX: " + str(maxFoldChange))
        #check if this maxFoldChange is above our threshold
        if (maxFoldChange >= fcThreshold):
            possiblePeaks.append((windowStart,windowEnd))

#function to process output BED file, write using python Output stream
# test -> paired_list = [(3,78),(99,174),(201,276)]
def makeBED(paired_list):
    file = open("myfile.bed", "w")
    with file as f:
            for pair in paired_list:
                f.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(1) + "\n")

makeBED(paired_list)

  
