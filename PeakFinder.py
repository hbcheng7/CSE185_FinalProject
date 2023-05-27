import re
#WE WANT TO FIND PEAKS
#PEAKS ARE AREAS WHERE COVERAGE IN SEQUENCING DATA IS MUCH HIGHER, HIGHER READ COUNTS
#we need to approach this as a interval scheduling problem ,keeping track of number of reads per base pair
#we need a metric to define how many counts and how an increasing count trend can act as a peak.


#OUTLINE:

#Function to read TSV Files from each chromosome, and select the chromosome of either user choice or by default select one with most data
def processTagTSVs():
  #process a dictionary list of start end intervals, of pair each. Key = start coordinate, value = (START, END)
  intervalList = {0:(0,1)}
  #for loop given amount of intervals to add to the list
  return intervalList

#Function to count number of reads per BP. 
def countReads(intervalList):
  #keep track of each BP and its read count in a dictionary. Key = int bp, Value = Count of Reads
  allCounts = {}
  #go through each interval and count up number of reads per BP. time compleixty might suck but start with this
  
  
 return
#function to process trends in data and idefntiy peaks given thresholds of significance
def findPeaks():
  #return a list of start end values of peaks
  return peaks

#function to process output BED file, write using python Output stream
# test -> paired_list = [(3,78),(99,174),(201,276)]
def makeBED(paired_list):
    file = open("myfile.bed", "w")
    with file as f:
            for pair in paired_list:
                f.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(1) + "\n")

makeBED(paired_list)

  
