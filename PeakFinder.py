import re
#WE WANT TO FIND PEAKS
#PEAKS ARE AREAS WHERE COVERAGE IN SEQUENCING DATA IS MUCH HIGHER, HIGHER READ COUNTS
#we need to approach this as a interval scheduling problem ,keeping track of number of reads per base pair
#we need a metric to define how many counts and how an increasing count trend can act as a peak.


#OUTLINE:

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

  
