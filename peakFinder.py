import os, math, sys, argparse, multiprocessing
from collections import deque
from scipy.stats import poisson


def processTagTSVs(dir_name):
    # Dictionary - key: chromsome ID, value: dictionary containing bp and reads at bp
    allCounts = {}
    
    # Checking if input file
    doFilt = True
    if dir_name == None:
        return None
    # Lists all files in given directory
    files = os.listdir(dir_name)
    
    
    for file_name in files:
        # Isolating only the tags files
        if file_name.endswith(".tags.tsv"):
            
            # Dictionary - key: base pair, value: reads associated with base pair
            fileCounts = {}
            
            with open(os.path.join(dir_name, file_name)) as tsv:
                for line in tsv:
                    # Turning line into list length 4 - [0] chr ID, [1] start bp of tag, [2] strand (0 forward, 1 reverse),
                    # [3] num reads for the tag, [4] length of tag
                    line = line.strip("\n")
                    line = line.split("\t")
                    line.pop(0)
                    
                    # Accounting for forward or reverse strand:
                    if line[2] == '1':
                        start = int(line[1]) - int(line[4])
                        end = int(line[1])
                    else:
                        start = int(line[1])
                        end = int(line[1]) + int(line[4])
                    
                    # Adding to dictionary, all bps within the tag
                    for i in range(start, end):
                        key = i
                        if key not in fileCounts.keys():
                            fileCounts[key] = int(float(line[3]))
                        else:
                            fileCounts.update({key : fileCounts[key]+int(float(line[3]))})
            if doFilt: fileCounts = {key : value for key, value in fileCounts.items() if value > 4}
            allCounts[file_name.split(".")[0]] = fileCounts
    
    sortKeys = list(allCounts.keys())
    sortKeys.sort()
    allCounts = {i : allCounts[i] for i in sortKeys}

    return allCounts

#processes taginfo stuff to return
def processTagInfo(file_path):
    peakSize = 0
    averageTagLength = 0
    tagsPerBP = 0
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("peakSizeEstimate="):
                p = line.strip().split("=")
                peakSize = float(p[1])
            if line.startswith("averageTagLength="):
                p = line.strip().split("=")
                averageTagLength = float(p[1])
            if line.startswith("tagsPerBP="):
                p = line.strip().split("=")
                tagsPerBP = float(p[1])
                    
    return peakSize, averageTagLength, tagsPerBP

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
    #find average height difference
    averageHeightDiff = (TFcounter - ControlCounter)/(windowEnd-windowStart)
    #if control or TF isn't 0, which it could be which would make value undefined
    if (ControlCounter != 0 and TFcounter != 0):
        return TFcounter/ControlCounter, maxPosition, averageHeightDiff
    return TFcounter/1, maxPosition, averageHeightDiff

def findPeaks(cr, cr_dict,input_cr_dict, fcThreshold, windowSize, estimatedPeakSize, averageTagLength, tagsPerBP):
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
        maxFoldChange, maxPosition, averageH = calcMaxFoldChange(windowStart, windowEnd, cr_dict, input_cr_dict)
        #check if this maxFoldChange is above our threshold
        if (maxFoldChange >= fcThreshold):
            count+=1
            #call pValue
            p_value = pValue(averageH, averageTagLength, tagsPerBP)
            #peak is the middle of the window +- half peak size on each end
            peak = (cr, int(maxPosition-(math.floor(estimatedPeakSize/2))),int(maxPosition+(math.floor(estimatedPeakSize/2))), p_value)
            possiblePeaks.append(peak)
    #print(possiblePeaks)        
    return possiblePeaks

def pValue(averageH, averageTagLength, tagsPerBP):
    #lambda value
    lambda_param = averageTagLength * tagsPerBP
    #observed value
    observed = averageH
    
    p_value = 1 - poisson.cdf(observed, lambda_param)
    
    #print ("Observed Value =", observed)
    #print ("Lambda Value =", lambda_param)
    #print("P-value =", p_value)
    return p_value

#function to process output BED file, write using python Output stream
# test -> paired_list = [(3,78),(99,174),(201,276)]
def makeBED(output_list):
    output_name = "myfile.bed"
    if (args.output is not None):
        output_name = args.output
    file = open(output_name, "w")
    with file as f:
            for output in output_list:
                f.write(str(output[0]) + "\t" + str(output[1]) + "\t" + str(output[2]) + "\t" + str(output[3]) + "\n")

def peakFinder(chr_tag, ctrl_tag, chr, peakSize, averageTagLength, tagsPerBP):
    if chr_tag == None or chr_tag == {}:
        return []
    windowSize = 2*float(averageTagLength)
    windowSizeRounded =  int(round(windowSize))
    return findPeaks(chr, chr_tag, ctrl_tag, 10, windowSizeRounded, peakSize, averageTagLength, tagsPerBP)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dir_name', type=str, help='input tags directory path')
    parser.add_argument('-c', '--control_tag_dir', type=str, default=None, help='input control tags directory path')
    parser.add_argument('-O', '--output', type=str, default=None, help='output file path')
    parser.add_argument('-t', '--threads', type=int, default=4, help='sets the amount of threads to use (default is 4)') 
    args = parser.parse_args()

    # Handling File not found exception, determining whether optional control tags given
    try:
        dir_files = os.listdir(args.dir_name)
        if args.control_tag_dir != None:
            input_tags = True
        else:
            input_tags = False
    except OSError as e:
        sys.exit("Invalid Directory Path ({})".format(e))


    tag_dir = args.dir_name
    ctrl_dir = args.control_tag_dir
    out_path = args.output
    peakSize, averageTagLength, tagsPerBP = processTagInfo(os.path.join(tag_dir, 'tagInfo.txt'))
    print("finished processing tagInfo.txt")
    
    pool = multiprocessing.Pool(processes=args.threads)
    directories = list(pool.map(processTagTSVs, [tag_dir, ctrl_dir]))
    print("finished processing tag directories")
    
    input = []
    for key, value in directories[0].items():
        try:
            input.append((value, directories[1][key], key, peakSize, averageTagLength, tagsPerBP))
        except TypeError as e:
            input.append((value, {}, key, peakSize, averageTagLength, tagsPerBP))
    results = list(pool.starmap(peakFinder, input))
    print("found peaks, printing to BED file")
    
    pool.close()
    possible_peaks = []
    
    for result in results:
        possible_peaks += result
    peaksSorted = sorted(possible_peaks, key=lambda x: x[3])
    makeBED(peaksSorted)
    #print(peaksSorted)

