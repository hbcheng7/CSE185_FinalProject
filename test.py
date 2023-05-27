import PeakFinder.py

def testProcessTagsTSVs():

    processTagTSVs("/tagdirs/test")
    i = 0
    print(len(allCounts))
    for chrom, cr_dict in allCounts.items():
        i = 0
        for x, y in cr_dict.items():
            print(x,y)
            i+=1
            if i > 100:
                break
testProcessTagsTSVs()

def testPeaks():
    cr_dict_sox2 = {}
    cr_dict_input = {}
    #get sox2 stuff for testing
    for chrom, cr_dict in allCounts.items():
        if chrom == "17":
            cr_dict_sox2 = cr_dict
    #get input stuff for testing
    for chrom, cr_dict in inputCounts.items():
        if chrom == "17":
            cr_dict_input = cr_dict
    findPeaks(cr_dict_sox2, cr_dict_input, fcThreshold = 2, windowSize = 5)
testPeaks()


paired_list = [(3,78),(99,174),(201,276)]
makeBED(paired_list)
