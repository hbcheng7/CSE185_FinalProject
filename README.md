# CSE185_FinalProject
We are making a tool to find peaks given transcription factor data. We will take an input of a tag directory of a certain input genome and a control genome for comparison, and then utilize an approach to estimate where the peaks are.   
Our output will be a .BED file written to the directory of the user's choice.  
Tag directories MUST be generated prior to running our program, with tag files in the TSV format.  
Tag directories MUST share file names for respective chromosomes and end with .tags.tsv (i.e tag file for chromosome 17 must be 17.tags.tsv in both directories). 

# Calling the Function:
python PeakFinder.py <tag_directory_path> <control_tag_directory_path> -O <optional_output_path> 
