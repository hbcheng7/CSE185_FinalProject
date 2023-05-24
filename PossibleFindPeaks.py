import numpy as np
import scipy.stats as stats
import os

def findPeaks(tag_directory, input_tag_directory, output_file, peak_size=None, frag_length=None, min_dist=None, fdr=None, genome_size=2e9):
    # Load tags from each chromosome
    tags = load_tags(tag_directory)
    input_tags = load_tags(input_tag_directory)

    # Adjust tags to the center of their fragments or half of the estimated fragment length in the 3' direction
    adjusted_tags = adjust_tags(tags, frag_length)

    # Perform peak calling algorithm to identify clusters
    clusters = find_clusters(adjusted_tags, peak_size, min_dist)

    # Calculate the tag threshold based on the desired FDR
    tag_threshold = calculate_tag_threshold(tags, input_tags, fdr, genome_size)

    # Filter peaks based on the tag threshold
    peaks = filter_peaks(clusters, tag_threshold)

    # Write the peak calling results to the output file
    write_peaks_to_file(peaks, output_file)

def load_tags(tag_directory):
    tags = []
    chromosomes = get_chromosomes(tag_directory)
    for chromosome in chromosomes:
        tags_file = os.path.join(tag_directory, f"{chromosome}.tags.tsv")
        with open(tags_file, 'r') as file:
            for line in file:
                tag = int(line.strip())
                tags.append(tag)
    return tags

def get_chromosomes(tag_directory):
    chromosomes = []
    files = os.listdir(tag_directory)
    for file in files:
        if file.endswith(".tags.tsv"):
            chromosome = file[:-10]  # Remove the ".tags.tsv" extension
            chromosomes.append(chromosome)
    return chromosomes

def adjust_tags(tags, frag_length):
    adjusted_tags = []
    if frag_length is None:
        frag_length = estimate_frag_length(tags)
    for tag in tags:
        adjusted_tag = tag + (frag_length / 2)  # Adjust tag position by half of the fragment length
        adjusted_tags.append(adjusted_tag)
    return adjusted_tags

def estimate_frag_length(tags):
    diff_tags = np.diff(tags)
    frag_length = int(np.median(diff_tags))
    return frag_length

def find_clusters(tags, peak_size, min_dist):
    clusters = []
    sorted_tags = sorted(tags)
    start = sorted_tags[0]
    cluster_start = start
    cluster_end = start
    for i in range(1, len(sorted_tags)):
        if sorted_tags[i] - cluster_end <= min_dist:
            cluster_end = sorted_tags[i]
        else:
            cluster = (cluster_start, cluster_end)
            if cluster_end - cluster_start >= peak_size:
                clusters.append(cluster)
            cluster_start = sorted_tags[i]
            cluster_end = sorted_tags[i]
    cluster = (cluster_start, cluster_end)
    if cluster_end - cluster_start >= peak_size:
        clusters.append(cluster)
    return clusters

def calculate_tag_threshold(tags, input_tags, fdr, genome_size):
    expected_peaks = []
    for threshold in range(1, max(tags) + 1):
        observed_peaks = len([tag for tag in tags if tag >= threshold])
        expected_peaks.append(observed_peaks * genome_size / len(tags))

    false_positive_rates = []
    for threshold in range(1, max(tags) + 1):
        false_positives = len([tag for tag in input_tags if tag >= threshold])
        false_positive_rate = false_positives / genome_size
        false_positive_rates.append(false_positive_rate)

    for i in range(len(expected_peaks)):
        if false_positive_rates[i] <= fdr:
            tag_threshold = i + 1
            break

    return tag_threshold

def filter_peaks(clusters, tag_threshold):
    filtered_peaks = []
    for cluster in clusters:
        if cluster[1] - cluster[0] >= tag_threshold:
            filtered_peaks.append(cluster)
    return filtered_peaks

def write_peaks_to_file(peaks, output_file):
    with open(output_file, 'w') as file:
        for peak in peaks:
            file.write(f"{peak[0]}\t{peak[1]}\n")
