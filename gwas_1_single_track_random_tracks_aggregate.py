# Purpose of file: To perform an analysis  to compare the performance of GWAS 1 single tracks to random 33-track sets for SNP selection.and
# Author: Sophie Sigfstead
# 
import random
def main():
    # Step 1. Define a dictionary that houses random track sets to run the single track analysis on. 
    # constants to set up the dictionary
    num_track_sets = 10
    num_tracks_per_set = 33
    track_set_dictionary = {}

    # brain tracks - to be excluded from the random track sets
    brain_tracks = {0, 1, 9, 76, 78, 80, 81, 172, 179, 216, 240, 261, 278,
            319, 326, 338, 355, 370, 403, 411, 421, 458, 462, 469,
            499, 524, 552, 580, 582, 602, 644, 669}
    
    # Create the keys - these will be the ids of the track sets
    keys = [chr(i) for i in range (ord('a'), ord('a')+num_track_sets)]

    # Create the list of tracks to chose from (non-brain)
    valid_numbers = [i for i in range(684) if i not in brain_tracks]

    for key in keys:
        track_set_dictionary[key] = random.sample(valid_numbers,num_tracks_per_set)

    print(track_set_dictionary)
    return 

if __name__ == "__main__":
    main()