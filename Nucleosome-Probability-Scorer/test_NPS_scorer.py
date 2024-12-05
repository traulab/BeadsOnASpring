import numpy as np
from scipy import stats

def precompute_distributions(nwps_frag_range, mode_DNA_length):
    """
    Precompute fragment score distributions for each fragment length in the range
    and store them in a dictionary.
    """
    distributions = {}
    
    for fragment_length in nwps_frag_range:
        # Create the distribution for the current fragment length
        if fragment_length < mode_DNA_length:
            total_length = mode_DNA_length + (mode_DNA_length - fragment_length)
        else:
            total_length = fragment_length

        # Calculate dynamic midpoint based on mode_DNA_length
        midpoint = (mode_DNA_length - 1) // 2
        second_half_start = midpoint + 1

        scores = np.zeros(total_length)

        for i in range(total_length):
            if i <= midpoint:
                scores[i] = i / midpoint  # Scale up the score to peak at 1
            elif i <= mode_DNA_length - 1:
                scores[i] = 1 - (i - second_half_start) / midpoint  # Adjust to descend back to 0

        # Second distribution: flip the first distribution
        end_scores = scores[::-1]  # Reverse the whole scores array

        # Sum the two distributions
        combined_scores = scores + end_scores

        central_score = np.mean(combined_scores)
        # print(central_score)

        centered_scores = [x - central_score for x in combined_scores]

        # Store the centered scores in the dictionary
        distributions[fragment_length] = centered_scores

    return distributions

length = 206

distributions = precompute_distributions([length], 166)

for value in distributions[length]:
    print(value)