# utilities/diversity_index.py
from math import log as ln

# utilities/diversity_index.py
from math import log as ln

def p(n, N):
    if n == 0:
        return 0
    return (float(n) / N) * ln(float(n) / N)

def shannon_diversity(data):
    N = sum(data.values()) if isinstance(data, dict) else sum(data)
    if N == 0:
        return 0
    return -sum(p(n, N) for n in (data.values() if isinstance(data, dict) else data) if n > 0)

def evenness(H, S):
    if S <= 1:
        return 0
    return H / ln(S)

def simpsons_diversity_index(counts):
    """Calculates Simpson's Diversity Index for a given list of species counts.

    Args:
        counts: A list of integers representing the number of individuals
                of each species.

    Returns:
        The Simpson's Diversity Index.
    """

    N = sum(counts)
    sum_of_squares = sum(n * (n - 1) for n in counts)

    if N <= 1:
        return 0  # No diversity if only one individual or no individuals
    else:
        return 1 - (sum_of_squares / (N * (N - 1)))

    # Example usage
    #species_counts = [10, 5, 2, 8, 15]
    #diversity = simpsons_diversity_index(species_counts)
    #print("Simpson's Diversity Index:", diversity)

def inverse_simpson_index(counts):
    """
    Calculates the Inverse Simpson's Diversity Index for a given list of counts.

    Args:
        counts: A list of counts representing the abundance of each species.

    Returns:
        The Inverse Simpson's Diversity Index.
    """

    total_count = sum(counts)
    if total_count <= 0:
        return 0
    simpson_index = sum((count / total_count)**2 for count in counts)
    if simpson_index <= 0:
        return 0
    return 1 / simpson_index
