def profile_matrix(sequences=list):
    """ Returns the profile matrix of given sequences.
        Extended with Laplace’s Rule of Succession. """
    profile = {} 
    characters = set() # characters in sequences (eg. A, G, T, C)
    for seq in sequences:
        for char in seq:
            characters.add(char)
    for item in characters:
        profile[item] = [1 for i in range(len(sequences[0]))]
    for seq in sequences:
        for ind,nucl in enumerate(seq):
            profile[nucl][ind] += 1
    normalized_profile = profile.copy()
    for key, value in profile.items():
        for idx, score in enumerate(value):
            normalized = score/(len(sequences)+4)
            normalized_profile[key][idx] =  round(normalized, 2)
    return normalized_profile


def kmer_mostprobable(dna=str, k=int, profile=dict):
    """ Returns the most probable k-mer in text, which is the most similar
        to 'profile' profile matrix."""
    max_score = -float('inf')
    most_probable_kmer = ''
    for idx in range(len(dna) - k + 1):
        motif = dna[idx:idx + k]
        motif_score = 1
        for idx2, nuc in enumerate(motif):
            motif_score *= profile[nuc][idx2]
        if motif_score > max_score:
            max_score = motif_score
            most_probable_kmer = motif
    return most_probable_kmer


def consensus(sequences=list):
    """ Returns consensus sequence created from sequences in list."""
    consensus = ""
    for i in range(len(sequences[0])):
        column_values = ""
        for seq in sequences:
            column_values += seq[i]
        freq = 0
        most_freq_letter = ""
        for nucl in set(column_values):
            if column_values.count(nucl) > freq:
                freq = column_values.count(nucl)
                most_freq_letter = nucl 
        consensus += most_freq_letter
    return consensus
