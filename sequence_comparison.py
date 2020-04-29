import random


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


def score(sequences=list):
    """ Sum of Hamming-distance scores at each nucleotide position."""
    final_score = 0
    for i in range(len(sequences[0])):
        column_values = ""
        for seq in sequences:
            column_values += seq[i] 
        freq = 0
        for nucl in set(column_values):
            if column_values.count(nucl) > freq:
                freq = column_values.count(nucl)   
        score = len(column_values) - freq
        final_score += score
    return final_score


def random_motif_search(dna=list, k=int):
    """ Returns kmers in dna list, which have the lowest score.
    (So they are more identical.)"""
    random_pos_start = random.randint(0, (len(dna[0])-k))
    random_pos_end = random_pos_start + k
    random_kmers = [seq[random_pos_start:random_pos_end] for seq in dna]
    while True:
        prof_matrix = profile_matrix(random_kmers)
        prob_kmers = []
        for seq in dna:
            most_prob_kmer = kmer_mostprobable(seq, k, prof_matrix)
            prob_kmers.append(most_prob_kmer)
        score1 = score(random_kmers)
        score2 = score(prob_kmers)
        if score1 > score2:
            random_kmers = prob_kmers[:]
        return random_kmers
    
    
def gibbs_sampling(dna=list, k=int, t=int, N=int): 
    """ Returns kmers in 'dna' list, which have the lowest score.
    (So they are the most identical.)"""
    random_pos_start = random.randint(0, (len(dna[0])-k))
    random_pos_end = random_pos_start + k
    random_kmers = [seq[random_pos_start:random_pos_end] for seq in dna] 
    for i in range(0,N):
        seq_out = random.choice(random_kmers)
        seq_out_ind = 0
        for i, s in enumerate(random_kmers):
            if s == seq_out:
                seq_out_ind = i
        cor_random_kmers = [i for i in random_kmers if i != seq_out]
        prof_matrix = profile_matrix(cor_random_kmers)
        i_kmer = kmer_mostprobable(seq_out, k, prof_matrix)
        new_motif = cor_random_kmers[:]
        new_motif[seq_out_ind] = i_kmer
        if score(random_kmers) > score(new_motif):
            random_kmers = new_motif
        return random_kmers
