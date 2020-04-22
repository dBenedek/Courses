def global_alignment_matrix(seq1=str, seq2=str, scoring=dict):
    """ Returns an n+1 * m+1 matrix for global alignment of seq1 and seq2 
    sequences. Scoring dict contains match, mismatch, and indel scores."""
    m = len(seq1)
    n = len(seq2)
    matrix = [[0 for x in range(n+1)] for x in range(m+1)] 
    for i in range(m+1): 
        for j in range(n+1): 
            if i == 0 and j == 0: 
                matrix[i][j] = 0
            elif i == 0 and j != 0: 
                matrix[i][j] = matrix[i][j-1] - scoring["indel"]
            elif i != 0 and j == 0: 
                matrix[i][j] = matrix[i-1][j] - scoring["indel"]
            elif seq1[i-1] == seq2[j-1]:
                matrix[i][j] = matrix[i-1][j-1] + scoring["match"]
            else:
                s = max(matrix[i-1][j]-scoring["indel"], matrix[i][j-1]-scoring["indel"], 
                        matrix[i-1][j-1]-scoring["mismatch"])
                if s == matrix[i-1][j]-scoring["indel"]:
                    matrix[i][j] = matrix[i-1][j] - scoring["indel"]
                elif s == matrix[i][j-1]-scoring["indel"]:
                    matrix[i][j] = matrix[i][j-1]-scoring["indel"]
                else:
                    matrix[i][j] = matrix[i-1][j-1] - scoring["mismatch"]
    return matrix


def global_alignment(seq1=str, seq2=str, scoring=dict):
    """ Returns the global alignment of seq1 and seq2 sequences."""
    matrix = global_alignment_matrix(seq1, seq2, scoring) # set up alignment matrix of seq1 and seq2
    i = len(seq1)
    j = len(seq2)
    alignment1 = "" # from seq1
    alignment2 = "" # from seq2
    while True: 
        if i == 0 and j == 0:
            break
        else:
            if seq1[i-1] == seq2[j-1] and i > 0 and j > 0: 
                i-=1
                j-=1
                alignment1 += seq1[i]
                alignment2 += seq2[j]
            elif i == 0 and j != 0:
                j -= 1
                alignment1 += "-"
                alignment2 += seq2[j]
            elif j == 0 and i != 0:
                i -= 1
                alignment1 += seq1[i]
                alignment2 += "-"
            else:
                diag = matrix[i-1][j-1] # mismatch
                hori = matrix[i][j-1] # indel
                verti = matrix[i-1][j] # indel
                best_step = max(diag, hori, verti)
                if diag == best_step: 
                    i -= 1
                    j -= 1
                    alignment1 += seq1[i]
                    alignment2 += seq2[j]
                elif hori == best_step:
                    j -= 1
                    alignment2 += seq2[j]
                    alignment1 += "-"
                else:
                    i -= 1
                    alignment1 += seq1[i]
                    alignment2 += "-"
    return alignment1[::-1], alignment2[::-1]
