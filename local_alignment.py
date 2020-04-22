def local_alignment_matrix(seq1=str, seq2=str, scoring=dict):
    """ Returns an n+1 * m+1 matrix for the local alignment
    of seq1 and seq2 sequences. """
    m = len(seq1)
    n = len(seq2)
    matrix = [[0 for x in range(n+1)] for x in range(m+1)] 
    for i in range(m+1): 
        for j in range(n+1): 
            if i == 0 or j == 0: 
                matrix[i][j] = 0
            elif seq1[i-1] == seq2[j-1]:
                matrix[i][j] = matrix[i-1][j-1] + scoring["match"]
            else:
                s = max(matrix[i-1][j]-scoring["indel"], matrix[i][j-1]-scoring["indel"], 
                        matrix[i-1][j-1]-scoring["mismatch"], 0)
                if s == matrix[i-1][j]-scoring["indel"]:
                    matrix[i][j] = matrix[i-1][j] - scoring["indel"]
                elif s == matrix[i][j-1]-scoring["indel"]:
                    matrix[i][j] = matrix[i][j-1]-scoring["indel"]
                elif s == matrix[i-1][j-1]-scoring["mismatch"]:
                    matrix[i][j] = matrix[i-1][j-1] - scoring["mismatch"]
                else:
                    matrix[i][j] = 0
    return matrix
    

def local_alignment(seq1=str, seq2=str, scoring=dict):
    """ Performs local alignment of seq1 and seq2 sequencies. 
    Scoring holds the match, mismatch, and indel scores."""
    alignment1 = "" # from seq1
    alignment2 = "" # from seq2
    best_pos = -float("inf")
    matrix = local_alignment_matrix(seq1, seq2, scoring) # set up alignment matrix
    best_pos = ""
    best = 0
    for i in range(len(matrix)-1):
        for j in range(len(matrix[i])-1):
            if matrix[i][j] > best:
                best = matrix[i][j]
                best_pos = (i,j)
    i = best_pos[0]
    j = best_pos[1]
    while (matrix[i][j] != 0):
        if seq1[i-1] == seq2[j-1]:
            i -= 1
            j -= 1
            alignment1 += seq1[i]
            alignment2 += seq2[j]
        else:
            hori = matrix[i][j-1]
            diag = matrix[i-1][j-1]
            verti = matrix[i-1][j]
            step = max(hori, diag, verti)
            if step == hori:
                j -= 1
                alignment2 += seq2[j]
                alignment1 += "-"
            elif step == verti:
                i -= 1
                alignment1 += seq1[i]
                alignment2 += "-"
            else:
                i -= 1
                j -= 1
                alignment1 += seq1[i]
                alignment2 += seq2[j]
    return alignment1[::-1], alignment2[::-1]