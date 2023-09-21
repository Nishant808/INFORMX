def align_sequences(sequence1, sequence2):
    # Initialize the scoring matrix and traceback matrix
    gap_penalty = -1
    match_score = 1
    mismatch_penalty = -1

    rows = len(sequence1) + 1
    cols = len(sequence2) + 1

    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    traceback_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initialize the first row and column of the scoring matrix
    for i in range(rows):
        score_matrix[i][0] = i * gap_penalty
    for j in range(cols):
        score_matrix[0][j] = j * gap_penalty

    # Fill in the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            if sequence1[i - 1] == sequence2[j - 1]:
                match = score_matrix[i - 1][j - 1] + match_score
            else:
                match = score_matrix[i - 1][j - 1] + mismatch_penalty

            gap1 = score_matrix[i - 1][j] + gap_penalty
            gap2 = score_matrix[i][j - 1] + gap_penalty

            # Choose the maximum score from the three possibilities
            score_matrix[i][j] = max(match, gap1, gap2)

            # Set the traceback direction (1 for diagonal, 2 for up, 3 for left)
            if score_matrix[i][j] == match:
                traceback_matrix[i][j] = 1
            elif score_matrix[i][j] == gap1:
                traceback_matrix[i][j] = 2
            else:
                traceback_matrix[i][j] = 3

    # Traceback to find the aligned sequences and compute the alignment score
    aligned_sequence1 = ""
    aligned_sequence2 = ""
    i, j = rows - 1, cols - 1

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 1:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 2:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = "-" + aligned_sequence2
            i -= 1
        else:
            aligned_sequence1 = "-" + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            j -= 1

    alignment_score = score_matrix[rows - 1][cols - 1]

    return aligned_sequence1, aligned_sequence2, alignment_score

if __name__ == "__main__":
    print("Welcome to Sequence Alignment Program")
    sequence1 = input("Enter the first sequence: ")
    sequence2 = input("Enter the second sequence: ")

    aligned_seq1, aligned_seq2, score = align_sequences(sequence1, sequence2)

    print("Aligned Sequence 1:", aligned_seq1)
    print("Aligned Sequence 2:", aligned_seq2)
    print("Alignment Score:", score)
