from Bio import Entrez, SeqIO
import numpy as np


def affine_gap_penalty_global_alignment(
    seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty
):
    m, n = len(seq1), len(seq2)

    # Initialize matrices
    M = np.zeros((m + 1, n + 1))  # Match/mismatch matrix
    X = np.zeros((m + 1, n + 1))  # Gap in sequence 1 matrix
    Y = np.zeros((m + 1, n + 1))  # Gap in sequence 2 matrix

    # Initialize first row and column with gap penalties
    for i in range(1, m + 1):
        X[i][0] = -gap_open_penalty - (i - 1) * gap_extension_penalty
    for j in range(1, n + 1):
        Y[0][j] = -gap_open_penalty - (j - 1) * gap_extension_penalty

    # Fill in matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = M[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else -mismatch_penalty
            )
            gap_x = (
                X[i - 1][j] - gap_extension_penalty
                if X[i - 1][j] > Y[i - 1][j]
                else Y[i - 1][j] - gap_open_penalty - gap_extension_penalty
            )
            gap_y = (
                Y[i][j - 1] - gap_extension_penalty
                if Y[i][j - 1] > X[i][j - 1]
                else X[i][j - 1] - gap_open_penalty - gap_extension_penalty
            )

            M[i][j] = max(match, gap_x, gap_y)
            X[i][j] = max(
                M[i - 1][j] - gap_open_penalty - gap_extension_penalty,
                X[i - 1][j] - gap_extension_penalty,
            )
            Y[i][j] = max(
                M[i][j - 1] - gap_open_penalty - gap_extension_penalty,
                Y[i][j - 1] - gap_extension_penalty,
            )

    # Traceback to find the alignment
    alignment_seq1 = ""
    alignment_seq2 = ""
    i, j = m, n
    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and M[i][j]
            == M[i - 1][j - 1]
            + (match_score if seq1[i - 1] == seq2[j - 1] else -mismatch_penalty)
        ):
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            i -= 1
            j -= 1
        elif i > 0 and M[i][j] == X[i][j]:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = "-" + alignment_seq2
            i -= 1
        elif j > 0 and M[i][j] == Y[i][j]:
            alignment_seq1 = "-" + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            j -= 1
        else:
            break  # Agregamos una condición de salida para evitar el bucle infinito

    total_score = M[m][n]
    return alignment_seq1, alignment_seq2, total_score


# Número de acceso
accession_number1 = input("\nNumero de Secuencia 1: ")
accession_number2 = input("\nNumero de Secuencia 2: ")
print("\n")

# Búsqueda en GenBank
Entrez.email = (
    "your-email@example.com"  # Reemplaza con tu dirección de correo electrónico
)
handle = Entrez.efetch(
    db="nucleotide", id=accession_number1, rettype="gb", retmode="text"
)
record = SeqIO.read(handle, "genbank")
seq1 = record.seq

handle = Entrez.efetch(
    db="nucleotide", id=accession_number2, rettype="gb", retmode="text"
)
record = SeqIO.read(handle, "genbank")
seq2 = record.seq

# Parámetros de la alineación ingresados por teclado
match_score = int(input("Ingrese el costo de una coincidencia: "))
mismatch_penalty = int(input("Ingrese el costo de una falta de coincidencia: "))
gap_open_penalty = int(input("Ingrese el costo de apertura de gap: "))
gap_extension_penalty = int(input("Ingrese el costo de extensión de gap: "))

# Realizar la alineación
alignment_seq1, alignment_seq2, total_score = affine_gap_penalty_global_alignment(
    seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty
)

# Imprimir resultados
print("Secuencia 1 alineada:", alignment_seq1)
print("Secuencia 2 alineada:", alignment_seq2)
print("Puntaje total de la alineación:", total_score)
