# I used ChatGPT to help me write me this complete code
def is_symmetric(matrix):
    """
    This function check if the matrix is symmetric or not.
    :param matrix: A square matrix represented as a list of lists.
    :return: True if the matrix is symmetric, False otherwise.
    """
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True


def is_positive_definite(matrix):
    n = len(matrix)
    try:
        for k in range(n):
            sub_matrix = [row[:k + 1] for row in matrix[:k + 1]]
            if minor(sub_matrix, k + 1) <= 0:
                return False
        return True
    except Exception as e:
        print(f"An error occurred: {e}")
        return False


def minor(matrix, n):
    """
    Compute the determinant of a matrix A using recursion and the cofactor expansion.
    :param A: Square matrix
    :return: Determinant of the matrix A
    """
    if n == 1:
        return matrix[0][0]
    det = 0
    for c in range(n):
        det += ((-1) ** c) * matrix[0][c] * minor([row[:c] + row[c + 1:] for row in matrix[1:]], n - 1)
    return det


def cholesky_decomposition(matrix):
    """
        Performs the Cholesky decomposition on a given square, symmetric, and positive definite matrix.
        Parameters:matrix (list of lists of floats): A square, symmetric, and positive definite matrix to be decomposed.
        Returns:lower (list of lists of floats): The lower triangular matrix resulting from the Cholesky decomposition.
        Note: The input matrix must be square (n x n), symmetric, and positive definite. This function does not perform checks for these properties and may produce incorrect results or errors if the input does not meet these criteria.The function uses only the lower triangular part of the input matrix, assuming it is symmetric.
    """
    n = len(matrix)
    lower = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            sum1 = sum(lower[i][k] * lower[j][k] for k in range(j))

            if i == j:  # Diagonal elements
                lower[i][j] = (matrix[i][i] - sum1) ** 0.5
            else:
                lower[i][j] = (matrix[i][j] - sum1) / lower[j][j]
    return lower


def doolittle_decomposition(matrix): # I did not do the Gauss_Seidel so I get help from ChatGPT to get this code
    n = len(matrix)
    L = [[0 if i != j else 1 for j in range(n)] for i in range(n)]
    U = [[0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for k in range(i, n):
            U[i][k] = matrix[i][k] - sum(L[i][j] * U[j][k] for j in range(i))
        for k in range(i, n):
            if i == k:
                L[i][i] = 1
            else:
                L[k][i] = (matrix[k][i] - sum(L[k][j] * U[j][i] for j in range(i))) / U[i][i]
    return L, U


def forward_substitution(L, b):
    n = len(b)
    y = [0 for _ in range(n)]
    for i in range(n):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]
    return y


def backward_substitution(U, y):
    n = len(y)
    x = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x


def solve_linear_system(matrix, b):
    if is_symmetric(matrix) and is_positive_definite(matrix):
        L = cholesky_decomposition(matrix)
        y = forward_substitution(L, b)
        x = backward_substitution([row[:] for row in map(list, zip(*L))], y)  # Transpose L to get U
    else:
        L, U = doolittle_decomposition(matrix)
        y = forward_substitution(L, b)
        x = backward_substitution(U, y)
    return x


# Problems
A = [
    [1, -1, 3, 2],
    [-1, 5, -5, -2],
    [3, -5, 19, 3],
    [2, -2, 3, 21]

]

b = [15, -35, 94, 1]
B = [
    [4, 2, 4, 0],
    [2, 2, 3, 2],
    [4, 3, 6, 3],
    [0, 2, 3, 9]

]

a = [20, 36, 60, 122]

solution1 = solve_linear_system(A, b)
solution2 = solve_linear_system(B, a)

print("Solution:", solution1)
print("Solution:", solution2)

