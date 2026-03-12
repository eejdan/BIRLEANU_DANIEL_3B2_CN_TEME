import numpy as np
from scipy.linalg import lu_solve, lu_factor, lu

def getIntInput(prompt):
    return int(input(prompt))


def main():
    matrix_size = getIntInput("Size: ")
    
    testMatrix = matrixCreate(matrix_size)
    testCoeffVector = coeffVectorCreate(matrix_size)

    demoLUD(testMatrix, testCoeffVector)

    precisionFactor = getIntInput("Precision factor: ")
    
    demoCholesky(testMatrix, 10**(-precisionFactor))



def matrixCreate(size):
    seed = np.random.randint(0, 31, size=(size, size))
    return (seed @ seed.T).astype(float) 

def coeffVectorCreate(size):
    return np.random.randint(0, 31, size=size).astype(float)

def demoLUD(matrix, coeffVector):
    _, L, U = lu(matrix)
    print("L:\n", L)
    print("U:\n", U)
    print("Computing x from Ax = b");
    xlib = lu_solve(lu_factor(matrix), coeffVector)
    print("x:\n", xlib)

def computeCholesky(matrix, epsilon):
    # sa se calculeze o desc LDLt // choleski 
    # T = LD = (t_ij); U = Lt = (u_ij);
    # t_ij = d_j * l_ij; u_ij = l_ji
    D = np.zeros(matrix.shape[0])
    for p in range(matrix.shape[0]):
        D[p] = matrix[p, p]
        for k in range(p):
            D[p] -= D[k] * matrix[p, k] * matrix[p, k] # punem valorile l in matrix in loc de L
        
        if abs(D[p]) < epsilon:
            print(f"Impartire cu 0,p={p}");
            raise ValueError("Matricea trebuie pozitiv definita")

        for i in range(p+1, matrix.shape[0]):
            for k in range(p):
                matrix[i, p] -= D[k] * matrix[i, k] * matrix[p, k]
            matrix[i, p] /= D[p]

    return D, matrix

def demoCholesky(matrix, epsilon):
    D, L = computeCholesky(matrix, epsilon)
    print("D:\n", D)
    print("L:\n", L)


if __name__ == "__main__":
    main();
