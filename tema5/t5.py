import numpy as np
from scipy.linalg import cholesky

programEpsilon = 1e-8

def getIntInput(prompt):
	return int(input(prompt))

def generateSymmetricMatrix(size):
    AMatrix = np.random.rand(size, size).astype(np.float64) * 31.0
    
    for i in range(size):
        for j in range(i + 1, size):
            AMatrix[j, i] = AMatrix[i, j]
    
    return AMatrix


def checkSymmetry(AMatrix):
    return np.allclose(AMatrix, AMatrix.T, atol=programEpsilon)

def rotationMatrix(size, p, q, c, s):
    RMatrix = np.identity(size, dtype=np.float64)
    RMatrix[p, p] = c
    RMatrix[q, q] = c
    RMatrix[p, q] = s
    RMatrix[q, p] = -s
    return RMatrix

def computePQ(AMatrix):
    n = AMatrix.shape[0]
    p, q = 1, 0
    max_val = abs(AMatrix[p, q])
    for i in range(n):
        for j in range(i):
            # print("i=", i, "j=", j, "A[i,j]=", AMatrix[i, j], "max_val=", max_val)
            if abs(AMatrix[i, j]) > max_val:
                max_val = abs(AMatrix[i, j])
                p, q = i, j
    return p, q


def computeCST(AMatrix, p, q):
    if abs(AMatrix[p, q]) < programEpsilon:
        return 1.0, 0.0, 0.0

    alpha = (AMatrix[p, p] - AMatrix[q, q]) / (2.0 * AMatrix[p, q])

    if abs(alpha) < programEpsilon:
        t = 1.0
    else:
        t = np.sign(alpha) / (abs(alpha) + np.sqrt(alpha * alpha + 1.0))

    c = 1.0 / np.sqrt(1.0 + t * t)
    s = t * c

    return c, s, t

# def checkDiagonality(AMatrix):
#     n = AMatrix.shape[0]
#     for i in range(n):
#         for j in range(i):
#             if abs(AMatrix[i, j]) > programEpsilon:
#                 return False
#     return True

def ex1_jacobi(AMatrix, max_iterations=10000):
    if not checkSymmetry(AMatrix):
        raise ValueError("Matricea trebuie sa fie simetrica pentru metoda Jacobi.")
    
    iterations = 0; U = np.eye(AMatrix.shape[0], dtype=np.float64);
    
    Ainit = AMatrix.copy()

    p, q = computePQ(AMatrix)
    c, s, _ = computeCST(AMatrix, p, q)

    RMatrix = np.eye(AMatrix.shape[0], dtype=np.float64);

    # while not ( checkDiagonality(AMatrix) or abs(AMatrix[p, q]) < programEpsilon) and iterations < max_iterations:
    while abs(AMatrix[p, q]) > programEpsilon and iterations < max_iterations:
        RMatrix = rotationMatrix(AMatrix.shape[0], p, q, c, s)
        AMatrix = RMatrix @ AMatrix @ RMatrix.T

        U = U @ RMatrix.T

        p, q = computePQ(AMatrix)
        if abs(AMatrix[p, q]) <= programEpsilon:
            break

        c, s, _ = computeCST(AMatrix, p, q)
        iterations += 1

    eigenvalues = np.diag(AMatrix)
    Lambda = np.diag(eigenvalues)

    jacobiErrorFro = np.linalg.norm(Ainit @ U - U @ Lambda, ord='fro')
    # jacobiError1 = np.linalg.norm(Ainit @ U - U @ Lambda, ord=1)
    print("Valorile proprii aproximative:\n", eigenvalues)
    print("Vectorii proprii aproximativi U:\n", U)
    print("Eroare Frobenius: ", jacobiErrorFro) 
    # print("Eroare 1: ", jacobiError1) 


def libComputeCholeskyLLT(AMatrix):
    L = cholesky(AMatrix, lower=True)
    return L

def generatePositiveDefiniteMatrix(size):
    M = np.random.rand(size, size).astype(np.float64) * 31.0
    A = M.T @ M
    A += size * np.eye(size)
    return A

def checkPositiveDefinite(AMatrix):
    eigvalues = np.linalg.eigvals(AMatrix)
    return np.all(eigvalues > 0)

def ex2_CholeskyMatrixSequence(AMatrix, max_iterations=10000):
    if not checkSymmetry(AMatrix):
        raise ValueError("Matricea trebuie sa fie simetrica pentru metoda Cholesky.")
    
    iterations = 1;
    
    # Ainit = AMatrix.copy()
    Aprev = AMatrix.copy()
    L = libComputeCholeskyLLT(AMatrix)
    while iterations < max_iterations:
        L = libComputeCholeskyLLT(Aprev)
        Acurrent = L.T @ L

        if np.linalg.norm(Acurrent - Aprev, ord='fro') < programEpsilon:
            break

        Aprev = Acurrent.copy()
        iterations += 1

    print("Matricea finala A: \n", Acurrent)
    print("Valorile proprii aproximative: \n", np.diag(Acurrent))
    print("Numar de iteratii: ", iterations)

def generateRandomMatrix(p, n):
    return np.random.rand(p, n).astype(np.float64) * 31.0

def ex3_svd(AMatrix):
    U, singular_values, Vt = np.linalg.svd(AMatrix, full_matrices=True, compute_uv=True)

    print("Valorile singulare ale matricei A: ", singular_values)

    positive_sigmas = singular_values[singular_values > programEpsilon] # e bine asa sau > 0.0?
    
    high_sigma = positive_sigmas.max() if positive_sigmas.size > 0 else 0.0;
    low_sigma = positive_sigmas.min() if positive_sigmas.size > 0 else 0.0;
    Arank = positive_sigmas.size;

    print("Rangul matricei A: ", Arank)

    condition_number = np.inf if low_sigma == 0 else high_sigma / low_sigma

    print("Numarul de conditionare al matricei A: ", condition_number)

    p, n = AMatrix.shape
    SIMatrix = np.zeros((n, p), dtype=np.float64) 

    for i, sigma in enumerate(singular_values):
        if sigma > programEpsilon: # sau > 0.0?
            SIMatrix[i, i] = 1.0 / sigma
        else:
            SIMatrix[i, i] = 0.0
        
    AIMatrix = Vt.T @ SIMatrix @ U.T

    print("Pseudoinversa Moore-Penrose A_I:\n", AIMatrix)

    AJMatrix = np.linalg.inv(AMatrix.T @ AMatrix) @ AMatrix.T

    norm_diff = np.linalg.norm(AIMatrix - AJMatrix, ord=1)
    print("Norma ||A_I - A_J||_1: ", norm_diff)    



def main():
    p = getIntInput("Dimensiunea matricei p=");
    n = getIntInput("Dimensiunea matricei n=");
    epsPower = getIntInput("Epsilon=10^(-precision), precision=");
    global programEpsilon
    programEpsilon = 10 ** (-epsPower)

    symmetricA = generateSymmetricMatrix(p)
    print("matrice generata simetric: ", checkSymmetry(symmetricA))
    print();
    print("Matricea A:\n", symmetricA)

    ex1_jacobi(symmetricA)

    positiveDefiniteA = generatePositiveDefiniteMatrix(p)
    print("matrice generata pozitiv definita: ", checkPositiveDefinite(positiveDefiniteA))
    print();
    print("Matricea A:\n", positiveDefiniteA)
    print();

    ex2_CholeskyMatrixSequence(positiveDefiniteA)

    if p <= n:
        print("Ex. 3 se ruleaza doar pentru p > n.")
        return;
    A3 = generateRandomMatrix(p, n)
    print("Matricea A:\n", A3)
    ex3_svd(A3)
    

if __name__ == "__main__":
    main();