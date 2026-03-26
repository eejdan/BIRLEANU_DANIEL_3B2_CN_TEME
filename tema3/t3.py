import numpy as np
from scipy.linalg import qr

def getIntInput(prompt):
    return int(input(prompt))

def matrixCreate(size):
    return np.random.randint(0, 31, size=(size, size)).astype(float)

def vectorCreate(size):
    return np.random.randint(0, 31, size=size).astype(float)

def main():
    matrix_size = getIntInput("Size: ");print();
    compute_precision = getIntInput("Precision factor: ");print();

    instanceMatrix = matrixCreate(matrix_size)
    instanceVector = vectorCreate(matrix_size)

    print("Instance Matrix:\n", instanceMatrix)
    print("Instance Vector:\n", instanceVector)

    epsilon = 10**(-compute_precision)
    print("Epsilon:", epsilon)
    ex1(instanceMatrix, instanceVector, epsilon)
    ex2(instanceMatrix, instanceVector, epsilon)
    ex3(instanceMatrix, instanceVector, epsilon)
    ex4(instanceMatrix, instanceVector, epsilon)
    ex5(instanceMatrix, instanceVector, epsilon)


def ex1(aMatrix, sVector, epsilon):
    print("\n\nEx 1")
    bVector = computeBVector(aMatrix, sVector)
    print("bVector:\n", bVector)

def ex2(aMatrix, sVector, epsilon):
    print("\n\nEx 2")
    bVector = computeBVector(aMatrix, sVector)

    qMatrix, rMatrix, _ = householder(aMatrix.copy(), bVector.copy(), epsilon)

    orthogonality = np.allclose(qMatrix.T @ qMatrix, np.eye(qMatrix.shape[0]), atol=1e-6)
    triangular = is_upper_triangular(rMatrix)
    reconstruction = np.allclose(aMatrix, qMatrix @ rMatrix, atol=1e-6)

    print(f"Orthogonal Q:        {orthogonality}")
    print(f"Upper triangular R:  {triangular}")
    print(f"A ≈ Q @ R:           {reconstruction}")
    
def ex3(aMatrix, sVector, epsilon):
    print("\n\nEx 3")
    bVector = computeBVector(aMatrix, sVector)

    qLib, rLib = qr(aMatrix)
    xQR = np.linalg.solve(rLib, qLib.T @ bVector)

    _, rHouseholder, bHouseholder = householder(aMatrix.copy(), bVector.copy(), epsilon)
    xHouseholder = backward_substitution(rHouseholder, bHouseholder, epsilon)

    diff_norm = np.linalg.norm(xQR - xHouseholder, ord=2)

    print("xQR:\n", xQR)
    print("xHouseholder:\n", xHouseholder)
    print("||xQR - xHouseholder||2 =", diff_norm)
    
def ex4(aMatrix, sVector, epsilon):
    print("\n\nEx 4")
    bVector = computeBVector(aMatrix, sVector)

    qLib, rLib = qr(aMatrix)
    xQR = np.linalg.solve(rLib, qLib.T @ bVector)

    _, rHouseholder, bHouseholder = householder(aMatrix.copy(), bVector.copy(), epsilon)
    xHouseholder = backward_substitution(rHouseholder, bHouseholder, epsilon)

    diff_norm = np.linalg.norm(aMatrix @ xHouseholder - bVector, ord=2)
    diff2_norm = np.linalg.norm(aMatrix @ xQR - bVector, ord=2)
    diff3_norm = np.linalg.norm(xHouseholder - sVector, ord=2) / np.linalg.norm(sVector, ord=2)
    diff3_norm = np.linalg.norm(xQR - sVector, ord=2) / np.linalg.norm(sVector, ord=2)


    print("xQR:\n", xQR)
    print("xHouseholder:\n", xHouseholder)
    print("||AxHouseholder - b||2 =", diff_norm)
    print("||AxQR - b||2 =", diff2_norm)
    print("||xHouseholder - s||2 / ||s||2 =", diff3_norm)
    print("||xQR - s||2 / ||s||2 =", diff3_norm)

def ex5(aMatrix, sVector, epsilon):
    print("\n\nEx 5");
    inverseHouseholder = computeInverse(aMatrix, sVector, epsilon)
    inverseLib = np.linalg.inv(aMatrix)

    print("Inverse Householder:\n", inverseHouseholder)
    print("Inverse Library:\n", inverseLib)
    
    diff_norm = np.linalg.norm(inverseHouseholder - inverseLib, ord=2)
    print("||InverseHouseholder - InverseLib||2 =", diff_norm)


def computeInverse(aMatrix, sVector, epsilon):
    n = aMatrix.shape[0]
    zeroVector = np.zeros(n)

    qHouseholder, rHouseholder, _ = householder(aMatrix.copy(), zeroVector, epsilon)

    for i in range(n):
        if abs(rHouseholder[i, i]) < epsilon:
            raise ValueError("Matrice singulara sau aproape singulara. Inversa nu poate fi calculata.")

    inverseMatrix = np.zeros((n, n))

    for j in range(n):
        eVector = np.zeros(n)
        eVector[j] = 1.0

        bVector = qHouseholder.T @ eVector
        xVector = backward_substitution(rHouseholder, bVector, epsilon)
        inverseMatrix[:, j] = xVector

    return inverseMatrix


def computeBVector(aMatrix, sVector):
    bVector = np.zeros(sVector.shape[0])
    for i in range(bVector.shape[0]):
        for j in range(aMatrix.shape[0]):
            bVector[i] += aMatrix[i, j] * sVector[j]
    return bVector

def identityMatrix(size):
    return np.identity(size)

def backward_substitution(rMatrix, bVector, epsilon):
    n = rMatrix.shape[0]
    xVector = np.zeros(n)

    for i in range(n - 1, -1, -1):
        if abs(rMatrix[i, i]) < epsilon:
            raise ValueError("Matrice singulara sau aproape singulara in substitutia inversa")

        sum_upper = 0
        for j in range(i + 1, n):
            sum_upper += rMatrix[i, j] * xVector[j]

        xVector[i] = (bVector[i] - sum_upper) / rMatrix[i, i]

    return xVector

def householder(aMatrix, bVector, epsilon):
    
    qMatrix = identityMatrix(aMatrix.shape[0])
    
    for r in range(aMatrix.shape[0]-1):


        sigma = 0;
        for j in range(r, aMatrix.shape[0]):
            sigma += aMatrix[j, r] * aMatrix[j, r]
        
        if abs(sigma) < epsilon:
            break

        k = np.sqrt(sigma)
        if aMatrix[r, r] >= 0:
            k = -k


        beta = sigma - k * aMatrix[r, r]

        u = np.zeros(aMatrix.shape[0])
        u[r] = aMatrix[r, r] - k
        for i in range(r+1, aMatrix.shape[0]):
            u[i] = aMatrix[i, r]

        for j in range(r+1, aMatrix.shape[0]):
            gamma = 0
            for i in range(r, aMatrix.shape[0]):
                gamma += u[i] * aMatrix[i, j]
            gamma /= beta

            for i in range(r, aMatrix.shape[0]):
                aMatrix[i, j] -= gamma * u[i]
        
        # transformarea coloanei r...
        aMatrix[r, r] = k
        for i in range(r+1, aMatrix.shape[0]):
            aMatrix[i, r] = 0
        # b = Pr @ b
        gamma = 0
        for i in range(r, aMatrix.shape[0]):
            gamma += u[i] * bVector[i]
        gamma /= beta

        for i in range(r, aMatrix.shape[0]):
            bVector[i] -= gamma * u[i]

        for j in range(aMatrix.shape[0]):
            gamma = 0
            for i in range(r, aMatrix.shape[0]):
                gamma += u[i] * qMatrix[i, j]
            gamma /= beta

            for i in range(r, aMatrix.shape[0]):
                qMatrix[i, j] -= gamma * u[i]
    
    return qMatrix.T, aMatrix, bVector

def is_upper_triangular(R, tol=1e-8):
    return np.allclose(R, np.triu(R), atol=tol)

if __name__ == "__main__":
    main();