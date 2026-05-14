import numpy as np

programEpsilon = 1e-8

def getIntInput(prompt):
    return int(input(prompt))

def HornerPolynomial(coefficients, x):
    result = np.zeros_like(x, dtype=np.float64)
    for k in range(len(coefficients)):
        result = result * x + coefficients[k]
    return result

def computeDerivativeCoefficients(coefficients):
    n = len(coefficients) - 1
    derivative_coefficients = np.zeros(n, dtype=np.float64)
    for i in range(n):
        derivative_coefficients[i] = coefficients[i] * (n - i)
    return derivative_coefficients

def IsDistinct(v1, v2):
    return abs(v1 - v2) > programEpsilon

def writeVectorToFile(filename, vector):
    if len(vector) == 0:
        return;
    with open(filename, 'w') as f:
        f.write(f"{vector[0]} ")
        for i in range(1, len(vector)):
            if IsDistinct(vector[i], vector[i-1]):
                f.write(f"{vector[i]} ")

def getRandomPointInInterval(low, high):
    return np.random.uniform(min(low, high), max(low, high))

def metNewton(coefficients, x, max_iterations=10000):
    iterations = 0;
    first_derivate_coefficients = computeDerivativeCoefficients(coefficients)
    
    while iterations < max_iterations:
        f_x = HornerPolynomial(coefficients, x)
        f_prime_x = HornerPolynomial(first_derivate_coefficients, x)

        if abs(f_prime_x) < programEpsilon:
            break

        deltax = f_x / f_prime_x

        x = x - deltax;

        if abs(deltax) < programEpsilon:
            return x, iterations, True
        
        if abs(deltax) > 10**8:
            break
        iterations += 1

    return x, iterations, False

def metOlver(coefficients, x, max_iterations=10000):
    iterations = 0;
    first_derivate_coefficients = computeDerivativeCoefficients(coefficients)
    second_derivate_coefficients = computeDerivativeCoefficients(first_derivate_coefficients)
    
    while iterations < max_iterations:
        f_x = HornerPolynomial(coefficients, x)
        f_prime_x = HornerPolynomial(first_derivate_coefficients, x)
        f_double_prime_x = HornerPolynomial(second_derivate_coefficients, x)

        if abs(f_prime_x) < programEpsilon:
            break

        c = ((f_x ** 2) * f_double_prime_x) / (f_prime_x ** 3) 
        deltax = f_x / f_prime_x + 0.5 * c
        x = x - deltax
        if abs(deltax) < programEpsilon:
            return x, iterations, True
        
        if abs(deltax) > 10**8:
            break

        iterations += 1

    return x, iterations, False

def main():
    coefficients = [1.0, -6.0, 11.0, -6.0]
    global programEpsilon
    programEpsilon = 10**(-getIntInput("epsilon=10**(-precizie), precizie="))
    numar_puncte = getIntInput("numar puncte x0=")

    root_radius = (abs(coefficients[0]) + np.max(np.abs(coefficients[1:]))) / abs(coefficients[0])

    distinct_roots = [];

    for i in range(numar_puncte):
        x0 = np.float64(getRandomPointInInterval(-root_radius, root_radius))
        root_newton, iter_newton, converged_newton = metNewton(coefficients, x0)
        root_olver, iter_olver, converged_olver = metOlver(coefficients, x0)

        print(f"Initial point: {x0:.6f}")
        print(f"Newton's method: Root = {root_newton:.6f}, Iterations = {iter_newton}, Converged = {converged_newton}")
        print(f"Olver's method: Root = {root_olver:.6f}, Iterations = {iter_olver}, Converged = {converged_olver}")
        print("-" * 40)

        if converged_newton:
            if all(IsDistinct(root_newton, r) for r in distinct_roots):
                distinct_roots.append(root_newton)
        
        if converged_olver:
            if all(IsDistinct(root_olver, r) for r in distinct_roots):
                distinct_roots.append(root_olver)

    print(f"Distinct roots found: {len(distinct_roots)}")
    for r in distinct_roots:
        print(f"{r:.6f}")
    
    writeVectorToFile("roots.txt", distinct_roots);


if __name__ == "__main__":
    main()