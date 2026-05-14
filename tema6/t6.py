import numpy as np
# from scipy.linalg import cholesky
import matplotlib.pyplot as plt

def getIntInput(prompt):
    return int(input(prompt))

def getFloatInput(prompt):
    return float(input(prompt))

def getSpaceSeparatedFloats(expected_count):
    print(f"Enter {expected_count} space-separated float values:")
    return np.array(list(map(float, input().split())))

def generateVector(size, low, high):
    return np.random.rand(size).astype(float) * (abs(high - low)) + np.average([high, low])

def generateRisingVector(size, low, high):
    return np.sort(generateVector(size, low, high))


def f(x):
    return np.exp(0.1 * x) * (x + np.sin(x))

def plot_function(func, x_start, x_end, known_x, known_y, func2=None):
    x_plot = np.linspace(x_start, x_end, 300)
    y_plot = func(x_plot)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, y_plot, 'b-', label='f(x)', linewidth=2)
    plt.scatter(known_x, known_y, color='red', label='Known points', zorder=5)
    if func2 is not None:
        y_plot2 = func2(x_plot)
        plt.plot(x_plot, y_plot2, 'g-', label='Interpolated f(x)', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Function Graph')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

def ex1(known_points_number, points, imaged_points, polynomial_degree): 
    
    B = np.zeros((known_points_number, polynomial_degree + 1))
    for i in range(known_points_number):
        for j in range(polynomial_degree + 1):
            B[i, j] = points[i] ** j

    coefficients = np.linalg.lstsq(B, imaged_points, rcond=None)[0]

    def horner_polynomial(x):
        result = np.zeros_like(x, dtype=float)
        for k in range(len(coefficients) - 1, -1, -1):
            result = result * x + coefficients[k]
        return result

    plot_function(f, points[0], points[-1], points, imaged_points, horner_polynomial)

def f_derivate1(x):
    return np.exp(0.1 * x) * (0.1 * (x + np.sin(x)) + 1 + np.cos(x))
   

def main():
    range_low = getFloatInput("x1 = a=")
    range_high = getFloatInput("xn = b=")
    known_points_number = getIntInput("Enter numar puncte stiute ale functiei, n+1=")
    
    print("Se va testa pe functia f(x) = exp(0.1 * x) * (x + sin(x))")
    points = generateRisingVector(known_points_number, range_low, range_high);
    imaged_points = f(points)
    polynomial_degree = getIntInput("Enter gradul polinomului de interpolare, m=")

    ex1(known_points_number, points, imaged_points, polynomial_degree)

    # print("Prima derivata a functiei f(x) este f'(x) = exp(0.1 * x) * (0.1 * (x + sin(x)) + 1 + cos(x))")



if __name__ == "__main__":
    main()