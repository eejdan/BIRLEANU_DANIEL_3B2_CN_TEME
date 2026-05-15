import numpy as np

def getIntInput(prompt):
	return int(input(prompt))

def generateRandomPoint(n, low=-10, high=10):
	return np.array([np.random.uniform(low, high) for _ in range(n)], dtype=np.float64)

programEpsilon = 1e-8
programFunction = lambda x: np.sum(x**2)
programGradFns = []

def constantLearningRate(_, _1, currentLearningRate = None):
	return 0.001;

def approximateLearningRate(currentPoint, gradientFn, currentLearningRate=0.1, beta=0.8):
	miu = 1;
	p = 1;
	gradient = gradientFn(currentPoint);
	while programFunction(currentPoint - miu * gradient) \
		> programFunction(currentPoint) - 0.5 * miu * np.linalg.norm(gradient)**2 \
		and p < 8:
		miu *= beta
		p = p + 1;

	return miu;

def approximateGradient(point, h=1e-5):
	gradient = np.zeros(len(point))

	for i in range(len(point)):
		ep = np.zeros_like(point)
		ep[i] = 2 * h;
		f1_x = programFunction(point + ep)
		ep[i] = h;
		f2_x = programFunction(point + ep)
		ep[i] = -h;
		f3_x = programFunction(point + ep)
		ep[i] = -2 * h;
		f4_x = programFunction(point + ep)

		gradient[i] = (-f1_x + 8 * f2_x - 8 * f3_x + f4_x) / (12 * h)

	return gradient

def analyticalGradient(point):
	return np.array([g(point) for g in programGradFns])
	


def descent(initialPoint, gradientFn, learnCalcFn, iterationLimit=30000):
	iteration = 0
	currentPoint = initialPoint.copy()
	currentLearningRate = learnCalcFn(currentPoint, gradientFn)

	cmpTerm = float('inf');

	while True:
		gradient = gradientFn(currentPoint)
		currentLearningRate = learnCalcFn(currentPoint, gradientFn, currentLearningRate)

		nextPoint = currentPoint - currentLearningRate * gradient
		currentPoint = nextPoint
		iteration += 1

		cmpTerm = currentLearningRate * np.linalg.norm(gradient)
		if not (cmpTerm >= programEpsilon and iteration <= iterationLimit and cmpTerm <= 1e10):
			break

	if cmpTerm <= programEpsilon:
		return currentPoint, iteration, True
	else:
		return currentPoint, iteration, False

def isDistinct(point1, point2, epsilon):
	return np.linalg.norm(point1 - point2) > epsilon

def addIfDistinct(points, newPoint, epsilon):
    for p in points:
        if np.linalg.norm(p - newPoint) <= epsilon:
            return
    points.append(newPoint)

def sigma(z):
    return 1 / (1 + np.exp(-z))

def main():
	testMathFunctions = [
		(
			2,
			lambda x: -np.log(1 - sigma(x[0] - x[1])) - np.log(sigma(x[0] + x[1])),
			[
				lambda x: sigma(x[0] - x[1]) + sigma(x[0] + x[1]) - 1,
				lambda x: sigma(x[0] + x[1]) - sigma(x[0] - x[1]) - 1
			]
		),(
			2,
			lambda x: x[0] ** 2 + x[1] ** 2 - 2 * x[0] - 4 * x[1] - 1,
			[
				lambda x: 2 * x[0] - 2,
				lambda x: 2 * x[1] - 4
			]
		),

		(
			2,
			lambda x: 3 * x[0] ** 2 - 12 * x[0] + 2 * x[1] ** 2 + 16 * x[1] - 10,
			[
				lambda x: 6 * x[0] - 12,
				lambda x: 4 * x[1] + 16
			]
		),

		(
			2,
			lambda x: x[0] ** 2 - 4 * x[0] * x[1] + 4.5 * x[1] ** 2 - 4 * x[1] + 3,
			[
				lambda x: 2 * x[0] - 4 * x[1],
				lambda x: -4 * x[0] + 9 * x[1] - 4
			]
		),

		(
			2,
			lambda x: x[0] ** 2 * x[1] - 2 * x[0] * x[1] ** 2 + 3 * x[0] * x[1] + 4,
			[
				lambda x: 2 * x[0] * x[1] - 2 * x[1] ** 2 + 3 * x[1],
				lambda x: x[0] ** 2 - 4 * x[0] * x[1] + 3 * x[0]
			]
		),
	]
	global programEpsilon
	programEpsilon = 10 ** (-getIntInput("Epsilon=10^(-precizie), precizie="));


	configurations_default = [
		('Constant Learning Rate, Analytical Gradient', constantLearningRate, analyticalGradient),
		('Constant Learning Rate, Approximate Gradient', constantLearningRate, approximateGradient),

		('Approximate Learning Rate, Analytical Gradient', approximateLearningRate, analyticalGradient),
		('Approximate Learning Rate, Approximate Gradient', approximateLearningRate, approximateGradient)
	]
	trial_runs = 20
	for testFn in testMathFunctions:

		constant_analytical_points = []
		constant_analytical_avg_iters = 0
		constant_analytical_convergence = 0

		constant_approximate_points = []
		constant_approximate_avg_iters = 0
		constant_approximate_convergence = 0

		approximate_analytical_points = []
		approximate_analytical_avg_iters = 0
		approximate_analytical_convergence = 0

		approximate_approximate_points = []
		approximate_approximate_avg_iters = 0
		approximate_approximate_convergence = 0

		global programFunction
		programFunction = testFn[1]
		global programGradFns
		programGradFns = testFn[2]

		for _ in range(trial_runs):
			point = generateRandomPoint(testFn[0])
			res_point, res_iters, res_convergence = descent(point, analyticalGradient, constantLearningRate)
			constant_analytical_avg_iters += res_iters
			if res_convergence:
				addIfDistinct(constant_analytical_points, res_point, programEpsilon)
				constant_analytical_convergence += 1

			res_point, res_iters, res_convergence = descent(point, approximateGradient, constantLearningRate)
			constant_approximate_avg_iters += res_iters
			if res_convergence:
				addIfDistinct(constant_approximate_points, res_point, programEpsilon)
				constant_approximate_convergence += 1

			res_point, res_iters, res_convergence = descent(point, analyticalGradient, approximateLearningRate)
			approximate_analytical_avg_iters += res_iters
			if res_convergence:
				addIfDistinct(approximate_analytical_points, res_point, programEpsilon)
				approximate_analytical_convergence += 1

			res_point, res_iters, res_convergence = descent(point, approximateGradient, approximateLearningRate)
			approximate_approximate_avg_iters += res_iters
			if res_convergence:
				addIfDistinct(approximate_approximate_points, res_point, programEpsilon)
				approximate_approximate_convergence += 1


		constant_analytical_convergence /= trial_runs
		constant_approximate_convergence /= trial_runs
		approximate_analytical_convergence /= trial_runs
		approximate_approximate_convergence /= trial_runs

		print("Results for function: ", testFn[1]);
		print(f"Ct Lr, An Gr \t\t| Ct Lr, Ap Gr \t\t| Ap Lr, An Gr \t\t| Ap Lr, Ap Gr")
		print("Convergence")
		print(f"{constant_analytical_convergence:.2f} \t\t\t| {constant_approximate_convergence:.2f} \t\t\t| {approximate_analytical_convergence:.2f} \t\t\t| {approximate_approximate_convergence:.2f}")
		print("Average iterations")
		print(f"{constant_analytical_avg_iters / trial_runs:.2f} \t\t| {constant_approximate_avg_iters / trial_runs:.2f} \t\t| {approximate_analytical_avg_iters / trial_runs:.2f} \t\t| {approximate_approximate_avg_iters / trial_runs:.2f}")





if __name__ == "__main__":
	main()