import numpy as np
from pathlib import Path

def getIntInput(prompt):
	return int(input(prompt))

def loadRareMatrixSystem(fileno, epsilon):
	data_dir = Path(__file__).resolve().parent / "data"

	def load_vector(prefix):
		file_path = data_dir / f"{prefix}_{fileno}.txt"
		if not file_path.exists():
			raise FileNotFoundError(f"Fisier inexistent: {file_path}")
		values = np.loadtxt(file_path, dtype=np.float64)
		return np.atleast_1d(values)

	b = load_vector("b")
	d0 = load_vector("d0")
	d1 = load_vector("d1")
	d2 = load_vector("d2")

	n = b.size
	if d0.size != n:
		raise ValueError("d0 trebuie sa aiba aceeasi dimensiune ca b")

	p = n - d1.size
	q = n - d2.size

	if p < 1 or p >= n:
		raise ValueError(f"Offset invalid p={p} pentru fileno={fileno}")
	if q < 1 or q >= n:
		raise ValueError(f"Offset invalid q={q} pentru fileno={fileno}")
	return {
		"epsilon": epsilon,
		"n": n,
		"p": p,
		"q": q,
		"b": b,
		"d0": d0,
		"d1": d1,
		"d2": d2,
	}

def ex1(system):
	print("Dimensiunea sistemului: n=", system["b"].size)
	
def ex2(system):
	print("p=", system["p"])
	print("q=", system["q"])

def ex3(system):
	# Verificat¸i c˘a toate elementele din d0 sunt nenule.
	for i in range(system["d0"].size):
		if abs(system["d0"][i]) < system["epsilon"]:
			raise ValueError(f"Elementul d0[{i}] este prea aproape de 0, consideram nul.");
			return
	print("Toate elementele din d0 sunt nenule.")

def ex4(system):	
	n = system["n"]
	p = system["p"]
	q = system["q"]
	b = system["b"]
	d0 = system["d0"]
	d1 = system["d1"] # +p a[i,i+p] = d1[i]
	d2 = system["d2"] # -q a[i,i-q] = d2[i-q]
	epsilon = system["epsilon"]

	x_current = np.zeros(system["n"]).astype(np.float64)
	x_prev = x_current.copy()
	
	#verifica ca nu exista elemente nule pe diag principala
	# ex3(system)
	iterationsLeft = 100
	while iterationsLeft > 0:
		x_prev = x_current.copy()
		for i in range(system["n"]):
			x_current[i] = b[i]
			if i - p >= 0:
				x_current[i] -= d1[i - p] * x_current[i - p]
			if i + p < n:
				x_current[i] -= d1[i] * x_prev[i + p]
			if i - q >= 0:
				x_current[i] -= d2[i - q] * x_current[i - q]
			if i + q < n:
				x_current[i] -= d2[i] * x_prev[i + q]
			x_current[i] /= d0[i]

		manhattan_distance = 0.0;
		for i in range(n):
			manhattan_distance += abs(x_current[i] - x_prev[i])
		if manhattan_distance < epsilon:
			print("Convergent in", 100 - iterationsLeft, "iteratii.")
			return x_current
		iterationsLeft -= 1


	raise ValueError("Metoda Gauss-Seidel este divergenta.")
	

def ex5(system, xgs):
	d0 = system["d0"]
	d1 = system["d1"]
	d2 = system["d2"]
	n = system["n"]
	p = system["p"]
	q = system["q"]

	y = np.zeros(n)
	for i in range(n):
		y[i] = d0[i] * xgs[i]
		if i - p >= 0:
			y[i] += d1[i - p] * xgs[i - p]
		if i + p < n:
			y[i] += d1[i] * xgs[i + p]
		if i - q >= 0:
			y[i] += d2[i - q] * xgs[i - q]
		if i + q < n:
			y[i] += d2[i] * xgs[i + q]

	return y;

def ex6(system, y):
	b = system["b"]
	n = system["n"]

	maxval = 0;
	
	for i in range(n):
		if abs(y[i] - b[i]) > maxval:
			maxval = abs(y[i] - b[i])
	
	print("Norma infinita a erorii: ", maxval)
	


	
# def saveVector(filename, vector):
# 	np.savetxt(filename, vector, fmt='%.6f')

def main():
	epsilon = 10**(-getIntInput(f"eps=10^(-p), p="))
	system = None
	for fileno in [1,2,3,4,5]:
		system = loadRareMatrixSystem(fileno, epsilon)
		
		ex1(system) # dimensiune
		ex2(system) # ordinul p si q
		try:
			ex3(system) # verificare d0 nenule
			xgs = ex4(system) # Ax=B cu GaussSeidel
			print("Solutia aproximativa xgs:")
			print(xgs)
			#saveVector(f"result/xgs_{fileno}.txt", xgs)
			
			y = ex5(system, xgs) # y = Axgs
			print("Produsul y=Axgs:")
			print(y)
			#saveVector(f"result/y_{fileno}.txt", y)

			ex6(system, y) # eroare maxima Linf
		except ValueError as e:
			print(e)

if __name__ == "__main__":
	main()