import random

N = 2300

matrix = list()
for i in range(N):
	matrix.append([0] * N)

with open("./data.txt", "w") as file:
	for i in range(N):
		for j in range(N):
			if j < i:
				continue
			randFloat = random.uniform(-100.0, 100.0)
			randFloat = round(randFloat, 2)
			if i == j:
				randFloat += 200.0
			matrix[i][j] = randFloat
			matrix[j][i] = randFloat
	for i in range(N):
		for j in range(N):
			file.write(str(matrix[i][j]) + " ")
		file.write("\n")

	for j in range(N):
		randFloat = random.uniform(-100.0, 100.0)
		randFloat = round(randFloat, 2)
		file.write(str(randFloat) + " ")
