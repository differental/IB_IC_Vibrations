import numpy as np
import matplotlib.pyplot as plt

# Constants
m = 1.83  # mass of one floor
L = 0.2  # length
N = 3  # number of degrees of freedom
b = 0.08  # width
E = 210e9  # Young's Modulus
d = 0.001  # thickness
I = b * d**3 / 12  # second moment of area
k = (24 * E * I) / (L**3)  # static stiffness for each floor

# Natural frequencies of the original structure without absorbers
natural_frequencies = np.array([3.3933, 9.5078, 13.7391])

# Number of absorbers and mass distribution
#n = 1000  # number of absorbers

plt.figure(figsize=(10, 6))

for n in [0, 3, 30, 300]:
    
    results_mass = []
    results_avg = []
    results_max = []

    for absorber_total_mass in np.arange(0.005, 0.5, 0.005):

    
        #absorber_total_mass = 0.3
        absorber_m = np.full(n, (absorber_total_mass / n) if n else 0)  # array of absorber mass
        absorber_location = np.full(n, 3)  # all absorbers located at floor 3

        # Initialize absorber stiffness array
        absorber_k = np.zeros(n)

        # Tuning absorbers to match natural frequencies
        for i in range(n):
            selected_frequency = natural_frequencies[i % 3]
            #selected_frequency = 9.5078
            absorber_k[i] = (2 * np.pi * selected_frequency)**2 * m
            r = 0.95 + (1.05 - 0.95) * np.random.rand()  # random number between 95% and 105%
            absorber_k[i] *= r

        # Mass matrix
        M = m * np.eye(N + n)

        # Assign values from absorber_m to diagonal positions starting from M[4,4]
        for i in range(n):
            M[i + N, i + N] = absorber_m[i]

        # Stiffness matrix
        K = np.zeros((N + n, N + n))

        K[0, 0] = 2 * k
        K[0, 1] = -k
        K[1, 0] = -k
        K[1, 1] = 2 * k
        K[2, 2] = k
        K[1, 2] = -k
        K[2, 1] = -k

        # Modify stiffness matrix to include absorbers
        for i in range(n):
            loc = absorber_location[i] - 1  # adjust for Python's 0-indexing
            K[loc, loc] += absorber_k[i]
            K[3 + i, 3 + i] = absorber_k[i]
            K[loc, 3 + i] = -absorber_k[i]
            K[3 + i, loc] = -absorber_k[i]

        # Solve the eigenvalue problem
        eigenvalues, eigenvectors = np.linalg.eig(np.linalg.inv(M) @ K)
        
        normalized_eigenvectors = eigenvectors / np.linalg.norm(eigenvectors, axis=0) # normalisation

        values = np.abs(normalized_eigenvectors[0:3, 0:3])
        
        print(values)
        
        print(f"Absorbers: {n}, Avg: {np.sum(values[0:3]) / 9}, Max: {np.max(values[0:3])}\n")
        #results.append((absorber_total_mass, np.sum(values[0:3]) / 9, np.max(values[0:3])))

        results_mass.append(absorber_total_mass)
        results_avg.append(np.sum(values[0:3]) / 9)
        results_max.append(np.max(values[0:3]))


    #plt.plot(results_mass, results_avg, label='Avg')
    plt.plot(results_mass, results_max, label=f"{n} absorbers")
    
plt.title(f"Harmonic Maximum Amplitude")

plt.legend()
plt.grid(True)
plt.show()
    