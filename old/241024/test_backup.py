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
n = 1000  # number of absorbers
absorber_total_mass = 0.3
absorber_m = np.full(n, absorber_total_mass / n)  # array of absorber mass
absorber_location = np.full(n, 3)  # all absorbers located at floor 3

# Initialize absorber stiffness array
absorber_k = np.zeros(n)

# Tuning absorbers to match natural frequencies
for i in range(n):
    selected_frequency = natural_frequencies[i % 3]
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

# Calculate natural frequencies in Hz
natural_frequencies_with_absorbers = np.sqrt(np.real(eigenvalues)) / (2 * np.pi)
print("Natural Frequencies (Hz):", natural_frequencies_with_absorbers)

# Force vector for harmonic analysis
F = np.zeros(N + n)
F[0] = 1  # unit force at floor 1

# Frequency response functions
frequencies = np.arange(1, 131)
all_disp = []

for w in frequencies:
    B = K - (w**2) * M
    disp = np.linalg.solve(B, F)
    #all_disp.append(np.abs(disp[:N]))  # store displacements of the floors
    all_disp.append(disp[:N])

all_disp = np.array(all_disp).T

# Plot frequency response for each floor
plt.figure(figsize=(10, 6))

# Plot for each floor
for i in range(N):
    plt.plot(frequencies, all_disp[i], label=f'Floor {i+1}')

plt.xlabel('Frequency (Hz)')
plt.ylabel('Displacement')
plt.title('Frequency Response of Each Floor')
plt.legend()
plt.grid(True)
plt.show()

