import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

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
n = 30  # number of absorbers
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

# Symbolic displacement
w = sp.symbols('w')
M_sym = sp.Matrix(M)
K_sym = sp.Matrix(K)

# Force vector for harmonic analysis
F = np.zeros(N + n)
F[0] = 1  # unit force at floor 1

F_sym = sp.Matrix(F)

# Solve for symbolic displacement
B_sym = K_sym - w**2 * M_sym
disp_sym = B_sym.inv() * F_sym  # symbolic displacement

# Plot the symbolic displacement for each floor
frequency_range = np.linspace(0.1, 130, 1000)  # avoid w=0 to prevent division by zero
disp_func = [sp.lambdify(w, disp_sym[i], 'numpy') for i in range(N)]  # displacement functions for each floor

# Plotting the symbolic expressions over the frequency range
plt.figure(figsize=(10, 6))

# Loop through each floor and plot
colors = ['b', 'g', 'r']
for i in range(N):
    disp_vals = disp_func[i](frequency_range)
    plt.plot(frequency_range, np.abs(disp_vals), label=f'Floor {i+1}', color=colors[i])

plt.xlabel('Frequency (Hz)')
plt.ylabel('Displacement')
plt.title('Frequency Response of Each Floor (Symbolic)')
plt.legend()
plt.grid(True)
plt.show()
