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

plot_colours = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']

all_all_disp = []


plt.figure(figsize=(10, 6))




'''
# Damped natural frequencies in radians
damped_natural_frequencies = np.array([3.16027156, 9.08205525, 13.51090025]) * 2 * np.pi

all_good_frequencies = np.array([3.16, 3.39, 9.08, 9.51, 13.51, 13.74]) * 2 * np.pi

# Initialize an empty list to store all frequency values
frequencies = []

# Define the frequency range
start_freq = 1
end_freq = 131

# Add fine intervals around each damped natural frequency
for freq in all_good_frequencies:
    fine_range = np.arange(freq - 0.3, freq + 0.3, 0.01)
    frequencies.extend(fine_range)

# Add broader intervals outside the fine intervals
outside_range = np.arange(start_freq, end_freq, 0.3)
frequencies.extend(outside_range)

# Convert to a unique sorted array
frequencies = np.unique(np.array(frequencies))
frequencies.sort()
'''

frequencies = np.arange(1, 130, 0.005)





for idx, n in enumerate([0, 3]):

    absorber_total_mass = 1.3
    absorber_m = np.full(n, (absorber_total_mass / n) if n else 0)  # array of absorber mass
    absorber_zeta = 0.3
    absorber_location = np.full(n, 3)  # all absorbers located at floor 3
    structure_zeta = 0.01

    # Initialize absorber stiffness array
    absorber_k = np.zeros(n)

    # Tuning absorbers to match natural frequencies
    for i in range(n):
        #selected_frequency = natural_frequencies[i % 3]
        selected_frequency = 3.3933
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

    # Damping factor matrix
    L = np.zeros((N + n, N + n))

    K[0, 0] = 2 * k
    K[0, 1] = -k
    K[1, 0] = -k
    K[1, 1] = 2 * k
    K[2, 2] = k
    K[1, 2] = -k
    K[2, 1] = -k
    
    
    structure_damping_ratio = 2 * np.sqrt(m * k) * structure_zeta
    
    L[0, 0] = 2 * structure_damping_ratio
    L[0, 1] = -structure_damping_ratio
    L[1, 0] = -structure_damping_ratio
    L[1, 1] = 2 * structure_damping_ratio
    L[2, 2] = structure_damping_ratio
    L[1, 2] = -structure_damping_ratio
    L[2, 1] = -structure_damping_ratio

    # Modify stiffness matrix to include absorbers
    for i in range(n):
        loc = absorber_location[i] - 1  # adjust for Python's 0-indexing

        K[loc, loc] += absorber_k[i]
        K[3 + i, 3 + i] = absorber_k[i]
        K[loc, 3 + i] = -absorber_k[i]
        K[3 + i, loc] = -absorber_k[i]

        absorber_damping_ratio = 2 * np.sqrt(absorber_m[i] * absorber_k[i]) * absorber_zeta

        L[loc, loc] += absorber_damping_ratio
        L[3 + i, 3 + i] = absorber_damping_ratio
        L[loc, 3 + i] = -absorber_damping_ratio
        L[3 + i, loc] = -absorber_damping_ratio

    # Solve the eigenvalue problem
    #eigenvalues, eigenvectors = np.linalg.eig(np.linalg.inv(M) @ K)


    #values = np.abs(eigenvectors[0:3][0:3])
    #print(f"Absorbers: {n}, Avg: {np.sum(values[0:3]) / 9}, Max: {np.max(values[0:3])}\n")

    # Calculate natural frequencies in Hz
    #natural_frequencies_with_absorbers = np.sqrt(np.real(eigenvalues)) / (2 * np.pi)
    #print("Natural Frequencies (Hz):", natural_frequencies_with_absorbers)

    # Force vector for harmonic analysis
    F = np.zeros(N + n)
    F[0] = 1  # unit force at floor 1

    # Frequency response functions
    all_disp = []

    for w in frequencies:
        B = - (w**2) * M + 1j * w * L + K
        disp = np.linalg.solve(B, F)
        all_disp.append(np.abs(disp[:N]))  # store displacements of the floors
        #all_disp.append(disp[:N])

    all_disp = np.array(all_disp).T

    # Plot frequency response for each floor
    #plt.figure(figsize=(10, 6))

    # Plot for each floor
    #for i in range(N):

    all_all_disp.append((all_disp[0], all_disp[1], all_disp[2], n, idx))

    #plt.plot(frequencies, np.log(all_disp[floor]), label=f'{n} Absorber', color=plot_colours[idx % len(plot_colours)])

for floor in range(3):
    for item in all_all_disp:
        plt.plot(frequencies, np.log(item[floor]), label=f'{item[3]} Absorber(s)', color=plot_colours[item[4] % len(plot_colours)])

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Log of Displacement Magnitude')
    plt.xlim(1, 131)
    #plt.ylim(-15, 5)
    plt.title(f"Frequency Response of Floor {floor+1}")

    plt.legend()
    plt.grid(True)
    plt.show()

