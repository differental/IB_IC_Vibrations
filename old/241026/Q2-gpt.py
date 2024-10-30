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
force_amplitude = 1.0  # amplitude of harmonic force on the first floor

# Natural frequencies of the original structure without absorbers
natural_frequencies = np.array([3.3933, 9.5078, 13.7391])

# Frequency range for harmonic response analysis
freq_range = np.linspace(0.1, 20, 500)  # Frequencies from 0.1 to 20 Hz

# Set up plot
plt.figure(figsize=(10, 6))

for n in [0, 3, 300]:  # Different numbers of absorbers
    max_response = []

    for absorber_total_mass in np.arange(0.1, 1, 0.1):
        # Define absorber mass and stiffness arrays
        absorber_m = np.full(n, (absorber_total_mass / n) if n else 0)  # absorber masses
        absorber_location = np.full(n, 3)  # all absorbers on floor 3
        absorber_k = np.zeros(n)  # absorber stiffness

        # Tune absorbers to natural frequencies
        for i in range(n):
            selected_frequency = natural_frequencies[i % 3]
            absorber_k[i] = (2 * np.pi * selected_frequency)**2 * m

        # Mass matrix
        M = m * np.eye(N + n)
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

        # Include absorbers in stiffness matrix
        for i in range(n):
            loc = absorber_location[i] - 1  # floor index for absorber
            K[loc, loc] += absorber_k[i]
            K[N + i, N + i] = absorber_k[i]
            K[loc, N + i] = -absorber_k[i]
            K[N + i, loc] = -absorber_k[i]

        # Force vector (applied only on the first floor)
        F = np.zeros(N + n)
        F[0] = force_amplitude

        # Harmonic response analysis over frequency range
        response_amplitudes = []
        for omega in 2 * np.pi * freq_range:  # angular frequency
            # Compute dynamic displacement
            try:
                X = np.linalg.solve(K - (omega**2) * M, F)
            except np.linalg.LinAlgError:
                X = np.zeros(N + n)  # In case of singular matrix issues

            # Store the response amplitude on the third floor
            response_amplitudes.append(np.abs(X[2]))

        # Track the maximum response amplitude
        max_response.append(max(response_amplitudes))

    # Plot maximum response vs. absorber mass
    plt.plot(np.arange(0.1, 1, 0.1), max_response, label=f"{n} absorbers")

# Finalize plot
plt.title("Maximum Harmonic Response Amplitude at Floor 3")
plt.xlabel("Total Absorber Mass (kg)")
plt.ylabel("Max Response Amplitude")
plt.legend()
plt.grid(True)
plt.show()
