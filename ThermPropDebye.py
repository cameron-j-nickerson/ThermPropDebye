# Program to calculate thermal properties of acoustic phonons from elastic constants
# Cameron Nickerson
# June 2025


import numpy as np
from scipy.integrate import quad
import glob


# READ IN DATA

in_files = glob.glob("*.clmz")

if len(in_files) != 1:
    print("should be exactly one file with '.clmz' extension in current directory")
    quit()

in_file_name = in_files[0]
in_file = open(in_file_name, 'r')

C = []  # GPa
for i in range(6):
    line = [float(i) for i in in_file.readline().split()]
    C = np.append(C, line)
C = C.reshape([6,6])

in_file.readline()

lattice = []    # Angstroms
for i in range(3):
    line = [float(i) for i in in_file.readline().split()]
    lattice = np.append(lattice, line)
lattice = lattice.reshape([3,3])

in_file.readline()

mass = float(in_file.readline()) # atomic mass units

in_file.readline()

Z = int(in_file.readline()) # molecules per cell

in_file.close()


# CREATE RECIP LATTICE

# IN ANGSTROM
a = lattice[0, :]
b = lattice[1, :]
c = lattice[2, :]

v_cross = np.cross(b, c)
v_dot = np.dot(a, v_cross)
volume = v_dot * np.power(1E-10, 3) # m^3

# Get the reciprocal Lattice
recip_lattice = []
astar = (2*np.pi) / v_dot * np.cross(b, c)
bstar = (2*np.pi) / v_dot * np.cross(c, a)
cstar = (2*np.pi) / v_dot * np.cross(a, b)
recip_lattice.append(astar)
recip_lattice.append(bstar)
recip_lattice.append(cstar)
recip_lattice = np.reshape(np.array(recip_lattice), ([3, 3]))  # rads / angstrom


# constants (SI)
hbar = 1.05457182e-34
kB = 1.380649e-23
N_A = 6.02214076e23

# unit conversions
mass_kg = 1.66054e-27
Pa = 1e9

mass_tot = mass * mass_kg  # In kg
density = mass_tot / volume # kg / m^3
C = C * Pa # Conversion from GPa to Pa


# SOLVE CHRISTOFFEL EQUATIONS

Christoffel_Mat = np.zeros([3, 3])

# thirteen symmetry-unique directions sampled
kpt_sampling_directions = np.array([0.5, 0, 0,
                                    0, 0.5, 0,
                                    0, 0, 0.5,
                                    0.5, 0.5, 0,
                                    0.5, 0, 0.5,
                                    0, 0.5, 0.5,
                                    0.5, 0.5, 0.5,
                                    0.5, -0.5, 0,
                                    0.5, 0, -0.5,
                                    0, 0.5, -0.5,
                                    0.5, 0.5, -0.5,
                                    0.5, -0.5, 0.5,
                                    -0.5, 0.5, 0.5]).reshape([-1, 3])

Wmax = [] # list of Debye frequencies (13 directions times 3 modes)
soundSpeed = []
for i in range(int(np.size(kpt_sampling_directions)/3)):
    d_cos = np.zeros(3)
    recip_sampling = np.matmul(kpt_sampling_directions[i], recip_lattice)

    # Take direction cosines wrt recip lattice system
    d_cos[0] = (recip_sampling[0] / np.linalg.norm(recip_sampling))
    d_cos[1] = (recip_sampling[1] / np.linalg.norm(recip_sampling))
    d_cos[2] = (recip_sampling[2] / np.linalg.norm(recip_sampling))

    # Building the Christoffel matrix, \Gamma_ik
    A_1 = np.square(d_cos[0])*C[0, 0] + np.square(d_cos[1])*C[5, 5] + np.square(d_cos[2])*C[4, 4] + \
        2*d_cos[1]*d_cos[2]*C[4, 5] + 2*d_cos[2] * \
        d_cos[0]*C[0, 4] + 2*d_cos[0]*d_cos[1]*C[0, 5]
    A_2 = np.square(d_cos[0])*C[5, 5] + np.square(d_cos[1])*C[1, 1] + np.square(d_cos[2])*C[3, 3] + \
        2*d_cos[1]*d_cos[2]*C[1, 3] + 2*d_cos[2] * \
        d_cos[0]*C[3, 5] + 2*d_cos[0]*d_cos[1]*C[1, 5]
    A_3 = np.square(d_cos[0])*C[4, 4] + np.square(d_cos[1])*C[3, 3] + np.square(d_cos[2])*C[2, 2] + \
        2*d_cos[1]*d_cos[2]*C[2, 3] + 2*d_cos[2] * \
        d_cos[0]*C[2, 4] + 2*d_cos[0]*d_cos[1]*C[3, 4]
    alpha_23 = np.square(d_cos[0])*C[4, 5] + np.square(d_cos[1])*C[1, 3] + np.square(d_cos[2])*C[2, 3] + 2*d_cos[1]*d_cos[2]*(
        1/2*(C[1, 2] + C[3, 3])) + 2*d_cos[2]*d_cos[0]*(1/2*(C[2, 5]+C[3, 4])) + 2*d_cos[0]*d_cos[1]*(1/2*(C[1, 4]+C[3, 5]))
    alpha_13 = np.square(d_cos[0])*C[0, 4] + np.square(d_cos[1])*C[3, 5] + np.square(d_cos[2])*C[2, 4] + 2*d_cos[1]*d_cos[2]*(
        1/2*(C[2, 5] + C[3, 4])) + 2*d_cos[2]*d_cos[0]*(1/2*(C[0, 2]+C[4, 4])) + 2*d_cos[0]*d_cos[1]*(1/2*(C[0, 3]+C[4, 5]))
    alpha_12 = np.square(d_cos[0])*C[0, 5] + np.square(d_cos[1])*C[1, 5] + np.square(d_cos[2])*C[3, 4] + 2*d_cos[1]*d_cos[2]*(
        1/2*(C[1, 4] + C[3, 5])) + 2*d_cos[2]*d_cos[0]*(1/2*(C[0, 3]+C[4, 5])) + 2*d_cos[0]*d_cos[1]*(1/2*(C[0, 1]+C[5, 5]))

    Christoffel_Mat[0, 0] = A_1
    Christoffel_Mat[0, 1] = alpha_12
    Christoffel_Mat[0, 2] = alpha_13
    Christoffel_Mat[1, 0] = alpha_12
    Christoffel_Mat[1, 1] = A_2
    Christoffel_Mat[1, 2] = alpha_23
    Christoffel_Mat[2, 0] = alpha_13
    Christoffel_Mat[2, 1] = alpha_23
    Christoffel_Mat[2, 2] = A_3

    PhaseVelocity = np.linalg.eig(Christoffel_Mat)[0]  # = pv^2
    PhaseVelocity = abs(PhaseVelocity)
    kzb = np.linalg.norm(np.matmul(kpt_sampling_directions[i], recip_lattice)) * 1E10  # kzb In rad / m

    PhaseVelocity = np.transpose(PhaseVelocity)
    PhaseVelocity = PhaseVelocity / density # -- kg/(m*s^2) / kg/m^3 -> m^2/s^2
    PhaseVelocity = np.sqrt(PhaseVelocity)  # -- m*s^-1 -> Velocity
    nan_ind = np.isnan(PhaseVelocity)
    PhaseVelocity[nan_ind] = 0

    Coeff = 2 * PhaseVelocity * kzb / np.pi  # -- m/s * rad/m -> rad/s = Angular frequency
    Wmax = np.append(Wmax, Coeff)
    soundSpeed = np.append(soundSpeed, PhaseVelocity)

print("Average speed of sound:", "{:.0f}".format(np.mean(soundSpeed)), "m/s")

Wmax = Wmax.reshape([-1, 3])

w_D = np.mean(Wmax)
w_D_THz = w_D / (10**12) / (2*np.pi)

print("Debye frequency:", "{:.4f}".format(w_D_THz), "THz")

def I(t):
    return t**3 / (np.exp(t)-1)

def D(x):
    return (3/(x**3)) * quad(I, 0, x)[0]


# calculate Fvib and write to file
T_arr = np.arange(1, 301, 1) # K
out_file_name = in_file_name[0:-4] + "therm"
out_file = open(out_file_name, 'w')
out_file.write("T_(K)    Fvib_(kJ/mol)    Svib_(kJ/mol/K)    Uvib_(kJ/mol)\n")
Fvib_arr = []
S_arr = []
k=0
for T in T_arr:
    Fvib = 9/8 * hbar*w_D + 3*kB*T * np.log(1-np.exp(-hbar*w_D/(kB*T))) - kB*T * D(hbar*w_D/(kB*T)) # J
    Fvib = Fvib / 1000 * N_A / Z # kJ/mol
    Fvib_arr.append(Fvib)
    if k>0:
        dF = Fvib - Fvib_arr[k-1]
        S = -dF
        U = Fvib + T*S
        line = "{:.1f}".format(T) + "      " + "{:.4f}".format(Fvib) + "           " + "{:.4f}".format(S) + "             " + "{:.4f}".format(U) + "\n"
        out_file.write(line)
    k+=1
out_file.close()

print("Calculation done.")
