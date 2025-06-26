# ThermPropDebye
Scripts to calculate thermal properties of acoustic phonons from the elastic constant tensor.

The main python program must be run in a directory which contains a single input file that has the suffix ".clmz".
Each letter in the suffix represents a field of data that needs to be in the input file. 'c' is the 6x6 matrix of
elastic constants which is to be given in GPa, 'l' is the set of three lattice vectors which are to be given in
angstroms, 'm' is the mass per unit cell which is to be given in atomic mass units, and 'z' is the number of molecules
per unit cell. An example input file named "test.clmz" is included.

The program will write output to a file with the suffix ".
