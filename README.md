# ThermPropDebye
Scripts to calculate thermal properties of acoustic phonons from the elastic constant tensor using the Debye
approximation. Intended primarily for use with molecular crystals.

In order to use this program, you first need to know the elastic constant tensor for your crystal. It doesn't
matter how you acquire this but I will suggest using a pair of octave scripts written by Dr. Alberto Otero-de-
la-Roza which generate strained unit cells and, after performing DFT calculations to get the resulting stress,
generates the elastic constant tensor. These scripts can be found here:
https://github.com/aoterodelaroza/tricks/tree/master/elastic

The main python program "ThermPropDebye.py" must be run in a directory which contains a single input file that has
the suffix ".clmz". Each letter in the suffix represents a field of data that needs to be in the input file. 'c' is
the 6x6 matrix of elastic constants which is to be given in GPa, 'l' is the set of three lattice vectors which are to
be given in angstroms, 'm' is the mass per unit cell which is to be given in atomic mass units, and 'z' is the number
of molecules per unit cell. Each of these fields is separated by a single blank line. An example input file named
"test.clmz" is included.

The program will write output to a file with the suffix ".therm". The first column is the temperature in K, the
second column is the free-energy in kJ/mol per molecule. The third column is the entropy in kJ/mol/K per molecule,
and the fourth column is the energy in kJ/mol per molecule.

An explanation of the relevant theory can be found in the following literature:
NOT YET AVAILABLE
