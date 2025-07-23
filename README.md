# Potential changes and notes to be considered:
1. Use free functions for element stress strain rather than a derived class. Use MFEM access functions instead. 
2. OpenMP for the for loops in GlobalStressStrain has not been successful. This is because the pointers are of type abstract base class, and each thread points to the same object. Since it is in an abstract base class, a copy cannot be created either.

# MFEMPLUS

*MFEMPLUS* is a set of cpp files that are extensions to the *MFEM* library.
It contains the following cpp files (with corresponding hpp files):

## 1. customintegrators.cpp

Includes the  following integrators:
- IsotropicElasticity (with elastic constants E and NU)
- AnisotropicElasticity (with stiffness matrix C)

Elasticity integrators use the B matrix for assembly. Refer Hughes' book for reference.

## 2. customfunctions.cpp

Includes the following classes and functions:
- ElementStressStrain
	- ComputeElementStrain
	- ComputeElementStress
	- ComputeBoundaryElementArea
- GlobalStressStrain
	- GlobalStrain
	- GlobalStress
	- ComputeBoundaryForce
