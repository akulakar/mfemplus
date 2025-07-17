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
