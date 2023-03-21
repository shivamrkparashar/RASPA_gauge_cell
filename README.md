Gauge Cell Simulations in RASPA
======

RASPA is a general purpose classical simulation package that can be used for the
simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields. 
We implemented the Gauge Cell Monte Carlo simulation method in RASPA. This method 
is used to generate the full canonical adsorption isotherm in the form of a Van der Walls
loop.


Installation
============

## using Cmake
```
cd src_gauge_cell
cmake .
make
```

Usage
=====

## To run Gauge Cell simulations

In simulation.input file, make the following changes-

1. Put `GaugCellSize  1000`  # Box lenght of gauge cell (Angstroms)
2. Put `ComputeNumberOfMoleculesHistogram yes` # Molecule histogram is required to calculate several points on isotherm from one simulation
3. Put `GaugeCellSwapProbability    1.0`  # Gauge cell exchange move for adsorbates 
4. Under Component 0, put `CreateNumberOfMolecules 100`  # This is Ntotal = Ngauge + Npore which will remain constant throughout the simulation

Output
=======

1. Each simulation generates `HistogramGaugeCell.dat` file in the simulation directory. 
This file contain histogram of particles in the gauge cell and pore cell along with the chemical potential.
These points can be directly plotted and are referred as the ideal gas gauge cell (IGGC) points.

2. In the Output/System_0/output\*data file, the `Average Ideal Gauge Cell chemical potential` is printed. This chemical potential 
is based on the average density of particles in the pore cell during the simulaton. This point is known as the mean density gauge cell (MDGC) point.


References for Gauge Cell Simulations
=====================================
1. A. V. Neimark and A. Vishnyakov. A Simulation Method for the Calculation of Chemical Potentials in Small, Inhomogeneous, and Dense Systems, Journal of Chemical Physics, 2005, V. 122, 234108.
2. A. V. Neimark and A. Vishnyakov. Phase transitions and criticality in small systems: vapor-liquid transition in nanoscale spherical cavities, Journal of Physical Chemistry, 2006, V.110, p.p.9403-9412.
3. A. V. Neimark and A. Vishnyakov. The Gauge Cell Method for Molecular Simulation Studies of Phase Transitions in Confined Fluids, Phys. Rev. E, 2000, V.62, p.p.4611-4622.
4. A. Vishnyakov and A. V. Neimark. Studies of Liquid-Vapor Equilibrium, Criticality and Spinodal Transitions in Nanopores by the Gauge Cell Monte Carlo Simulation Method, Journal of Physical Chemistry B, 2001, V. 105, p.p. 7009-7020.
