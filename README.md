# Molecular Dynamics Simulations

* Motivation

Master's computational project to learn the basic knowledege of parallel computing applied to molecular dynamics. We develop programs where the tasks of each loop are distributed in different processors. Also, we check the optimization characteristics for each one. The execution was made at the BSC Mare Nostrum where large number of CPUs where used and we worked almost with 400 processors. Another objective of this project is analysing the speed-up and comparing time execution with different number of particles in sequential and parallel program.

* The system

We work with a Van der Waals gas, concretly helium gas, with pair Lenard-Johnes interactions. We consider a system of N particles in a canonical ensemble (NVT ensemble). It lets us change the number of particles, density and diferent steps used during the simulation. Firstly, the system has a FCC structure and it is under periodic bounday conditions, but after a melting process, it becomes a fluid. In this point, the system starts to be studied. The kinetics and thermodynamic parameters are analysed.

* Molecular dynamics

We use the Velocity Verlet algorithm to integrate the equations and set if the user neeeds, switch on the Andersen thermostat. During the simulation with this algorithms, we recalculate the positions and the velocities of the system so many times. If we are working with a heat bath, the expected results are the constant temperature during the simulation as same as the total energy, while the kinetic and potential energy are going to fluctuate.

## First steps 💡
Information to install and execute the programs

### Pre-requisites 📋

Working environment:

```
Linux Shell and Bash
```

Sequential compilers:

```
ifort (Default)
gfortran (Must configure Makefile)
```

Parallel compilers:
```
intel openmpi (Default)
```

### Instalation 🔧

The programs are ready-to-use. The user have to download the repository in a local computer folder or computing cluster, configure the compiler and flags options in the Makefile

Sequential program
```
Makefile:  configure compiler and flags variables (ifort by default)

```
Paralel program (compututing cluster)
```
Makefile:  configure the compiler and flags variables (mpifort by default)
"run_sub.sh" (1): Check the execution order ( mpirun by default)
"run_sub.sh" (2): Configure the submit options ( BSC by default)
"run.sh" (1): Check Makefile flags for ifort or gfortran.
"run.sh" (2): Configure "run.sh" number of processors.

```

## Execution 🚀

Sequential program
```
(1) Configure the simulation parameters (INPUT folder)
(2) Execute the "run.sh" script.
(3) Collect results in the OUTPUT folder.
    The results folder name is the date-time when the task was submitted.

```
Paralel program (compututing cluster)
```
(1) Configure the symulation parameters (INPUT folder)
(2) Execute the "run_sub.sh" script.
(3) Collect results in the OUTPUT folder.
    The results folder name is the date-time when the task was executed.
```
### Program-check 🔎

In the OUTPUT foler is provided a run_check subfolder with input configuration parameters and graphs. 
Put the same parameters in the INPUT files, run the program and compare the graphs. 
They should be similar except for a random factor

### Main theoretical characteristics ⌨️


```
- Inital FCC structure in a cubic volume.
- Uniform random initial velocities.
- Melting and equilibration at a customizable temperature.
- Velocity Verlet algorithm to integrate the equations.
- Andersen Thermostat to control the bath temperature.
- Pair interactions with Lennard-Johnes potential.
- Periodic boundary conditions.
- Thermodinamic results in real and reduced units.
```

## Technologies 🛠️

```
- Fortran
- Open MPI subrouines
- Random numbers: CALL RANDOM_NUMBER(x) (no explicit seed)
- GNU Plot
- Bash shell scripts
- Computing Cluster
```

## Version 📌

Outcome : 21 / 04 / 2020 (version 1.0)

Last moifyed:  NONE (version --)

## Authors ✒️

* **Alexandre Sureda**
* **Elena Ricart**
* **Laia Navarro**
* **Oriol Cabanas**
* **Silvia Àlvarez**


## Acknowledgments 🎁

* Comenta a otros sobre este proyecto 📢
* Invita una cerveza 🍺 o un café ☕ a alguien del equipo. 
* Da las gracias públicamente 🤓.
* etc.

# Appendix
* Input parameters
* Speedup and runnung time recomendations
## A1: Input parameters
parametters.dat
```
particles        # Number of particles (x^3 *4 ; with x natural and positive)
density          # Density (reduced units)
time             # Symulation time (reduced units)
h                # Time step (reduced units)
sigma            # Sigma of the gas (Angstroms)
epsilon          # Epsilon of the gas (kJ/mol)
mass             # Mass (g/mol)
(boolean)        # To add a thermostat
temperature      # If true, temperature of the thermostat (reduced units)
dx               # Precision for the radial distribution function (reduced units)
```
config.dat
```
temperature      # Temperature of the initial melting (reduced units)
iterations       # Melting Velo Verlet Intagration steps
(boolean)        # Print thermodynamic magnitudes
iterations       # Delta iterations to measure thermodynamic magnitudes
(boolean)        # Compute the radial distribution function
iterations       # Delta iterations to compute the Rad. Dist. Func
(boolean)        # Time-positions of the particles (.xyz file)
iterations       # Delta iterations to save the positions
```
constants.dat
```
0.008314462      # Boltzman constant in kJ/mol.K
6.022d23         # Avogadro number
```
## A2: Speedup and runnung time recomendations
MDP-Double Work
```
- Same numer of particles for each processator
- For each particle is computed the interactions with the others. 
     No symetric reduction is maid to comupute the half of the matrix.
- Every processator have the same work
```
MDP- Pair
```
- Same number of interactions for each patricle
- It is computed the direct and the symetric term. We loop over half of the matrix
- Every processator have the same work
```
MDP- Symetric Matrix
```
- Same numer of particles for each processator
- It is computed the direct and the symetric term. We loop over half of the matrix
- Different distribution of the work
```

