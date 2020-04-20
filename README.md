# Molecular Dynamics Symulations

Lennard-Jones potential with parir interactions.

Computational progect to develop sequenctial and paralelized programs to symulate molecular dynamics. Ready to compute the energies, temperature, pressure and radial distribution function in real and reduced units. For an N patricles system (compatible with FCC structure) the integration is done with the Velocity Verlet algorigthm. The program also come with an optional Andersen thermostat.
Also, pair interactions with Lennard-Johnes potential and periodic boundary conditions are included

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

The programs are ready-to-use. The user have to download the repository in a local computer folder or computing cluster, configure the compiller and flags options in the Makefile

Sequential program
```
Makefile:  configure compiller and flags variables (ifort by default)

```
Paralel program (compututing cluster)
```
Makefile:  configure the compiler anf flags variables (mpifort by default)
"run_sub.sh" (1): Check the execution order ( mpirun by default)
"run_sub.sh" (2): Configure the submit options ( BSC by default)
```

## Execution 🚀

Sequential program
```
(1) Configure the symulation parameters (INPUT folder)
(2) Execute the "run.sh" script.
(3) Collect the results in the OUTPUT folder
    The results folder name is the date-time whem the task was executed

```
Paralel program (compututing cluster)
```
(1) Configure the symulation parameters (INPUT folder)
(2) Execute the "run_sub.sh" script.
(3) Collect the results in the OUTPUT folder
    The results folder name is the date-time whem the task was executed
```
### Program-check 🔎

In the OUTPUT foler is provided a run_check subfolder with input configuration parameters and graphs. 
Put the same paremeters in the INPUT files, run the program and compare the graphs. 
They should be similar except for a random factor

### Main theoretical characterisitics ⌨️


```
- Inital FCC structure
- Melting and equilibration
- Velocity Verlet
- Andersen Thermostat
- Pair interaction with Lenard-Johnes
- Periodic boundary conditions

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

## Licence 📄

Blblablabla

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
- 
-
```
config.dat
```
- 
-
```
constants.dat
```
- 
-
```
## A2: Speedup and runnung time recomendations
MDP-Double Work
```
- 
-
```
MDP- Pair
```
- 
-
```
MDP- Symetric Matrix
```
- 
-
```
MDP- Distributed Work
```
- 
-
```

