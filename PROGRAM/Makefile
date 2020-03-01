##Makefile per una simulacio de Dinamica Molecular


# Compilador
F90=gfortran

# Main program
TARGET=main

# Optimitzador
OPT=-O

# Flags
FLAGS=-Wall  -fbounds-check

# Grafics d'energia pressio, temperatura i moment total
energy.eps : Results.txt
    @echo "Generants plots amb els resultats..."
    gnuplot Scripts_GNUPlot/plot_Energy_Raquel.gnu
    gnuplot Scripts_GNUPlot/plot_momentum.gnu
    gnuplot Scripts_GNUPlot/plot_rdf.gnu
    @echo Copying to Results folder...
    mkdir -p Results
    cp *.txt *.dat *.eps *.xyz Results
    @echo "Fet!"

# Execucio del main program
Results.txt : $(TARGET).x INPUT/config.dat constants.dat parameters.dat
    @echo "Executant el programa amb els valors de l'input..." 
    ./$(TARGET).x < INPUT/config.dat constants.dat parameters.dat

# Compilacio del main program
$(TARGET).x : $(TARGET).o energy.o PBC.o verlet.o vel_initial.o Full_modul_Inicialitzar.o Full_Modul_Reescalar.o Full_Modul_Distrib_Radial.o L_J.o read_data.o main_vars.o
    @echo "Compilant el programa main.x..."
    $(F90) $(OPT) $(FLAGS) -o $(TARGET).x $(TARGET).o energy.o PBC.o verlet.o vel_initial.o Full_modul_Inicialitzar.o Full_Modul_Reescalar.o Full_Modul_Distrib_Radial.o L_J.o read_data.o main_vars.o

# Tots els fitxers amb extensio .f90 son compilats a objectes .o
%.o : %.f90
    @echo "Compilant les subrutines necessaries..."
    $(F90) $(OPT) $(FLAGS) -c $<


# ## statistics : binning of the time series for different magnitudes
# .PHONY : statistics
# statistics :
#     $(F90) -o binning.x Code/binning2.f90
#     ./binning.x
#     gnuplot Scripts_GNUPlot/plot_binning.gnu
#     @echo Moving results to Results_binning
#     mkdir -p Results_binning
#     cp *binning.txt *mit.txt *bin.eps Results_binning

## help: instruccions sobre l'us d'aquest Makefile
.PHONY : help
help :
    @sed -n 's/^##//p' Makefile


## backup : fer una copia comprimida del codi basic
.PHONY : backup
backup:
    tar -czvf "backup.tar.gz" PROGRAM/*.f90

## clean_all : regla per fer un clean dels objectes executables, resultats i imatges
.PHONY : clean_all
clean_all:
    @echo Removing compiled objects and results
    rm *.o *.dat *.eps *.xyz *.x

## clean_exe : regla per fer un clean nomes dels objectes executables
.PHONY : clean_exe
clean_exe:
    @echo Removing compiled objects only
    rm *.o *.x