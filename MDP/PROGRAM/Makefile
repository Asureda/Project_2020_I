##Makefile per una simulacio de Dinamica Molecular


# Compilador
comp=mpifort
compf=mpif90

# Main program
TARGET=main

# Optimitzador
OPT=-O3 -g 
# Flags-traceback -check all -fp-stack-check

FLAGS=-Wall  -fbounds-check

# Grafics d'energia, pressio, temperatura i funcio de distribució radial
#energy.png : termodynamics_reduced.dat termodynamics_real.dat distrib_func.dat
#    @echo "Generants plots amb els resultats..."
#    gnuplot red.gnu
#    gnuplot real.gnu
#    gnuplot gr.gnu

# Execucio del main program
#termodynamics_reduced.dat : $(TARGET) config.dat constants.dat parameters.dat
#    @echo "Executant el programa amb els valors de l'input..." 
#    ./$(TARGET)

# Compilacio del main program
#$(TARGET) : read_data.f90 parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90 main_vars.f90 verlet.f90 $(TARGET).f90 
#	@echo "From MAKEFILE. Compiling main program"
#	$(comp)  $(TARGET).f90 read_data.f90 main_vars.f90 parallelization.f90 energy.f90 PBC.f90 verlet.f90 vel_initial.f90 Full_Modul_Andersen.f90 Full_modul_Inicialitzar.f90 Full_Modul_Reescalar.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 -o $(TARGET)

$(TARGET) : read_data.o parallelization.o PBC.o Full_Modul_Reescalar.o Full_modul_Inicialitzar.o Full_Modul_Andersen.o Full_Modul_Distrib_Radial.o L_J.o vel_initial.o energy.o main_vars.o verlet.o $(TARGET).o
	@echo "From MAKEFILE. Compiling main program"
	$(comp) $(TARGET).o read_data.o main_vars.o parallelization.o energy.o PBC.o verlet.o vel_initial.o Full_Modul_Andersen.o Full_modul_Inicialitzar.o Full_Modul_Reescalar.o Full_Modul_Distrib_Radial.o L_J.o -o $(TARGET)

# Tots els fitxers amb extensio .f90 son compilats a objectes .o
#%.o : read_data.f90 parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90 main_vars.f90 verlet.f90
#	@echo "From MAKEFILE. Compiling all modules"
#	$(compf) -c read_data.f90 parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90 main_vars.f90 verlet.f90

$(TARGET).o : $(TARGET).f90 read_data.f90 parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90 main_vars.f90 verlet.f90
	$(comp) $(OPT) -c $(TARGET).f90
read_data.o : read_data.f90 main_vars.f90 parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c read_data.f90
main_vars.o :  main_vars.f90 parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c main_vars.f90
parallelization.o : parallelization.f90 PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c parallelization.f90
PBC.o : PBC.f90 Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c PBC.f90
Full_Modul_Reescalar.o : Full_Modul_Reescalar.f90 Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c Full_Modul_Reescalar.f90
Full_modul_Inicialitzar.o : Full_modul_Inicialitzar.f90 Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c Full_modul_Inicialitzar.f90
Full_Modul_Andersen.o : Full_Modul_Andersen.f90 Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c Full_Modul_Andersen.f90
Full_Modul_Distrib_Radial.o : Full_Modul_Distrib_Radial.f90 L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c Full_Modul_Distrib_Radial.f90
L_J.o : L_J.f90 vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c L_J.f90
vel_initial.o : vel_initial.f90 energy.f90  verlet.f90
	$(comp) $(OPT) -c vel_initial.f90
energy.o : energy.f90  verlet.f90
	$(comp) $(OPT) -c energy.f90
verlet.o : verlet.f90
	$(comp) $(OPT) -c verlet.f90




## help: instruccions sobre l'us d'aquest Makefile
.PHONY : help
help :
	@sed -n 's/^##//p' Makefile


## backup : fer una copia comprimida del codi basic
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" *.f90

## clean_all : regla per fer un clean dels objectes executables, resultats i imatges
.PHONY : clean_all
clean_all:
	@echo Removing compiled objects and results
	rm *.o *.dat *.eps *.xyz *.x

## clean_exe : regla per fer un clean nomes dels objectes executables
.PHONY : clean_exe
clean_exe:
	@echo Removing compiled objects only
	rm *.o *.mod