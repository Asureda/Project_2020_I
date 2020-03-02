#Execute file
echo "Deleting Data files from PROGRAM"
rm /PROGRAM/*.dat
echo "Copy files from INPUT to PROGRAM"
cp /INPUT/*.dat /PROGRAM/
echo "Inside /PROGRAM"
cd PROGRAM/
echo "Executing MakeFile"
make
echo "Make done. Deleting object files"
make clean_exe
echo "Executing program"
./main
echo "From bash. Main program ended"
echo "Getting out from PROGRAM/"
cd ..
echo "Moving files from PROGRAM to OUTPUT"
mv /PROGRAM/thermodynamics_reduced.dat /OUTPUT/thermodynamics_reduced.dat
mv /PROGRAM/thermodynamics_real.dat /OUTPUT/thermodynamics_real.dat
mv /PROGRAM/distriv_funct.dat /OUTPUT/distriv_funct.dat
mv /PROGRAM/positions.xyz /OUTPUT/positions.xyz
echo "Deleting all data files from PROGRAM"
rm /PROGRAM/*.dat
echo "From bash. End of script"

