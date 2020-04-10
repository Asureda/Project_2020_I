#Execute file
#Getting the temporal file name
tmp_folder=$(date +'%Y_%m_%d_%H_%M_%S')

# Creating the tmp file in PROGRM

mkdir PROGRAM/$tmp_folder/

# Copy input parametters to tmp file
rm PROGRAM/*.dat
cp INPUT/*.dat PROGRAM/$tmp_folder/

# Enter in PROGRAM folder

cd PROGRAM/

# Execute Makefile to get the last version

make

# Copy the programs and scripts to the tmp folder
# The thing we need start with "r_*"

cp main $tmp_folder/
cp *.gnu $tmp_folder/

#Enter in tmp folder

cd $tmp_folder/

#Execute the progrm

./main
gnuplot gr.gnu
gnuplot real.gnu
gnuplot red.gnu

# After the program finalizes, we go back

#rm r_*

cd ..
cd ..
mkdir OUTPUT/$tmp_folder/
# We crate the output folder in OUTPUT and copy results

cp PROGRAM/$tmp_folder/* OUTPUT/$tmp_folder/

#Delete the tmp folder from PROGRAM

 cd PROGRAM/
 rm -r $tmp_folder/
# cd ..


# END     ##################################################


# echo "Deleting Data files from PROGRAM"
# rm PROGRAM/*.dat
# echo "Copy files from INPUT to PROGRAM"
# cp INPUT/*.dat PROGRAM/
# echo "Inside /PROGRAM"
# cd PROGRAM/
# echo "Executing MakeFile"
# make
# echo "Make done"
# #. Deleting object files"
# #make clean_exe
# echo "Executing program"
# ./main
# echo "From bash. Main program ended"
# echo "Executing Plot scripts"
# gnuplot red.gnu
# gnuplot real.gnu
# gnuplot gr.gnu
# echo "Getting out from PROGRAM/"
# cd ..
# echo "Moving files from PROGRAM to OUTPUT"
# mv PROGRAM/thermodynamics_reduced.dat OUTPUT/thermodynamics_reduced.dat
# mv PROGRAM/thermodynamics_real.dat OUTPUT/thermodynamics_real.dat
# mv PROGRAM/distrib_funct.dat OUTPUT/distrib_funct.dat
# mv PROGRAM/positions.xyz OUTPUT/positions.xyz
# echo "Moving graphs"
# mv PROGRAM/*.png OUTPUT/
# echo "Deleting all data files from PROGRAM"
# rm PROGRAM/*.dat
# echo "From bash. End of script"

