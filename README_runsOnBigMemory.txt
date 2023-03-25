1. connect to a campus-bigmem machine:
srun -A ACF-UTK0011 --partition=campus-bigmem  --qos=campus-bigmem  --pty /bin/bash

2. Compile the code
make ARCH=KL_config_isaac

3. 
exit

4. Run the following to generate the script files
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 700 -en 899 -vst 0 -ven 1214 -np 56

5. edit slscript_big_mem.sh as needed.
For example
the line 
#SBATCH --nodelist=ilm0835
may be changed  based on the particular available nodes or changed to 
#SBATCH --exclude=il[1230-1235]

6. run the batch file
sbatch slscript_big_mem.sh


Please make any changes as needed (for example, machine specific comiles from items 1, 3) and changes needed in 5 for the slscript file

Also, the comman on line 4 may change depending on what we want to run


