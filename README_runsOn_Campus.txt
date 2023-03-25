I assume in addition to 1 long bigmem, you can run on 4 (or 5?) campus nodes. 
Please follow the same process for 4 to 5 nodes

1. connect to a campus node:
srun -A ACF-UTK0011 --partition=campus  --qos=campus  --pty /bin/bash

2. Compile the code
make ARCH=KL_config_isaac

3. 
exit

4. Run the following to generate the script files
For node 1 to 5 that you'll have access to, please in order run the following
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 900 -en 949 -vst 0 -ven 1214 -np 11
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 950 -en 999 -vst 0 -ven 1214 -np 11
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 1000 -en 1049 -vst 0 -ven 1214 -np 11
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 1050 -en 1099 -vst 0 -ven 1214 -np 11
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 1100 -en 1149 -vst 0 -ven 1214 -np 11
./solver -sp -mc config/MainConfig/mainConfig_axt_coh.txt -st 1150 -en 1199 -vst 0 -ven 1214 -np 11

....
and continue for subsequent dates (we can discuss this later after you do it for day 1)

5. edit slscript_campus.sh as needed.
similar to bigmem is node may be specified or nodes excluded ...

6. run the batch file
sbatch slscript_campus.sh


Please make any changes as needed (for example, machine specific comiles from items 1, 3) and changes needed in 5 for the slscript file

Also, the comman on line 4 may change depending on what we want to run


