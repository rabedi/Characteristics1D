Refer to lines 616-621 of file rdb.py
# and change these path as needed before running the program:

    #README:
    # 1. set the path to "InhomogeneousFiles" folder
    InpF.setInputMeshRootFolder("../InhomogeneousFiles")
    # 2. set path to "2023_03_24". I put it inside a "data" folder for myself. You can adjust this path
    folderSource = "../../data/2023_03_24"
    # 3. If needed, adjust output path were the files are generated
    folderDest = "../../data/Characteristics_data"



--------------------------------
Data Transfer

If you just want to work with the averages of responses for about 2300 combination of scalar values inp_s_ + (la, dd2, lc, and ldelc) you can directly get this file:
https://rezaabedi.com/share/HTMDEC/stat_mean.csv

If you want to also for individual runs with their corresponding random field realizations, we have two options:
A) Get the raw data from me and generate the final files under B. I prefer this one as the spatial fields csv file that is under B is about 15 GB (unzipped) whereas the "InhomogeneousFiles" folder is only 3.7 GB and you only need to download it once, with the catch that everytime I send you raw data, you need to run the python file rdb.py for about a couple of hours.

Files for A:
A1. Inhomogeneous "base" files (about 8-9k of meshes with 2^14+1 points). Only needs to be downloaded once.
https://rezaabedi.com/share/HTMDEC/InhomogeneousFiles.zip

A2. Runs raw data: about 52k runs individual data stored in multiple csv files
https://rezaabedi.com/share/HTMDEC/2023_03_24.zip

A3. Python code to turn A1 and A2 to B1, B2 (all run files)
https://rezaabedi.com/share/HTMDEC/python.zip

------------------------------------------------------
B) 
B1 (scalar in/out) and B2 (field in) are only needed if individual runs are to be used. 
B1. The scalar inputs and outputs for all runs (about 52k now):
https://rezaabedi.com/share/HTMDEC/all_wxs.zip

B2. fields
https://rezaabedi.com/share/HTMDEC/all_fx.zip


