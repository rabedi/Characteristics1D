configName = 'test.txt';
configName = 'configMaker_fine_axt_gen_Tri.txt';
%configName = 'configMaker_axt_resolution_x_Fracture_New.txt';
versionOffset = 0;

i1df = Instruction_1DFragmentation;
fid = fopen(configName, 'r');
success = 1;
while (success)
    i1df = i1df.read(fid, versionOffset);
    success = i1df.success;
end
fclose(fid);