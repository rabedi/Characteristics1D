configName = 'test.txt';
configName = 'configMaker_fine_axt_gen_Tri.txt';
%configName = 'configMaker_axt_resolution_x_Fracture_New.txt';
versionOffset = 0;
configName = 'configMaker_fine_axt_gen_Tri_lam2.txt';
configName = 'configMaker_fine_axt_gen_Tri_lam2p5.txt';
configName = 'configMaker_fine_axt_gen_Tri_lam3.txt';
if (strcmp(configName, 'configMaker_fine_axt_gen_Tri_lam2.txt') == 1)
    versionOffset = 1620;
end
if (strcmp(configName, 'configMaker_fine_axt_gen_Tri_lam2p5.txt') == 1)
    versionOffset = 1755;
end
if (strcmp(configName, 'configMaker_fine_axt_gen_Tri_lam3.txt') == 1)
    versionOffset = 1890;
end

i1df = Instruction_1DFragmentation;
fid = fopen(configName, 'r');
success = 1;
while (success)
    i1df = i1df.read(fid, versionOffset);
    success = i1df.success;
end
fclose(fid);