
#include "globalMacros.h"
#include "TestSuites.h"
#include "SL_OneInterfaceAllTimes.h"
#include "InhomogeousField.h"
#include "InhomogeneousElasticFracture.h"
#include "Domain_AllInterfacesAllTimes.h"
#include "SimpleConfigMaker.h"
#include "globalFunctions.h"

#include "DomainPostProcessingS2.h"
#include "DomainPostProcessingS3.h"
#include "TSR1D.h"


// returns true if this particular serial number was ran
bool Solve_1_serialNumber(SolveParameters& solvePara, unsigned int serNumIn, int versionNumIn);
void Solve_all_serialNumbers(SolveParameters& solvePara);

int main(int argc, char *argv[])
{
//	Test_TSR_Ortiz_1D();
//	return 0;

	setGlobalMembers();
	SolveParameters solvePara;

	if (argc > 1)
	{
		string opt = (string)argv[1];
		name2Type(opt, solvePara.sOpt);

		for (int i = 2; i < argc; ++i)
		{
			if (strcmp(argv[i], "-mc") == 0)
			{
				string solveParametersConfigName;
				solveParametersConfigName = (string)argv[++i];
				solvePara.Read_SolveParameters(solveParametersConfigName);
			}
			else if (strcmp(argv[i], "-c") == 0)
			{
				solvePara.configName = (string)argv[++i];
			}
			else if ((strcmp(argv[i], "-cp") == 0) || (strcmp(argv[i], "-pc") == 0))
			{
				solvePara.configPPName = (string)argv[++i];
			}
			else if (strcmp(argv[i], "-st") == 0)
			{
				solvePara.serialNumber_st = (unsigned int)atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-en") == 0)
			{
				solvePara.serialNumber_en = (unsigned int)atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-vst") == 0)
			{
				solvePara.versionNumber_st = (unsigned int)atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-ven") == 0)
			{
				solvePara.versionNumber_en = (unsigned int)atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-cv") == 0)
			{
				solvePara.version_configMakerGenName = (string)argv[++i];
			}
			else if (strcmp(argv[i], "-np") == 0)
			{
				solvePara.numParallelRuns = (unsigned int)atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-ls") == 0)
			{
				solvePara.low_disk_space = true;
			}
			else if (strcmp(argv[i], "-cmg") == 0)
			{
				string configMakerInstructionName = "config/_configMakerInstructionNNAR3.txt";
				i++;
				if (i < argc)
					configMakerInstructionName = argv[i];
				bool forceRewrite = true;
				i++;
				if (i < argc)
					forceRewrite = (bool)argv[i];
				sfcm_gen.CreateSimpleFormatConfigMakerFromInstructions(configMakerInstructionName, forceRewrite);
				EXIT;
			}
			// reading config makers simple version
			else if (strcmp(argv[i], "-cmr") == 0)
			{
				i++;
				int cntr = atoi(argv[i]);
				string configMaker = "config/_configMakerNNAR3.txt";
				i++;
				if (i < argc)
					configMaker = argv[i];
				bool success = sfcm.Read(configMaker, cntr);
				if (!success)
				{
					cout << "for configMaker\t" << configMaker << "\tcntr\t" << cntr << "\tdoes not exist => Exiting the code\n";
					EXIT
				}
			}
			else
			{
				cout << "i\t" << i << '\n';
				cout << "argv[i]\t" << argv[i] << '\n';
				THROW("ERROR: Invalid flag--value option");
			}
		}
		if (opt == "wn")
		{
			int i = 2;
			int numVertices = 16385; //1025, 
			int numRealizations = 5000;
			i++;
			if (i < argc)
				numVertices = (unsigned int)atoi(argv[i]);
			i++;
			if (i < argc)
				numRealizations = (unsigned int)atoi(argv[i]);
			GenerateWhiteRandomFile(numVertices, numRealizations);
			return 0;
		}
	}

	Solve_all_serialNumbers(solvePara);
#if VCPP
	cout << "press Any key\n";
	getchar();
#endif
	if (g_slf_conf != NULL)
		delete g_slf_conf;
	return 0;
}

bool Solve_1_serialNumber(SolveParameters& solvePara, unsigned int serNumIn, int versionNumIn)
{
	g_versionNumber = versionNumIn;
	if (g_versionNumber >= 0)
	{
		g_versionNumber_wOffset = g_versionNumber + solvePara.versionOffset;
		toString(g_versionNumber_wOffset, g_versionNumber_str);
		g_versionNumber_str = "_V_" + g_versionNumber_str;
		sfcm_gen.CreateSimpleFormatConfigMakerFromInstructions(solvePara.version_configMakerGenName, solvePara.version_configMaker_forceRewrite);
		bool success = sfcm.Read(sfcm_gen.fileNameOut, g_versionNumber);
		if (success == false)
			return true;
		Configure_sfcm_sfcm_gen();
	}
	else
		g_versionNumber_str = "";

	g_serialNumber = serNumIn;
	string str_serial;
	toString(g_serialNumber, str_serial);

	g_prefileName = "run" + g_versionNumber_str + "_" + str_serial;
	if (solvePara.isDomain)
	{
		vector<string> parts = getFileParts(solvePara.configName);
		string cn = parts[1];

		string lg = "_log";
		fileOperation(makeD, lg);

		string fileName_log = lg + "/" + cn + "_" + g_prefileName + ".txt";
		fstream in(fileName_log.c_str(), ios::in);
		if (in.is_open())
			return false;

		cout << "=======================================================\n";
		cout << "============ version " << versionNumIn << " serial number " << serNumIn << " ==========\n";
		cout << "=======================================================\n";
		if (g_prefileName != "")
		{
			fileOperation(makeD, g_prefileName);
			g_prefileNamePPS2 = "_PPS2_";
			g_prefileNamePPS2 += g_prefileName;
			fileOperation(makeD, g_prefileNamePPS2);
			g_prefileName += "/";
			g_prefileNamePPS2 += "/";
		}
		if (solvePara.b_solveParametersConfigName == true)
			CopyFile2OutputDirectory(solvePara.solveParametersConfigName);
		if (versionNumIn >= 0)
			CopyFile2OutputDirectory(solvePara.version_configMakerGenName);

		if ((solvePara.sOpt == so_domain_sp) || (solvePara.sOpt == so_domain_s))
		{
			int success = MAIN_Domain(solvePara.configName, g_serialNumber);
			cout << "Domain: solve: success:\t" << success << '\n';
		}
		if ((solvePara.sOpt == so_domain_sp) || (solvePara.sOpt == so_domain_p))
		{
			MAIN_DomainPostProcessS3(solvePara.configPPName);
			cout << "Domain: PPS3: success:\t" << 1 << '\n';
		}
		fstream out(fileName_log.c_str(), ios::out);
		out << 1 << '\n';
		return true;
	}
	if (solvePara.sOpt == so_interface_s)
	{
		////////////////////////////////////////////////////////////////////////////////
		// 1-space, all-time points
		MAIN_SL_OneInterfaceAllTimes_ONE_Interface(solvePara.configName);
		return true;
	}
	if (solvePara.sOpt == so_onePoint_s)
	{
		////////////////////////////////////////////////////////////////////////////////
		// 1-space, 1-time point
		Test_OnePointRiemannSolutions(solvePara.configName);
		return true;
	}
	if (solvePara.sOpt == so_elfrac_fields)
	{
		TestInhomogeneousElasticFractorField(solvePara.configName, &g_serialNumber, &solvePara.isPeriodic, &solvePara.xM, &solvePara.xm);
		cout << "success\t1\n";
		return true;
	}
	if (solvePara.sOpt == so_one_field)
	{
		TestInhomogeneousField(solvePara.configName);
		cout << "success\t1\n";
		return true;
	}
	if (solvePara.sOpt == so_domain_p2)
	{
		MAIN_DomainPostProcessS2(solvePara.configName);
		cout << "success\t1\n";
		return true;
	}
	cout << "solvePara.sOpt\t" << solvePara.sOpt << '\n';
	THROW("Invlaid option\n");
}

void Solve_all_serialNumbers(SolveParameters& solvePara)
{
	g_low_disk_space = solvePara.low_disk_space;
	solvePara.InitializeAfterSetup();
	if (solvePara.numParallelRuns <= 0)
	{
		if (solvePara.lv1s2)
		{
			for (int versionNumber = solvePara.versionNumber_st; versionNumber <= solvePara.versionNumber_en; ++versionNumber)
			{
				for (int serNum = solvePara.serialNumber_st; serNum <= solvePara.serialNumber_en; ++serNum)
				{
					if (g_slf_conf != NULL)
						delete g_slf_conf;
					g_slf_conf = new SLFractureGlobal_Configuration;
					Solve_1_serialNumber(solvePara, serNum, versionNumber);
				}
			}
		}
		else
		{
			for (int serNum = solvePara.serialNumber_st; serNum <= solvePara.serialNumber_en; ++serNum)
			{
				for (int versionNumber = solvePara.versionNumber_st; versionNumber <= solvePara.versionNumber_en; ++versionNumber)
				{
					if (g_slf_conf != NULL)
						delete g_slf_conf;
					g_slf_conf = new SLFractureGlobal_Configuration;
					Solve_1_serialNumber(solvePara, serNum, versionNumber);
				}
			}
		}
	}
	else
	{
		vector<fstream*> outPtrs(solvePara.numParallelRuns);
		for (int pi = 0; pi < solvePara.numParallelRuns; ++pi)
		{
			string str;
			toString(pi, str);
			str = "script" + str + ".sh";
			outPtrs[pi] = new fstream();
			outPtrs[pi]->open(str.c_str(), ios::out);
		}
		unsigned int cntr = 0, remainder;
		if (solvePara.lv1s2)
		{
			for (int versionNumber = solvePara.versionNumber_st; versionNumber <= solvePara.versionNumber_en; ++versionNumber)
			{
				for (int serNum = solvePara.serialNumber_st; serNum <= solvePara.serialNumber_en; ++serNum)
				{
					remainder = cntr++ % solvePara.numParallelRuns;
					(*(outPtrs[remainder])) << "./solver -sp - mc " << solvePara.solveParametersConfigName << " -st " << serNum << " -en " << serNum << " -vst " << versionNumber << " -ven " << versionNumber << '\n';
				}
			}
		}
		else
		{
			for (int serNum = solvePara.serialNumber_st; serNum <= solvePara.serialNumber_en; ++serNum)
			{
				for (int versionNumber = solvePara.versionNumber_st; versionNumber <= solvePara.versionNumber_en; ++versionNumber)
				{
					remainder = cntr++ % solvePara.numParallelRuns;
					(*(outPtrs[remainder])) << "./solver -sp - mc " << solvePara.solveParametersConfigName << " -st " << serNum << " -en " << serNum << " -vst " << versionNumber << " -ven " << versionNumber << '\n';
				}
			}
		}
		for (int pi = 0; pi < solvePara.numParallelRuns; ++pi)
			delete outPtrs[pi];
	}
}

