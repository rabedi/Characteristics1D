#include "globalTypesClasses.h"
#include "globalFunctions.h"
#include "commonMacros.h"

string getName(BoundaryConditionT dat)
{
	if (dat == bct_Undecided)
		return "Undecided";
	if (dat == bct_Dirichlet)
		return "Dirichlet";
	if (dat == bct_Neumann)
		return "Neumann";
	if (dat == bct_Characteristics)
		return "Characteristics";
	if (dat == bct_Unspecified)
		return "Unspecified";
	if (dat == bct_Symmetric)
		return "Symmetric";
	if (dat == bct_AntiSymmetric)
		return "AntiSymmetric";
	if (dat == bct_PeriodicOrBloch)
		return "PeriodicOrBloch";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, BoundaryConditionT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= BoundaryConditionT_SIZE)
			THROW("too large of a number\n");
		typeVal = (BoundaryConditionT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = -1; i < BoundaryConditionT_SIZE; ++i)
	{
		typeVal = (BoundaryConditionT)i; // casting integer to BoundaryConditionT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading BoundaryConditionT\n");
}

//operator for output
ostream& operator<<(ostream& out, BoundaryConditionT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, BoundaryConditionT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

AdaptivityS::AdaptivityS()
{
	a_flag = a_unassigned;
	a_delt = 0.0;
}

void AdaptivityS::Update_AdaptFlag(AdaptivityF a_flagIn)
{
	a_flag = (AdaptivityF)MAX((int)a_flagIn, (int)a_flag);
}

ostream & operator<<(ostream & out, const AdaptivityS & dat)
{
	out << (int)dat.a_flag << '\t' << dat.a_delt;
	return out;
}

/////////////////
string getName(errorCheckT dat)
{
	if (dat == ect_notActive)
		return "notActive";
	if (dat == ect_Dimensional)
		return "Dimensional";
	if (dat == ect_nonDimensional)
		return "nonDimensional";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, errorCheckT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= errorCheckT_SIZE)
			THROW("too large of a number\n");
		typeVal = (errorCheckT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < errorCheckT_SIZE; ++i)
	{
		typeVal = (errorCheckT)i; // casting integer to errorCheckT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading errorCheckT\n");
}

//operator for output
ostream& operator<<(ostream& out, errorCheckT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, errorCheckT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getName(IOF_type dat)
{
	if (dat == iof_none)
		return "none";
	if (dat == iof_ascii)
		return "ascii";
	if (dat == iof_binary)
		return "binary";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, IOF_type& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= IOF_type_SIZE)
			THROW("too large of a number\n");
		typeVal = (IOF_type)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < IOF_type_SIZE; ++i)
	{
		typeVal = (IOF_type)i; // casting integer to IOF_type, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading IOF_type\n");
}

//operator for output
ostream& operator<<(ostream& out, IOF_type dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, IOF_type& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

bool getExt(IOF_type dat, string& ext)
{
	if (dat == iof_ascii)
	{
		ext = "txt";
		return true;
	}
	if (dat == iof_binary)
	{ 
		ext = "bin";
		return true;
	}
	if (dat == iof_none)
	{
		ext = "none";
		return false;
	}
		cout << dat;
		THROW("Invalid dat\n");
}

/////////////////
string getName(solveOptions dat)
{
	if (dat == so_domain_sp)
		return "sp";
	if (dat == so_domain_s)
		return "s";
	if (dat == so_domain_p)
		return "p";
	if (dat == so_domain_p2)
		return "p2";
	if (dat == so_interface_s)
		return "i";
	if (dat == so_onePoint_s)
		return "o";
	if (dat == so_elfrac_fields)
		return "eff";
	if (dat == so_one_field)
		return "f";
	if (dat == so_configGen)
		return "cmg";
	if (dat == so_configRead)
		return "cmr";
	if (dat == so_wnRandomFieldGen)
		return "wn";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, solveOptions& typeVal)
{
	unsigned int sz = name.size();
	if ((sz > 0) && (name[0] = '-'))
		name = name.substr(1, sz);

	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= solveOptions_SIZE)
			THROW("too large of a number\n");
		typeVal = (solveOptions)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < solveOptions_SIZE; ++i)
	{
		typeVal = (solveOptions)i; // casting integer to solveOptions, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading solveOptions\n");
}

//operator for output
ostream& operator<<(ostream& out, solveOptions dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, solveOptions& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

void PrintHeader(const vector<HeaderLabels>& hLabels, ostream& out, bool isCommaSeparated)
{
	unsigned int sz = hLabels.size();
	if (sz == 0)
		return;
	out << hLabels[0].group;
	string sep = "\t";
	if (isCommaSeparated)
		sep = "\t,\t";
	for (unsigned int i = 1; i < sz; ++i)
		out << sep << hLabels[i].group;
	out << '\n';

	out << hLabels[0].subgroup;
	for (unsigned int i = 1; i < sz; ++i)
		out << sep << hLabels[i].subgroup;
	out << '\n';

	out << hLabels[0].textLabel;
	for (unsigned int i = 1; i < sz; ++i)
		out << sep << hLabels[i].textLabel;
	out << '\n';

	out << hLabels[0].latexLabel;
	for (unsigned int i = 1; i < sz; ++i)
		out << sep << hLabels[i].latexLabel;
	out << '\n';
}

void GetSubdomainIndexed_TimeIndexed_FileName(string& fileName, int subdomainNo, int timeIndex, const string& specificName, string ext, bool addUnderline)
{
	string posStr = "", nameBase, subDomainStr = "";
	if (timeIndex >= 0)
	{
		toString(timeIndex, posStr);
		posStr = "_tPos_" + posStr;
	}
	string root, sd = "";
	root = g_prefileName;
	if (subdomainNo >= 0)
	{
		toString(subdomainNo, subDomainStr);
		if (!addUnderline)
			sd = "sd_";
		else
			sd = "_sd_";
	}
	//	nameBase = root + sd + subDomainStr + posStr + "_" + specificName;
	nameBase = root + sd + subDomainStr + "_" + specificName + posStr;
	fileName = nameBase + "." + ext;
}

void GetSubdomainIndexed_SpaceIndexed_FileName(string& fileName, int spaceIndexOr_domainIndexLeft, const string& specificName, string ext, int domainIndexRight)
{
	string posStr = "", nameBase, subDomainStr;
	if (spaceIndexOr_domainIndexLeft >= 0)
	{
		if (domainIndexRight < 0)
		{
			toString(spaceIndexOr_domainIndexLeft, posStr);
			posStr = "xPos_" + posStr;
		}
		else
		{
			string posStrL, posStrR;
			toString(spaceIndexOr_domainIndexLeft, posStrL);
			toString(domainIndexRight, posStrR);
			posStr = "xPos_dl_" + posStrL + "_dR_" + posStrR;
		}
	}
	nameBase = g_prefileName + posStr + "_" + specificName;
//	nameBase = g_prefileName + "_" + specificName + posStr;
	fileName = nameBase + "." + ext;
}

SolveParameters::SolveParameters()
{
	sOpt = so_domain_sp;
	configName = "none";
	configPPName = "none";

	serialNumber_st = 0;
	serialNumber_en = -1;

	isPeriodic = false;
	xm = 0.0, xM = 10.0;

	versionNumber_st = -1;
	versionNumber_en = -1;
	version_configMakerGenName = "config/ConfigMakerVersion/configMaker_axt.inst";
	version_configMaker_forceRewrite = true;
	versionOffset = 0;
	lv1s2 = false;
	numParallelRuns = -1;
	solveParametersConfigName = "config/mainConfig_axt_coh.txt";
	b_solveParametersConfigName = false;
#if VCPP
	low_disk_space = false;
#else
	low_disk_space = true;
#endif
}

void SolveParameters::Read_SolveParameters(string solveParametersConfigNameIn)
{
	solveParametersConfigName = solveParametersConfigNameIn;
	fstream in(solveParametersConfigName.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "solveParametersConfigName\t" << solveParametersConfigName << '\n';
		THROW("Cannot open file\n");
	}
	string buf;///z
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "configName")
		{
			READ_NSTRING(in, buf, configName);
		}
		else if (buf == "configPPName")
		{
			READ_NSTRING(in, buf, configPPName);
		}
		else if (buf == "version_configMakerGenName")
		{
			READ_NSTRING(in, buf, version_configMakerGenName);
		}
		else if (buf == "version_configMaker_forceRewrite")
		{
			READ_NBOOL(in, buf, version_configMaker_forceRewrite);
		}
		else if (buf == "versionOffset")
		{
			READ_NINTEGER(in, buf, versionOffset);
		}
		else if (buf == "lv1s2")
		{
			READ_NBOOL(in, buf, lv1s2);
		}
		else if (buf == "low_disk_space")
		{
			READ_NBOOL(in, buf, low_disk_space);
		}
		else if (buf == "serialNumber_st")
		{
			READ_NINTEGER(in, buf, serialNumber_st);
		}
		else if (buf == "serialNumber_en")
		{
			READ_NINTEGER(in, buf, serialNumber_en);
		}
		else if (buf == "versionNumber_st")
		{
			READ_NINTEGER(in, buf, versionNumber_st);
		}
		else if (buf == "versionNumber_en")
		{
			READ_NINTEGER(in, buf, versionNumber_en);
		}
		else if (buf == "np")
		{
			READ_NBOOL(in, buf, numParallelRuns);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	b_solveParametersConfigName = true;
}

void SolveParameters::InitializeAfterSetup()
{
	isDomain = ((sOpt == so_domain_sp) || (sOpt == so_domain_s) || (sOpt == so_domain_p));
	if (configName == "none")
	{
		if (isDomain)
		{
			configName = "config/Domain/Inhomogeneous_axt_fracture_Coh.txt";
			configPPName = "config/Domain/PPS3Config_Coh.txt";
		}
		else if (sOpt == so_domain_p2)
			configName = "config/Domain/PPS2ConfigTest.txt";
		else if (sOpt == so_onePoint_s)
		{
			configName = "config/OneInterface/SampleConfig.txt";
			configName = "config/OneInterface/SampleConfig_TSR_Ortiz.txt";
		}
		else if (sOpt == so_onePoint_s)
		{
			//	configName= "TestFiles/TestPoint_Fracture";
			configName = "TestFiles/TestPoint_No_Fracture";
		}
		else if (sOpt == so_elfrac_fields)
		{
			configName = "TestFiles/TestInhomogeneousElasticFractorField_config.txt";
		}
		else if (sOpt == so_one_field)
		{
			//	configName = "TestFiles/file_v1_sz0_1";
			//	configName = "TestFiles/file_v0_sz0_1";
			//	configName = "TestFiles/file_v1_sz0_x_1";
			configName = "TestFiles/file_v0_sz0_x_1";
			configName = "TestFiles/file0";
			configName = "TestFiles/file_nox_0";
		}
	}
	if ((sOpt == so_onePoint_s) || (sOpt == so_one_field))
		configName = removeExtension(configName);
	if (serialNumber_en < 0)
		serialNumber_en = serialNumber_st;
}
