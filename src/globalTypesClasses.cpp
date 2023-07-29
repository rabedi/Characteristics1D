#include "globalTypesClasses.h"
#include "globalFunctions.h"
#include "commonMacros.h"
#include "RandomVariable.h"

SolveParameters solvePara;

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
	if (g_versionNumber < 0)
		nameBase = g_prefileName + posStr + "_" + specificName;
	else
	{
		string fldr = "../" + posStr;
		MakeDir(fldr);
		nameBase = fldr + "/" + g_prefileNameWOSlash + "_" + posStr + "_" + specificName;
	}
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
	PPS2_outside = true;
	delete_runFolders = 2;
	vis_outside = 2;
#if VCPP
	low_disk_space = 0;
#else
	low_disk_space = 1;
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
			READ_NINTEGER(in, buf, low_disk_space);
			if (low_disk_space < 0)
			{
				if (low_disk_space == -2)
					low_disk_space = 2;
				else if (low_disk_space == -1)
				{
#if VCPP
					low_disk_space = 0;
#else
					low_disk_space = 1;
#endif
				}
			}
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
		else if (buf == "PPS2_outside")
		{
			READ_NBOOL(in, buf, PPS2_outside);
		}
		else if (buf == "vis_outside")
		{
			READ_NINTEGER(in, buf, vis_outside);
		}
		else if (buf == "delete_runFolders")
		{
			READ_NINTEGER(in, buf, delete_runFolders);
			if (delete_runFolders == 3)
			{
#if VCPP
				delete_runFolders = 0;
#else
				delete_runFolders = 1;
#endif
			}
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


///
ostream& operator<<(ostream &output, extermumhold &exter)
{
	output << setprecision(22);
	output << "val\t";
	output << exter.value,
	output << "loc\t";
	output << exter.x_loc;
	return output;
}

istream& operator>>(istream &input, extermumhold &exter)
{
	char buff[25];
	input >> buff;
	input >> exter.value;
	input >> buff;
	input >> exter.x_loc;
	exter.e_initialized = true;
	return input;
}

extermumhold::extermumhold()
{
	value = 0.0;
	x_loc = 0.0;
	e_initialized = false;
}

bool extermumhold::updateMin(double value_in, double x_loc_in)
{
	if (e_initialized == false)
	{
		value = value_in;
		x_loc = x_loc_in;
		e_initialized = true;
		return true;
	}
	if (value_in < value)
	{
		value = value_in;
		x_loc = x_loc_in;
		return true;
	}
	return false;
}

bool extermumhold::updateMax(double value_in, double x_loc_in)
{
	if (e_initialized == false)
	{
		value = value_in;
		x_loc = x_loc_in;
		e_initialized = true;
		return true;
	}
	if (value_in > value)
	{
		value = value_in;
		x_loc = x_loc_in;
		return true;
	}
	return false;
}

void extermumhold::MergeMin(const extermumhold& exte_in)
{
	if (exte_in.e_initialized == false)
		return;
	if ((e_initialized == false) || (exte_in.value < value))
	{
		*this = exte_in; 
		return;
	}
}

void extermumhold::MergeMax(const extermumhold& exte_in)
{
	if (exte_in.e_initialized == false)
		return;
	if ((e_initialized == false) || (exte_in.value > value))
	{
		*this = exte_in;
		return;
	}
}

///////////////////////
ostream& operator<<(ostream &output, const statHolder &stat)
{
	output << "useMeasure" << stat.useMeasure << "\tmeasuresum\t" << stat.sumMeasure << '\n';
	output << setprecision(22);
	output << stat.name << '\t' << stat.nameLatex << '\t' << "count\t";
	output << stat.counter << '\n';
	if (stat.counter == 0)
		return output;

	output << "sum\t" << stat.sum << "\tAve\t" << stat.getAverage();
	output << "\tmax\t" << stat.max.value << "\tx_loc\t" << stat.max.x_loc;
	output << "\tmin\t" << stat.min.value << "\tx_loc\t" << stat.min.x_loc;
	output << "\tSndrdDvt\t" << stat.getStandardDeviation() << endl;
	return output;
}

istream& operator>>(istream &input, statHolder &stat)
{
	char buff[255];
	input >> buff;
	stat.useMeasure = false;
	stat.sumMeasure = 0.0;
	if (strcmp(buff, "useMeasure0") == 0)
	{
		stat.useMeasure = false;
		input >> stat.sumMeasure;
		input >> stat.name;
		input >> stat.nameLatex;

	}
	else if (strcmp(buff, "useMeasure1") == 0)
	{
		stat.useMeasure = true;
		input >> stat.sumMeasure;
		input >> stat.name;
		input >> stat.nameLatex;
	}
	else
	{
		stat.name = buff;
		stat.nameLatex = buff;
	}
	input >> buff;
	input >> stat.counter;
	if (stat.counter == 0)
		return input;

	double average, sdiv;
	input >> buff >> stat.sum >> buff >> average;
	input >> buff >> stat.max.value >> buff >> stat.max.x_loc;
	stat.max.e_initialized = true;
	input >> buff >> stat.min.value >> buff >> stat.min.x_loc;
	stat.min.e_initialized = true;
	input >> buff >> sdiv;
	stat.setSquareSumFromStandardDeviation(sdiv);
	return input;
}

statHolder::statHolder(bool useMeasureIn)
{
	setName("none", "none");
	counter = 0;
	sumMeasure = 0.0;
	useMeasure = useMeasureIn;
	sum = 0.0;
	sumSquares = 0.0;

	rN = 0.0;
	sdiv_saved = 0.0;
	b_sdiv_saved = false;
}

statHolder::statHolder(string& nameIn, string& nameLatexIn, bool useMeasureIn)
{
	setName(nameIn, nameLatexIn);
	counter = 0;
	sumMeasure = 0.0;
	useMeasure = useMeasureIn;
	sum = 0.0;
	sumSquares = 0.0;
}

void statHolder::setName(const string& nameIn, const string& nameLatexIn)
{
	name = nameIn;
	nameLatex = nameLatexIn;
}

double statHolder::getAverage() const
{
	if (counter > 0)
		return (double)(sum / get_measure());
	return 1e40;
}

double statHolder::getStandardDeviation() const
{
	if (b_sdiv_saved)
		return sdiv_saved;
	if (counter > 0)
	{
		double measureSum = (double)get_measure();
		return sqrt(fabs(sumSquares - sum * sum / measureSum) / measureSum);
	}
	return 1e40;
}

double statHolder::getCOV() const
{
	if (counter > 0)
		return computeRatio(getStandardDeviation(), getAverage());
	return 1e40;
}

void statHolder::update(double value_in, double x_loc_in, double weightIn)
{
	max.updateMax(value_in, x_loc_in); min.updateMin(value_in, x_loc_in);
	++counter;
	if (!useMeasure)
	{
		sum += value_in;
		sumSquares += (value_in * value_in);
	}
	else
	{
		sumMeasure += weightIn;
		double tmp = value_in * weightIn;
		sum += tmp;
		sumSquares += (tmp * value_in);
	}
}

void statHolder::MergeStatHolder(statHolder& stat)
{
	max.MergeMax(stat.max);
	min.MergeMin(stat.min);
	counter += stat.counter;
	sumMeasure += stat.sumMeasure;
	useMeasure = stat.useMeasure;
	if (counter == 0)
		return;
	sum += stat.sum;
	sumSquares += stat.sumSquares;
}

double statHolder::getValue(int setStatOp_type_i)
{
	setStatOp_type stat_sv_type = (setStatOp_type)setStatOp_type_i;
	if (stat_sv_type == sso_number)
		return getCount();
	if (stat_sv_type == sso_mean_arithmetic)
		return getAverage();
	if (stat_sv_type == sso_sdiv)
		return getStandardDeviation();
	if (stat_sv_type == sso_min)
		return getMin();
	if (stat_sv_type == sso_max)
		return getMax();
	if (stat_sv_type == sso_normalized_number)
		return rN;
	if (stat_sv_type == sso_index_min)
		return getMinLoc();
	if (stat_sv_type == sso_index_max)
		return getMaxLoc();
	if (stat_sv_type == sso_cov)
		return getCOV();
	cout << "stat_sv_type\t" << stat_sv_type << '\n';
	THROW("stat_sv_type\n");
}

void statHolder::setSquareSumFromStandardDeviation(double sDeviation)
{
	if (counter > 0)
	{
		if (!useMeasure)
			sumSquares = (counter * sDeviation *  sDeviation + sum * sum / counter);
		else
			sumSquares = (sumMeasure * sDeviation *  sDeviation + sum * sum / sumMeasure);
	}
	else
		sumSquares = 0.0;
}
