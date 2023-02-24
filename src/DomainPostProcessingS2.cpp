#include "DomainPostProcessingS2.h"
#include "globalFunctions.h"
#include "SLDescriptorData.h"

/////////////////
string getName(loadingStagesT dat)
{
	if (dat == loadingStagesT_none)
		return "none";
	if (dat == lst_actualTime)
		return "actualTime";
	if (dat == lst_lin)
		return "l";
	if (dat == lst_max)
		return "m";
	if (dat == lst_final)
		return "f";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, loadingStagesT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= loadingStagesT_SIZE)
			THROW("too large of a number\n");
		typeVal = (loadingStagesT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < loadingStagesT_SIZE; ++i)
	{
		typeVal = (loadingStagesT)i; // casting integer to loadingStagesT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading loadingStagesT\n");
}

//operator for output
ostream& operator<<(ostream& out, loadingStagesT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, loadingStagesT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getName(lstClassificationT dat)
{
	if (dat == lstClassificationT_none)
		return "none";
	if (dat == ct_time)
		return "time";
	if (dat == ct_epssig_bc)
		return "epssig_bc";
	if (dat == ct_epssig_bulk)
		return "epssig_bulk";
	if (dat == ct_interfaceDissEne)
		return "interfaceDissEne";
	if (dat == ct_Damage)
		return "Damage";
	if (dat == ct_U)
		return "U";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, lstClassificationT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= lstClassificationT_SIZE)
			THROW("too large of a number\n");
		typeVal = (lstClassificationT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < lstClassificationT_SIZE; ++i)
	{
		typeVal = (lstClassificationT)i; // casting integer to lstClassificationT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading lstClassificationT\n");
}

//operator for output
ostream& operator<<(ostream& out, lstClassificationT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, lstClassificationT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getLatexName(BrittlenessIndicatorRatioType dat)
{
	if (dat == BrittlenessIndicatorRatioType_none)
		return "none";
	if (dat == birt_i2m)
		return "{\\mathrm{i/m}}";
	if (dat == birt_i2f)
		return "{\\mathrm{i/f}}";
	if (dat == birt_m2f)
		return "{\\mathrm{m/f}}";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

string getName(BrittlenessIndicatorRatioType dat)
{
	if (dat == BrittlenessIndicatorRatioType_none)
		return "none";
	if (dat == birt_i2m)
		return "i,m";
	if (dat == birt_i2f)
		return "i,f";
	if (dat == birt_m2f)
		return "m,f";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, BrittlenessIndicatorRatioType& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= BrittlenessIndicatorRatioType_SIZE)
			THROW("too large of a number\n");
		typeVal = (BrittlenessIndicatorRatioType)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < BrittlenessIndicatorRatioType_SIZE; ++i)
	{
		typeVal = (BrittlenessIndicatorRatioType)i; // casting integer to BrittlenessIndicatorRatioType, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading BrittlenessIndicatorRatioType\n");
}

//operator for output
ostream& operator<<(ostream& out, BrittlenessIndicatorRatioType dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

istream& operator>>(istream& in, BrittlenessIndicatorRatioType& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getLatexName(BrittlnessIndicatorFieldType dat)
{
	if (dat == BrittlnessIndicatorFieldType_none)
		return "none";
	if (dat == birt_time)
		return "t";
	if (dat == birt_eps)
		return "\\epsilon";
	if (dat == birt_psi)
		return "\\psi";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

string getName(BrittlnessIndicatorFieldType dat)
{
	if (dat == BrittlnessIndicatorFieldType_none)
		return "none";
	if (dat == birt_time)
		return "time";
	if (dat == birt_eps)
		return "strain";
	if (dat == birt_psi)
		return "energy";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, BrittlnessIndicatorFieldType& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= BrittlnessIndicatorFieldType_SIZE)
			THROW("too large of a number\n");
		typeVal = (BrittlnessIndicatorFieldType)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < BrittlnessIndicatorFieldType_SIZE; ++i)
	{
		typeVal = (BrittlnessIndicatorFieldType)i; // casting integer to BrittlnessIndicatorFieldType, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading BrittlnessIndicatorFieldType\n");
}

//operator for output
ostream& operator<<(ostream& out, BrittlnessIndicatorFieldType dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, BrittlnessIndicatorFieldType& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getLatexName(FragmentationCriterionT dat)
{
	if (dat == FragmentationCriterionT_none)
		return "none";
	if (dat == fct_Damage)
		return "D";
	if (dat == fct_maxEffDelU)
		return "\\Delta{u}_{\\mathrm{max}}";
	if (dat == fct_DelU)
		return "\\Delta{u}";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

string getName(FragmentationCriterionT dat)
{
	if (dat == FragmentationCriterionT_none)
		return "none";
	if (dat == fct_Damage)
		return "D";
	if (dat == fct_maxEffDelU)
		return "Max_DelU";
	if (dat == fct_DelU)
		return "DelU";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, FragmentationCriterionT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= FragmentationCriterionT_SIZE)
			THROW("too large of a number\n");
		typeVal = (FragmentationCriterionT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < FragmentationCriterionT_SIZE; ++i)
	{
		typeVal = (FragmentationCriterionT)i; // casting integer to FragmentationCriterionT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading FragmentationCriterionT\n");
}

//operator for output
ostream& operator<<(ostream& out, FragmentationCriterionT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, FragmentationCriterionT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getLatexName(RunNormalizationQuantT dat)
{
	if (dat == RunNormalizationQuantT_none)
		return "none";
	if (dat == RunNormalizationQuantT_default)
		return "1";
	if (dat == ene_KInitial)
		return "K_0";
	if (dat == ene_UInitial)
		return "U_0";
	if (dat == ene_PhiInitial)
		return "\\phi_0";
	if (dat == ene_KFinal)
		return "K_f";
	if (dat == ene_UFinal)
		return "U_f";
	if (dat == ene_PhiFinal)
		return "\\phi_f";
	if (dat == eneBCFinal)
		return "\\mathcal{E}_{\\mathrm{bc},f}";
	if (dat == eneInputFinal)
		return "\\mathcal{E}_{\\mathrm{input},f}";
	if (dat == rnq_loadTimeScale)
		return "{t_l}";
	if (dat == rnq_EScale)
		return "E";
	if (dat == rnq_rhoScale)
		return "{\rho}";
	if (dat == rnq_dampingScale)
		return "{D}";
	if (dat == rnq_sigmaCScale)
		return "{\\sigma_c}";
	if (dat == rnq_deltaCScale)
		return "{\\delta_c}";
	if (dat == rnq_energyCScale)
		return "{G_c}";
	if (dat == rnq_cScale)
		return "c";
	if (dat == rnq_ZScale)
		return "Z";
	if (dat == rnq_tauCScale)
		return "{\tau_c}";
	if (dat == rnq_strainCScale)
		return "{\\epsilon_c}";
	if (dat == rnq_velCScale)
		return "{v_c}";
	if (dat == rnq_psiCScale)
		return "\\psi";
	if (dat == RunNormalizationQuantT_SIZE)
		return "size";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

string getName(RunNormalizationQuantT dat)
{
	if (dat == RunNormalizationQuantT_none)
		return "none";
	if (dat == RunNormalizationQuantT_default)
		return "default";
	if (dat == ene_KInitial)
		return "ene_KInitial";
	if (dat == ene_UInitial)
		return "ene_UInitial";
	if (dat == ene_PhiInitial)
		return "ene_PhiInitial";
	if (dat == ene_KFinal)
		return "ene_KFinal";
	if (dat == ene_UFinal)
		return "ene_UFinal";
	if (dat == ene_PhiFinal)
		return "ene_PhiFinal";
	if (dat == eneBCFinal)
		return "eneBCFinal";
	if (dat == eneInputFinal)
		return "eneInputFinal";
	if (dat == rnq_loadTimeScale)
		return "loadTime";
	if (dat == rnq_EScale)
		return "E";
	if (dat == rnq_rhoScale)
		return "rho";
	if (dat == rnq_dampingScale)
		return "damping";
	if (dat == rnq_sigmaCScale)
		return "sigmaC";
	if (dat == rnq_deltaCScale)
		return "deltaC";
	if (dat == rnq_energyCScale)
		return "energyC";
	if (dat == rnq_cScale)
		return "c";
	if (dat == rnq_ZScale)
		return "Z";
	if (dat == rnq_tauCScale)
		return "tauC";
	if (dat == rnq_strainCScale)
		return "strainC";
	if (dat == rnq_velCScale)
		return "velC";
	if (dat == rnq_psiCScale)
		return "psiC";
	if (dat == RunNormalizationQuantT_SIZE)
		return "size";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, RunNormalizationQuantT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= RunNormalizationQuantT_SIZE)
			THROW("too large of a number\n");
		typeVal = (RunNormalizationQuantT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < RunNormalizationQuantT_SIZE; ++i)
	{
		typeVal = (RunNormalizationQuantT)i; // casting integer to RunNormalizationQuantT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading RunNormalizationQuantT\n");
}

//operator for output
ostream& operator<<(ostream& out, RunNormalizationQuantT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, RunNormalizationQuantT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getName(timeIndexType dat)
{
	if (dat == timeIndexType_none)
		return "none";
	if (dat == ti_AllTimeSteps)
		return "ti_AllTimeSteps";
	if (dat == ti_DSU_Fragment)
		return "ti_DSU_Fragment";
	if (dat == ti_BulkInterfacePoints)
		return "ti_BulkInterfacePoints";
	if (dat == ti_RawFinalSln_Scalars)
		return "ti_RawFinalSln_Scalars";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, timeIndexType& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= timeIndexType_SIZE)
			THROW("too large of a number\n");
		typeVal = (timeIndexType)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < timeIndexType_SIZE; ++i)
	{
		typeVal = (timeIndexType)i; // casting integer to timeIndexType, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading timeIndexType\n");
}

//operator for output
ostream& operator<<(ostream& out, timeIndexType dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, timeIndexType& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

istream& operator>>(istream& in, PPS2_TimeStamp& dat)
{
	dat.PPS2_TimeStamp_Read(in);
	return in;
}

PPS2_TimeStamp::PPS2_TimeStamp()
{
	MakeVoid_PPS2_TimeStamp();
}

void PPS2_TimeStamp::MakeVoid_PPS2_TimeStamp()
{
	timeStampType = loadingStagesT_none;
	actualTime_indexType = timeIndexType_none;
	actualTime_val = 0.0;
	actualTime_index = 0;
	lct = lstClassificationT_none;
}

void PPS2_TimeStamp::PPS2_TimeStamp_Read(istream& in)
{
	MakeVoid_PPS2_TimeStamp();
	string buf;
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
		if (buf == "timeStampType")
		{
			in >> timeStampType;
		}
		else if (buf == "actualTime_indexType")
		{
			in >> actualTime_indexType;
		}
		else if (buf == "actualTime_val")
		{
			in >> actualTime_val;
		}
		else if (buf == "actualTime_index")
		{
			in >> actualTime_index;
		}
		else if (buf == "lct")
		{
			in >> lct;
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

PPS2_dataPointer::PPS2_dataPointer()
{
	MakeVoid_PPS2_dataPointer();
}

void PPS2_dataPointer::MakeVoid_PPS2_dataPointer()
{
	name_In_CSV_file = "none";
	name_Latex = "none";
	overwriterName = "";

	isActive = true;
	subdomainNo = -1;
	stat_sv_type = sso_none;
	timeStamp.MakeVoid_PPS2_TimeStamp();
	fldName = "";
	fragmentationCriterion = FragmentationCriterionT_none;
	brittlenessRatioT = BrittlenessIndicatorRatioType_none;
	brittlenessFieldT = BrittlnessIndicatorFieldType_none;
	normalizationMode = RunNormalizationQuantT_none;
	normalizationMode = RunNormalizationQuantT_none;
	oprType = vpo_none;
}

bool PPS2_dataPointer::PPS2_dataPointer_Read(istream& in)
{
	MakeVoid_PPS2_dataPointer();
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "isActive")
		{
			READ_NBOOL(in, buf, isActive);
		}
		else if (buf == "name_In_CSV_file")
		{
			READ_NSTRING(in, buf, name_In_CSV_file);
		}
		else if (buf == "name_Latex")
		{
			READ_NSTRING(in, buf, name_Latex);
		}
		else if (buf == "overwriterName")
		{
			READ_NSTRING(in, buf, overwriterName);
		}
		else if (buf == "fldName")
		{
			READ_NSTRING(in, buf, fldName);
		}
		else if (buf == "subdomainNo")
		{
			READ_NINTEGER(in, buf, subdomainNo);
		}
		else if (buf == "stat_sv_type")
			in >> stat_sv_type;
		else if (buf == "timeStamp")
			in >> timeStamp;
		else if (buf == "brittlenessRatioT")
		{
			in >> brittlenessRatioT;
			if (brittlenessFieldT == BrittlnessIndicatorFieldType_none)
				brittlenessFieldT = birt_eps;
		}
		else if (buf == "brittlenessFieldT")
		{
			in >> brittlenessFieldT;
			if (brittlenessRatioT == BrittlenessIndicatorRatioType_none)
				brittlenessRatioT = birt_m2f;
		}
		else if (buf == "fragmentationCriterion")
			in >> fragmentationCriterion;
		else if (buf == "normalizationMode")
			in >> normalizationMode;
		else if (buf == "oprType")
			in >> oprType;
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	return true;
}

//////////////////////////////////////////////////////////////
PPS2_dataPointer_StageOverwritter::PPS2_dataPointer_StageOverwritter()
{
	MakeVoid_PPS2_dataPointer_StageOverwritter();
}

void PPS2_dataPointer_StageOverwritter::MakeVoid_PPS2_dataPointer_StageOverwritter()
{
	post_name_In_CSV_file = "";
	post_name_Latex = "";

	isActive = -1;
	timeStamp.MakeVoid_PPS2_TimeStamp();
}

bool PPS2_dataPointer_StageOverwritter::PPS2_dataPointer_StageOverwritter_Read(istream& in)
{
	MakeVoid_PPS2_dataPointer_StageOverwritter();
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			cout << "buf\t" << buf << '\n';
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "isActive")
		{
			READ_NINTEGER(in, buf, isActive);
		}
		else if (buf == "post_name_In_CSV_file")
		{
			READ_NSTRING(in, buf, post_name_In_CSV_file);
		}
		else if (buf == "post_name_Latex")
		{
			READ_NSTRING(in, buf, post_name_Latex);
		}
		else if (buf == "timeStamp")
			in >> timeStamp;
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	return true;
}

void PPS2_dataPointer_StageOverwritter::Overwrite_PPS2_dataPointer(PPS2_dataPointer & datPointer)
{
	if ((isActive >= 0) && (datPointer.isActive != false))
		datPointer.isActive = (bool)isActive;
	datPointer.timeStamp = timeStamp;
	datPointer.name_In_CSV_file += post_name_In_CSV_file;
	datPointer.name_Latex += post_name_Latex;
}

//////////////////////////////////////////////////////////////

DataW4LineHeader::DataW4LineHeader()
{
	numDataPoints = 0;
	numFlds = 0;
}

void DataW4LineHeader::Initialize_DataW4LineHeader(string& nameIn, DataW4LineHeader& statOut, string& nameOut)
{
	DataW4LineHeader_ReadFileName(nameIn);
	Compute_Stat(statOut);
	fstream out(nameOut.c_str(), ios::out);
	statOut.DataW4LineHeader_Write(out);
}

void DataW4LineHeader::GetMaxAbs_FinalValue(const DataW4LineHeader& timeSequenceSummaryStat, unsigned int fldNo, bool classification_usePositive4MaxValue, bool forceTerminalvalue2ProvidedValue, double providedTerminalValue, bool forwardFinalValueSearch, double tol4FinalValCheck,
	bool &bValuesPositive, int& index_maxAbsVal, double& maxAbsVal, int& index_finalVal, double& finalVal)
{
	// value with maximum absolute value
	double minV = timeSequenceSummaryStat.data_vals[sso_min][fldNo], maxV = timeSequenceSummaryStat.data_vals[sso_max][fldNo];
	double maxAbs;
//	if (fabs(minV) < fabs(maxV))
	if (classification_usePositive4MaxValue)
	{
		bValuesPositive = true;
		index_maxAbsVal = (int)floor(timeSequenceSummaryStat.data_vals[sso_index_max][fldNo] + 0.2);
		maxAbsVal = maxV;
		maxAbs = maxV;
	}
	else
	{
		bValuesPositive = false;
		index_maxAbsVal = (int)floor(timeSequenceSummaryStat.data_vals[sso_index_min][fldNo] + 0.2);
		maxAbsVal = minV;
		maxAbs = -minV;
	}

	//// terminal value
	double absTol;
	if (tol4FinalValCheck < 0)
		absTol = -tol4FinalValCheck * maxAbs;
	else
		absTol = tol4FinalValCheck;

	double terminalVal = providedTerminalValue;
	unsigned int index_lastPt = numDataPoints - 1;
	if (!forceTerminalvalue2ProvidedValue)
	{
		terminalVal = data_vals[index_lastPt][fldNo];
		--index_lastPt;
		if (index_lastPt < 0)
		{
			finalVal = terminalVal;
			index_finalVal = -1;
		}
	}

	if (index_maxAbsVal == numDataPoints - 1)
	{
		index_finalVal = -1;
		finalVal = data_vals[numDataPoints - 1][fldNo];
		return;
	}
	if (forwardFinalValueSearch)
	{
		for (index_finalVal = index_maxAbsVal + 1; index_finalVal <= (int)index_lastPt; ++index_finalVal)
		{
			finalVal = data_vals[index_finalVal][fldNo];
			if (fabs(finalVal - terminalVal) < absTol)
				return;
		}
		index_finalVal = -1;
		return;
	}
	///// going backward ... no need for else as previous if returns
//	bool success = false;
	for (index_finalVal = index_lastPt; index_finalVal > index_maxAbsVal; --index_finalVal)
	{
		finalVal = data_vals[index_finalVal][fldNo];
		if (fabs(finalVal - terminalVal) > absTol)
		{
//			success = true;
			break;
		}
	}
	if (index_finalVal == index_lastPt)
	{
		index_finalVal = -1;
		return;
	}
	++index_finalVal;
	finalVal = data_vals[index_finalVal][fldNo];
}

bool DataW4LineHeader::GetFirstLocationGreaterThanValue(const DataW4LineHeader& timeSequenceSummaryStat, unsigned int fldNo, double providedTerminalValue, double tol4FinalValCheck,
	int& index)
{
	index = -1;
	double maxV2Check = providedTerminalValue - tol4FinalValCheck;
	double maxV = timeSequenceSummaryStat.data_vals[sso_max][fldNo];
	if (maxV < maxV2Check)
		return false;
	double val;
	for (index = 0; index <= (int)numDataPoints; ++index)
	{
		val = data_vals[index][fldNo];
		if (val >= maxV2Check)
			return true;
	}
	return false;
}

void DataW4LineHeader::GetLinValue(const DataW4LineHeader& timeSequenceSummaryStat, unsigned int fldNo, double tol4ZeroValCheck, int& index_LinVal)
{
	double absTol;
	if (tol4ZeroValCheck < 0)
	{
		double maxAbs = MAX(fabs(timeSequenceSummaryStat.data_vals[sso_min][fldNo]), fabs(timeSequenceSummaryStat.data_vals[sso_max][fldNo]));
		absTol = -tol4ZeroValCheck * maxAbs;
	}
	else
		absTol = tol4ZeroValCheck;

	if (absTol < DBL_MIN)
	{
		index_LinVal = 0;
		return;
	}
	double val;
	for (index_LinVal = 0; index_LinVal < (int)numDataPoints; ++index_LinVal)
	{
		val = data_vals[index_LinVal][fldNo];
		if (fabs(val) > absTol)
			break;
	}
	if ((index_LinVal == 0) || (index_LinVal >= (int)(numDataPoints - 1)))
		index_LinVal = -1;
	else
		--index_LinVal;
}

void DataW4LineHeader::GetIndexClosestSortingVal(unsigned int sortingFldNo, double sortingValIn, unsigned int& closestSortingIndex, double& closestSortingValOut)
{
	if (sortingValIn <= data_vals[0][sortingFldNo])
	{
		closestSortingIndex = 0;
		closestSortingValOut = data_vals[0][sortingFldNo];
		return;
	}
	if (sortingValIn >= data_vals[numDataPoints - 1][sortingFldNo])
	{
		closestSortingIndex = numDataPoints - 1;
		closestSortingValOut = data_vals[closestSortingIndex][sortingFldNo];
		return;
	}

	for (closestSortingIndex = 0; closestSortingIndex < numDataPoints; ++closestSortingIndex)
	{
		closestSortingValOut = data_vals[closestSortingIndex][sortingFldNo];
		if (closestSortingValOut > sortingValIn)
			break;
	}
	if (closestSortingIndex == numDataPoints)
	{
		closestSortingIndex = numDataPoints - 1;
		closestSortingValOut = data_vals[closestSortingIndex][sortingFldNo];
		return;
	}
	double diff = closestSortingValOut - sortingValIn;
	double diff2 = sortingValIn - data_vals[closestSortingIndex - 1][sortingFldNo];
	if (diff2 < diff)
	{
		--closestSortingIndex;
		closestSortingValOut = data_vals[closestSortingIndex][sortingFldNo];
	}
}

bool DataW4LineHeader::DataW4LineHeader_ReadFileName(string& fileNameIn)
{
	fstream in(fileNameIn.c_str(), ios::in);
	if (!in.is_open())
		return false;
	getline(in, header1);
	numFlds = BreakString(header1, headerStrs1);
	getline(in, header2);
	BreakString(header2, headerStrs2);
	getline(in, header3);
	BreakString(header3, headerStrs3);
	getline(in, header4);
	BreakString(header4, headerStrs4);

	double tmp;
	in >> tmp;
	unsigned int cntr = 0;
	while (!in.eof())
	{
		vector<double> val(numFlds);
		val[0] = tmp;
		for (unsigned int i = 1; i < numFlds; ++i)
			in >> val[i];

		if (!in.eof())
		{
			++cntr;
			data_vals.push_back(val);
		}
		in >> tmp;
	}
	numDataPoints = data_vals.size();
	return true;
}
void DataW4LineHeader::Compute_Stat(DataW4LineHeader& statSummary)
{
	statSummary.numDataPoints = setStatOp_type_SIZE_SHORT;
	statSummary.data_vals.resize(setStatOp_type_SIZE_SHORT);
	statSummary.header1 = header1;
	statSummary.header2 = header2;
	statSummary.header3 = header3;
	statSummary.header4 = header4;
	statSummary.headerStrs1 = headerStrs1;
	statSummary.headerStrs2 = headerStrs2;
	statSummary.headerStrs3 = headerStrs3;
	statSummary.headerStrs4 = headerStrs4;
	statSummary.numDataPoints = setStatOp_type_SIZE_SHORT;
	statSummary.numFlds = numFlds;
	for (unsigned int i = 0; i < statSummary.numDataPoints; ++i)
		statSummary.data_vals[i].resize(numFlds);

	vector<double> *means, *mins, *maxs, *vars, *sdivs, *covs, *index_mins, *index_maxs;
	statSummary.data_vals[sso_valStart] = data_vals[0];
	statSummary.data_vals[sso_valEnd] = data_vals[numDataPoints - 1];
	means = &statSummary.data_vals[sso_mean_arithmetic];
	mins = &statSummary.data_vals[sso_min];
	maxs = &statSummary.data_vals[sso_max];
	vars = &statSummary.data_vals[sso_var];
	sdivs = &statSummary.data_vals[sso_sdiv];
	covs = &statSummary.data_vals[sso_cov];
	index_mins = &statSummary.data_vals[sso_index_min];
	index_maxs = &statSummary.data_vals[sso_index_max];

	fill(means->begin(), means->end(), 0.0);
	fill(vars->begin(), vars->end(), 0.0);
	fill(sdivs->begin(), sdivs->end(), 0.0);
	fill(covs->begin(), covs->end(), 0.0);
	fill(mins->begin(), mins->end(), 0.9 * DBL_MAX);
	fill(maxs->begin(), maxs->end(), -0.9 * DBL_MAX);
	fill(index_mins->begin(), index_mins->end(), -1);
	fill(index_maxs->begin(), index_maxs->end(), -1);

	double tmp;
	for (unsigned int di = 0; di < numDataPoints; ++di)
	{
		vector<double>* tmpV = &data_vals[di];
		for (unsigned int fi = 0; fi < numFlds; ++fi)
		{
			tmp = (*tmpV)[fi];
			if (tmp < (*mins)[fi])
			{
				(*mins)[fi] = tmp;
				(*index_mins)[fi] = (double)di;
			}
			if (tmp > (*maxs)[fi])
			{
				(*maxs)[fi] = tmp;
				(*index_maxs)[fi] = (double)di;
			}
			(*means)[fi] += tmp;
		}
	}
	double sizeFactor = 1.0 / (double)numDataPoints;
	for (unsigned int fi = 0; fi < numFlds; ++fi)
		(*means)[fi] *= sizeFactor;
	if (numDataPoints > 1)
	{
		sizeFactor = 1.0 / (double)(numDataPoints - 1.0);
		for (unsigned int di = 0; di < numDataPoints; ++di)
		{
			vector<double>* tmpV = &data_vals[di];
			for (unsigned int fi = 0; fi < numFlds; ++fi)
			{
				tmp = (*tmpV)[fi] - (*means)[fi];
				(*vars)[fi] += tmp * tmp;
			}
		}
		for (unsigned int fi = 0; fi < numFlds; ++fi)
		{
			(*vars)[fi] *= sizeFactor;
			(*sdivs)[fi] = sqrt((*vars)[fi]);
			(*covs)[fi] = computeRatio((*sdivs)[fi], (*means)[fi]);
		}
	}
}

void DataW4LineHeader::DataW4LineHeader_Write(ostream& out)
{
	out << header1 << '\n';
	out << header2 << '\n';
	out << header3 << '\n';
	out << header4 << '\n';
	unsigned int numFldsm1 = numFlds - 1;
	vector<double> *vals;
	for (unsigned int i = 0; i < numDataPoints; ++i)
	{
		vals = &data_vals[i];
		for (unsigned int fi = 0; fi < numFlds; ++fi)
		{
			out << (*vals)[fi];
			if (fi != numFldsm1)
				out << '\t';
			else
				out << '\n';
		}
	}
}

string getBrittleness_IndicatorName(bool isLatex, BrittlnessIndicatorFieldType bf, BrittlenessIndicatorRatioType br)
{
	string name = "", name_bf, name_br;
	if (isLatex)
	{
		name_bf = getLatexName(bf);
		name_br = getLatexName(br);
		name = "B^" + name_bf + "_" + name_br;
		return name;
	}
	name_bf = getName(bf);
	name_br = getName(br);
	name = "B_" + name_bf + "_" + name_br;
	return name;
}

OneTimeValuePPS2Data::OneTimeValuePPS2Data()
{
#if 0
	subdomainNo = 0;
	ct = lstClassificationT_none;
	lst = loadingStagesT_none;
#endif

	timeIndex_4Stage = -1;
	timeValue_4Stage = -1.0;

	timeIndex_DSU_Fragment4Stage = -1;
	timeValue_DSU_Fragment4Stage = -1.0;

	timeIndex_BulkInterfacePoints4Stage = -1;
	timeValue_BulkInterfacePoints4Stage = -1.0;
}

void OneTimeValuePPS2Data::SetAllTimeValsIndices_FromActualTime(OneSubdomainPPS2Data_runInfo segmentInfo)
{
	if (timeValue_4Stage < 0.0)
	{
		timeIndex_4Stage = -1;
		timeIndex_DSU_Fragment4Stage = -1;
		timeValue_DSU_Fragment4Stage = -1.0;

		timeIndex_BulkInterfacePoints4Stage = -1;
		timeValue_BulkInterfacePoints4Stage = -1.0;
		return;
	}
//	timeIndex_DSU_Fragment4Stage
	timeIndex_DSU_Fragment4Stage = segmentInfo.get_Interface_DSU_Fragment_TimeIndexValue(timeValue_4Stage, timeValue_DSU_Fragment4Stage);
	timeIndex_BulkInterfacePoints4Stage = segmentInfo.get_Interface_BulkInterfacePoints_TimeIndexValue(timeValue_4Stage, timeValue_BulkInterfacePoints4Stage);
}

OneClassificationPPS2Data::OneClassificationPPS2Data()
{
#if 0
	subdomainNo = 0;
	ct = lstClassificationT_none;
#endif

	for (unsigned int bf = 0; bf < BrittlnessIndicatorFieldType_SIZE; ++bf)
		for (unsigned int br = 0; br < BrittlenessIndicatorRatioType_SIZE; ++br)
			brittlenessIndicators[bf][br] = invalidNonnegativeNum;
}

void OneClassificationPPS2Data::BrittlenessIndicators_Write(ostream& out)
{
	for (unsigned int br = 0; br < BrittlenessIndicatorRatioType_SIZE; ++br)
		out << "\t" << getName((BrittlenessIndicatorRatioType)br);
	out << '\n';
	for (unsigned int bf = 0; bf < BrittlnessIndicatorFieldType_SIZE; ++bf)
	{
		out	<< getName((BrittlnessIndicatorFieldType)bf);
		for (unsigned int br = 0; br < BrittlenessIndicatorRatioType_SIZE; ++br)
			out << '\t' << brittlenessIndicators[bf][br];
		out << '\n';
	}
	bool isLatex = false;
	for (unsigned int bf = 0; bf < BrittlnessIndicatorFieldType_SIZE; ++bf)
	{
		for (unsigned int br = 0; br < BrittlenessIndicatorRatioType_SIZE; ++br)
		{
			string tmp = getBrittleness_IndicatorName(isLatex, (BrittlnessIndicatorFieldType)bf, (BrittlenessIndicatorRatioType)br);
			out << tmp << '\t';
		}
		out << '\n';
	}
	isLatex = true;
	for (unsigned int bf = 0; bf < BrittlnessIndicatorFieldType_SIZE; ++bf)
	{
		for (unsigned int br = 0; br < BrittlenessIndicatorRatioType_SIZE; ++br)
		{
			string tmp = getBrittleness_IndicatorName(isLatex, (BrittlnessIndicatorFieldType)bf, (BrittlenessIndicatorRatioType)br);
			out << tmp << '\t';
		}
		out << '\n';
	}
}

void OneClassificationPPS2Data::BrittlenessIndicators_Read(istream& in)
{
	string buf;
	getline(in, buf);
	for (unsigned int bf = 0; bf < BrittlnessIndicatorFieldType_SIZE; ++bf)
	{
		in >> buf;
		for (unsigned int br = 0; br < BrittlenessIndicatorRatioType_SIZE; ++br)
			in >> brittlenessIndicators[bf][br];
	}
}

OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData()
{
	mean_fragmentSize = 0.0;
	min_fragmentSize = 0.0;
	max_fragmentSize = 0.0;
	sdiv_fragmentSize = 0.0;
	cov_fragmentSize = 0.0;
	numFragments = 0;
}

void OneTimeOneCriterionFragmentationRawData::ComputeFragmentStatsFromFragmentSizes()
{
	mean_fragmentSize = getStatvalue(fragment_lengths, sso_mean_arithmetic);
	min_fragmentSize = getStatvalue(fragment_lengths, sso_min);
	max_fragmentSize = getStatvalue(fragment_lengths, sso_max);
	sdiv_fragmentSize = getStatvalue(fragment_lengths, sso_sdiv);
	cov_fragmentSize = sdiv_fragmentSize / mean_fragmentSize;
}

bool OneTimeOneCriterionFragmentationRawData::Get_Scalar_Or_Vector_Output(setStatOp_type statType, double& scalarVal, vector<double>& vecVal, bool& outputIsScalar)
{
	if (statType == sso_all_field_val)
	{
		vecVal = fragment_lengths;
		outputIsScalar = false;
		return (numFragments > 0);
	}
	if (numFragments <= 0)
	{
		scalarVal = invalidNum;
		return false;
	}
	outputIsScalar = true;
	if (statType == sso_number)
		scalarVal = numFragments;
	else if (statType == sso_min)
		scalarVal = min_fragmentSize;
	else if (statType == sso_mean_arithmetic)
		scalarVal = mean_fragmentSize;
	else if (statType == sso_max)
		scalarVal = max_fragmentSize;
	else if (statType == sso_cov)
		scalarVal = cov_fragmentSize;
	else if (statType == sso_sdiv)
		scalarVal = sdiv_fragmentSize;
	return true;
}

void OneTimeOneCriterionFragmentationRawData::Fragmentation_Sizes_Write(ostream& out, int timeIndex, double timeVal)
{
	out << timeIndex << '\t' << timeVal << '\t';
	out << numFragments;
	for (unsigned int i = 0; i < numFragments; ++i)
		out << '\t' << fragment_lengths[i];
	out << '\n';
}

void OneTimeOneCriterionFragmentationRawData::Fragmentation_Sizes_Read(istream& in)
{
	string buf;
	in >> buf >> buf;
	in >> numFragments;
	fragment_lengths.resize(numFragments);
	for (unsigned int i = 0; i < numFragments; ++i)
		in >> fragment_lengths[i];
}

void OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Write_Data(ostream& out, int timeIndex, double timeVal)
{
	out << timeIndex << '\t';
	out << timeVal << '\t';
	out << numFragments << '\t';
	out << mean_fragmentSize << '\t';
	out << min_fragmentSize << '\t';
	out << max_fragmentSize << '\t';
	out << sdiv_fragmentSize << '\t';
	out << cov_fragmentSize << '\n';
}

void OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Read_Data(istream& in)
{
	in >> numFragments;
	in >> mean_fragmentSize;
	in >> min_fragmentSize;
	in >> max_fragmentSize;
	in >> sdiv_fragmentSize;
	in >> cov_fragmentSize;
}


void OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Write_Header(ostream& out)
{
	vector<HeaderLabels> hLabels;
	OneTimeOneCriterionFragmentationRawData_Get_Header_Labels(hLabels);
	PrintHeader(hLabels, out);
}

void OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Read_Header(istream& in)
{
	string buf;
	for (unsigned int i = 0; i < 4; ++i)
		getline(in, buf);
}

void OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Write_Data_void(ostream& out, int timeIndex, double timeVal)
{
	out << timeIndex << '\t';
	out << timeVal << '\t';
	out << invalidNum << '\t'; // numFragments 
	out << invalidNum << '\t'; // mean_fragmentSize 
	out << invalidNum << '\t'; // min_fragmentSize 
	out << invalidNum << '\t'; // max_fragmentSize 
	out << invalidNum << '\t'; // sdiv_fragmentSize 
	out << invalidNum << '\n'; // cov_fragmentSize 
}

void OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Get_Header_Labels(vector<HeaderLabels>& hLabels)
{
	HeaderLabels hLabel;
	hLabel.group = "index";
	hLabel.subgroup = "time";
	hLabel.textLabel = "timeIndex";
	hLabel.latexLabel = "ti";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "time";
	hLabel.latexLabel = "t";
	hLabels.push_back(hLabel);

	hLabel.group = "stat";
	hLabel.subgroup = "stat";
	hLabel.textLabel = "number_fragments";
	hLabel.latexLabel = "N";
	hLabels.push_back(hLabel);

	hLabel.group = "stat";
	hLabel.subgroup = "stat";
	hLabel.textLabel = "lmean";
	hLabel.latexLabel = "\\bar{l}";
	hLabels.push_back(hLabel);

	hLabel.group = "stat";
	hLabel.subgroup = "stat";
	hLabel.textLabel = "lmin";
	hLabel.latexLabel = "l_{\\mathrm{min}}";
	hLabels.push_back(hLabel);

	hLabel.group = "stat";
	hLabel.subgroup = "stat";
	hLabel.textLabel = "lmax";
	hLabel.latexLabel = "l_{\\mathrm{max}}";
	hLabels.push_back(hLabel);

	hLabel.group = "stat";
	hLabel.subgroup = "stat";
	hLabel.textLabel = "lstd";
	hLabel.latexLabel = "\\varsigma_l";
	hLabels.push_back(hLabel);

	hLabel.group = "stat";
	hLabel.subgroup = "stat";
	hLabel.textLabel = "lcov";
	hLabel.latexLabel = "c_l";
	hLabels.push_back(hLabel);
}

void OneTimeInterfaceFlds_FragmentationPPS2::Read_OneTimeInterfaceFlds_FragmentationPPS2(istream& in, unsigned int numInterfaces)
{
	string header1;
	getline(in, header1);
	D.resize(numInterfaces);
	maxEffDelU.resize(numInterfaces);
	for (unsigned int di = 0; di < DiM; ++di)
	{
		uL[di].resize(numInterfaces);
		uR[di].resize(numInterfaces);
		sigma[di].resize(numInterfaces);
	}
	for (unsigned int i = 0; i < numInterfaces; ++i)
	{
		in >> D[i] >> maxEffDelU[i];
		for (unsigned int di = 0; di < DiM; ++di)
			in >> uL[di][i];
		for (unsigned int di = 0; di < DiM; ++di)
			in >> uR[di][i];
		for (unsigned int di = 0; di < DiM; ++di)
			in >> sigma[di][i];
	}
}

void OneTimeInterfaceFlds_FragmentationPPS2::Empty_DUSigma(unsigned int numInterfaces)
{
	D.resize(numInterfaces);
	maxEffDelU.resize(numInterfaces);
	for (unsigned int di = 0; di < DiM; ++di)
	{
		uL[di].resize(numInterfaces);
		uR[di].resize(numInterfaces);
		sigma[di].resize(numInterfaces);
	}
	for (unsigned int i = 0; i < numInterfaces; ++i)
	{
		D[i] = invalidNum;
		maxEffDelU[i] = invalidNum;
		for (unsigned int di = 0; di < DiM; ++di)
			uL[di][i] = invalidNum;
		for (unsigned int di = 0; di < DiM; ++di)
			uR[di][i] = invalidNum;
		for (unsigned int di = 0; di < DiM; ++di)
			sigma[di][i] = invalidNum;
	}
}

bool OneTimeInterfaceFlds_FragmentationPPS2::Get_Vector_SolutionField(const OneSubdomainPPS2Data_runInfo& segmentInfo, const string& fldName, unsigned int defaultDir, vector<double>& vecVal)
{
	unsigned int sz = sigma->size();
	if (sz == 0)
		return false;
	if (fldName == "D")	{		vecVal = D;			return true;	}
	if (fldName == "sigma") { vecVal = sigma[defaultDir];			return true; }
	if (fldName == "sigma2sigmaC")
	{
		if (segmentInfo.sigmaCs.size() != sz)
			return false;
		vecVal.resize(sz);
		for (unsigned int i = 0; i < sz; ++i)
			vecVal[i] = sigma[defaultDir][i]  / segmentInfo.sigmaCs[i];
		return true;
	}
	if (fldName == "maxEffDelU") { vecVal = maxEffDelU;			return true; }
	if (fldName == "DelU") 
	{ 
		vecVal.resize(sz);
		for (unsigned int i = 0; i < sz; ++i)
			vecVal[i] = uR[defaultDir][i] - uL[defaultDir][i];
		return true; 
	}
	if (fldName == "uR") { vecVal = uR[defaultDir];			return true; }
	if (fldName == "uL") { vecVal = uL[defaultDir];			return true; }
	cout << "fldName\t" << fldName << '\n';
	THROW("Invalid field\n");
}

istream& operator>>(istream& in, OneSubdomainPPS2Data_runInfo& dat)
{
	string buf;
	in >> buf >> dat.isPeriodic;
	in >> buf >> dat.maxTime;
	in >> buf >> dat.timeStep;
	in >> buf >> dat.totalTimeSteps;
	in >> buf >> dat.numTimeStep_Interface_DSU_Fragment_Print_4PP;
	in >> buf >> dat.numTimeStep_BulkInterfacePoints_Print_4PP;
	in >> buf >> dat.numSpatialSubsegments_BulkInterfacePoints_Print_4PP;
	in >> buf >> dat.numTimeStep_InterfaceRawFinalSlnScalars_AllSpace_Print_4PP;

	in >> buf >> dat.loadTimeScale;

	in >> buf >> dat.EScale;
	in >> buf >> dat.rhoScale;
	in >> buf >> dat.dampingScale;

	in >> buf >> dat.sigmaCScale;
	in >> buf >> dat.deltaCScale;
	in >> buf >> dat.energyCScale;
	dat.ComputeAllScales();

	unsigned int sz;
	in >> buf >> sz;
	dat.lengths.resize(sz); for (unsigned int i = 0; i < sz; ++i) in >> dat.lengths[i];
	in >> buf >> sz;
	dat.deltaCs.resize(sz); for (unsigned int i = 0; i < sz; ++i) in >> dat.deltaCs[i];
	in >> buf >> sz;
	dat.sigmaCs.resize(sz); for (unsigned int i = 0; i < sz; ++i) in >> dat.sigmaCs[i];
	in >> buf >> dat.length;
	in >> buf >> dat.xm;
	in >> buf >> dat.xM;
	in >> buf >> dat.numSegments;
	in >> buf >> dat.bulk_cn_st;
	in >> buf >> dat.bulk_cn_en;
	in >> buf >> sz;
	dat.xs.resize(sz); for (unsigned int i = 0; i < sz; ++i) in >> dat.xs[i];

	dat.numInterfaces = dat.numSegments;
	dat.interface_offset = 0;
	if (!dat.isPeriodic)
	{
		++dat.numInterfaces;
		dat.interface_offset = 1;
	}

	dat.timeStep_BulkInterfacePoints_Print = dat.timeStep * dat.numTimeStep_BulkInterfacePoints_Print_4PP;
	dat.maxIndex_BulkInterfacePoints_Print = (unsigned int)floor(dat.maxTime / dat.timeStep_BulkInterfacePoints_Print + 1e-7);

	dat.timeStep_Interface_DSU_Fragment_Print = dat.timeStep * dat.numTimeStep_Interface_DSU_Fragment_Print_4PP;
	dat.maxIndex_Interface_DSU_Fragment_Print = (unsigned int)floor(dat.maxTime / dat.timeStep_Interface_DSU_Fragment_Print + 1e-7);
	return in;
}

void OneSubdomainPPS2Data_runInfo::ComputeAllScales()
{
	cScale = sqrt(EScale / rhoScale);
	ZScale = cScale * rhoScale;
	velCScale = sigmaCScale / ZScale;
	tauCScale = deltaCScale / velCScale;
	strainCScale = deltaCScale / EScale;
	psiCScale = strainCScale * sigmaCScale;
}

OneSubdomainPPS2Data_runInfo::OneSubdomainPPS2Data_runInfo()
{
	maxTime = -1.0;
	timeStep = -1.0;
}

int OneSubdomainPPS2Data_runInfo::get_Interface_DSU_Fragment_TimeIndexValue(double timeValue, double& mod_timeValue)
{
	mod_timeValue = timeValue;
	if (timeValue <= 0)
		return 0;
	if (timeValue >= maxTime)
		return totalTimeSteps;
	double index = timeValue / timeStep_Interface_DSU_Fragment_Print;
	unsigned int retInd = (unsigned int)round(index);
	retInd = MIN(retInd, maxIndex_Interface_DSU_Fragment_Print);
	mod_timeValue = retInd * timeStep_Interface_DSU_Fragment_Print;
	return retInd * numTimeStep_Interface_DSU_Fragment_Print_4PP;
}

int OneSubdomainPPS2Data_runInfo::get_Interface_BulkInterfacePoints_TimeIndexValue(double timeValue, double& mod_timeValue)
{
	mod_timeValue = timeValue;
	if (timeValue <= 0)
		return 0;
	if (timeValue - maxTime >= -1e-8 * maxTime)
		return totalTimeSteps;
	double index = timeValue / timeStep_BulkInterfacePoints_Print;
	unsigned int retInd = (unsigned int)round(index);
	retInd = MIN(retInd, maxIndex_BulkInterfacePoints_Print);
	mod_timeValue = retInd * timeStep_BulkInterfacePoints_Print;
	return retInd * numTimeStep_BulkInterfacePoints_Print_4PP;
}

bool OneSubdomainPPS2Data_runInfo::Read_OneTime_Interface_DSU_Fragment(unsigned int subdomainNo, double timeValue, OneTimeInterfaceFlds_FragmentationPPS2& fieldFragInfo, bool readFromRunOutputFolder, bool copyRunOutputFile2PPS2Folder, loadingStagesT lsi, lstClassificationT ci)
{
	double mod_timeValue;
	int timeIndex = get_Interface_DSU_Fragment_TimeIndexValue(timeValue, mod_timeValue);
	string specificName = "Interface_DSU_Fragment";
	string fileNameOutputFolder = "", fileNamePPS2 = "";

	fstream infrag;
	if (readFromRunOutputFolder)
	{
		GetSubdomainIndexed_TimeIndexed_FileName(fileNameOutputFolder, subdomainNo, timeIndex, specificName);
		infrag.open(fileNameOutputFolder.c_str(), ios::in);
		if (infrag.is_open() == false)
		{
			cout << "timeValue\t" << timeValue << '\n';
			cout << fileNameOutputFolder << '\n';
			THROW("Cannot open file\n");
		}

		if (copyRunOutputFile2PPS2Folder)
		{
			GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileNamePPS2, subdomainNo, timeIndex, specificName, lsi, ci);
			fileOperation(copyF, fileNameOutputFolder, fileNamePPS2);
		}
	}
	else
	{
		GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileNamePPS2, subdomainNo, timeIndex, specificName, lsi, ci);
		infrag.open(fileNamePPS2.c_str(), ios::in);
		if (infrag.is_open() == false)
		{
			cout << fileNamePPS2 << '\n';
			THROW("Cannot open file\n");
		}
	}
	fieldFragInfo.Read_OneTimeInterfaceFlds_FragmentationPPS2(infrag, numInterfaces);
	return true;
}

OneSubdomainPPS2Data::OneSubdomainPPS2Data()
{
}

void OneSubdomainPPS2Data::Compute_OneSubdomain_PPS2()
{

	//! 1. Computing all indices and times
	Set_All_Stage_Indices_Times();
	//! 2. Field extraction and fragmentation analysis
	ComputeAllFramentationAlaysis();
	StageRelatedData_PPS2_Write();
}

void OneSubdomainPPS2Data::StageRelatedData_PPS2_Write()
{
	string fileName, specificName;
	specificName = "TimeIndicesValues";
	GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName, loadingStagesT_none, lstClassificationT_none, FragmentationCriterionT_none, "txt", true);
	fstream out(fileName, ios::out);
	out << (int)lstClassificationT_SIZE << '\t';
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct)
		out << (lstClassificationT)ct << '\t';
	out << '\n';
	out << (int)loadingStagesT_SIZE << '\t';
	for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)
		out << (loadingStagesT)ls << '\t';
	out << '\n';
	out << "timeIndex\n";
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct)	{
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			out << diffClassificationPtr->data4Stages[ls].timeIndex_4Stage << '\t';
		out << '\n';	}
	out << "timeValue\n";
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			out << diffClassificationPtr->data4Stages[ls].timeValue_4Stage << '\t';
		out << '\n';	}
	out << "timeIndex_DSU_Fragment\n";
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			out << diffClassificationPtr->data4Stages[ls].timeIndex_DSU_Fragment4Stage << '\t';
		out << '\n';	}
	out << "timeValue_DSU_Fragment\n";
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			out << diffClassificationPtr->data4Stages[ls].timeValue_DSU_Fragment4Stage << '\t';
		out << '\n';	}
	out << "timeIndex_BulkInterfacePoints\n";
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			out << diffClassificationPtr->data4Stages[ls].timeIndex_BulkInterfacePoints4Stage << '\t';
		out << '\n';	}
	out << "timeValue_BulkInterfacePoints\n";
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			out << diffClassificationPtr->data4Stages[ls].timeValue_BulkInterfacePoints4Stage << '\t';
		out << '\n';	}
	out.close();

	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct)
	{
		if (!configPPS2->classification_Active[ct])
			continue;
		fstream* stage_frag_stat_outs[FragmentationCriterionT_SIZE];
		fstream* stage_frag_sizes_outs[FragmentationCriterionT_SIZE];
		fstream  out_integral, *out_frag_stat, *out_frag_sizes;

		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];

		specificName = "_Summary";
		GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
			loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT_none);

		out_integral.open(fileName, ios::out);
		Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_Data_Write_Header(out_integral);


		specificName = "Brittleness";
		GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
			loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT_none);
		fstream brittleout(fileName, ios::out);
		diffClassificationPtr->BrittlenessIndicators_Write(brittleout);

		if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_Stages)
		{
			for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
			{
				stage_frag_stat_outs[fc] = NULL;
				stage_frag_sizes_outs[fc] = NULL;

				if (configPPS2->factors2DecideFragmented[fc] > 0.0)
				{
					out_frag_stat = new fstream();
					stage_frag_stat_outs[fc] = out_frag_stat;

					specificName = "StatFragmentation";
					GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
						loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT(fc));

					out_frag_stat->open(fileName.c_str(), ios::out);
					OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Write_Header(*out_frag_stat);

					out_frag_sizes = new fstream();
					stage_frag_sizes_outs[fc] = out_frag_sizes;
					specificName = "SizeFragmentation";
					GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
						loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT(fc));

					out_frag_sizes->open(fileName.c_str(), ios::out);
				}
			}
		}

		///// now printing data
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)
		{
			loadingStagesT lst = (loadingStagesT)ls;
			OneTimeValuePPS2Data* ptvpd = &diffClassificationPtr->data4Stages[lst];
			int timeIndex = ptvpd->timeIndex_4Stage;
			vector<double>* vecPtr = &ptvpd->space_spacetime_integrals;// &timeSequenceSummary.data_vals[timeIndex];
			//				out_integral << timeIndex << '\t' << timeVal;
			for (unsigned int j = 0; j < vecPtr->size(); ++j)
				out_integral << (*vecPtr)[j] << '\t';
			out_integral << '\n';

			if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_Stages)
			{
				if (timeIndex > 0)
				{
					double timeVal = ptvpd->timeValue_4Stage;
					for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
					{
						if (configPPS2->factors2DecideFragmented[fc] <= 0.0)
							continue;
						OneTimeOneCriterionFragmentationRawData* fragmentationDatPtr = &ptvpd->fragmentation4Stages.fragmentationDat[fc];
						out_frag_stat = stage_frag_stat_outs[fc];
						fragmentationDatPtr->OneTimeOneCriterionFragmentationRawData_Write_Data(*out_frag_stat, timeIndex, timeVal);

						out_frag_sizes = stage_frag_sizes_outs[fc];
						fragmentationDatPtr->Fragmentation_Sizes_Write(*out_frag_sizes, timeIndex, timeVal);
					}
				}
				else
				{
					for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
					{
						if (configPPS2->factors2DecideFragmented[fc] <= 0.0)
							continue;
						OneTimeOneCriterionFragmentationRawData* fragmentationDatPtr = &ptvpd->fragmentation4Stages.fragmentationDat[fc];
						out_frag_stat = stage_frag_stat_outs[fc];
						OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Write_Data_void(*out_frag_stat, -1, -1.0);

						OneTimeOneCriterionFragmentationRawData tmpFrag;
						out_frag_sizes = stage_frag_sizes_outs[fc];
						tmpFrag.Fragmentation_Sizes_Write(*out_frag_sizes, -1, -1.0);
					}
				}
			}
		}

		for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
		{
			if (stage_frag_stat_outs[fc] != NULL)
				delete stage_frag_stat_outs[fc];
			if (stage_frag_sizes_outs[fc] != NULL)
				delete stage_frag_sizes_outs[fc];
		}
	}
}

bool OneSubdomainPPS2Data::StageRelatedData_PPS2_Read(int subdomainNoIn, DomainPostProcessS2* configPPS2In)
{
	configPPS2 = configPPS2In;
	string buf, fileName, specificName;
	subdomainNo = subdomainNoIn;
	if (segmentInfo.maxTime < 0)
	{
		bool addUnderline = true;
		specificName = "keyParameters";
		//		GetSubdomainIndexed_TimeIndexed_FileName(fileName, subdomainNo, -1, specificName, "txt", addUnderline);
		GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
			loadingStagesT_none, lstClassificationT_none, FragmentationCriterionT_none, "txt", true);
		fstream inkey(fileName.c_str(), ios::in);
		if (!inkey.is_open())
			return false;
		inkey >> segmentInfo;
	}

	specificName = "_Summary_stat";
	string fileNameIn;
	GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileNameIn, subdomainNo, -1, specificName);
	if (timeSequenceSummaryStat.DataW4LineHeader_ReadFileName(fileNameIn) == false)
		return false;

	specificName = "TimeIndicesValues";
	GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName, loadingStagesT_none, lstClassificationT_none, FragmentationCriterionT_none, "txt", true);
	fstream in(fileName, ios::in);
	getline(in, buf);
	getline(in, buf);
	in >> buf;
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			in >> diffClassificationPtr->data4Stages[ls].timeIndex_4Stage;
	}
	in >> buf;
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			in >> diffClassificationPtr->data4Stages[ls].timeValue_4Stage;
	}
	in >> buf;
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			in >> diffClassificationPtr->data4Stages[ls].timeIndex_DSU_Fragment4Stage;
	}
	in >> buf;
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			in >> diffClassificationPtr->data4Stages[ls].timeValue_DSU_Fragment4Stage;
	}
	in >> buf;
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			in >> diffClassificationPtr->data4Stages[ls].timeIndex_BulkInterfacePoints4Stage;
	}
	in >> buf;
	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct) {
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)			in >> diffClassificationPtr->data4Stages[ls].timeValue_BulkInterfacePoints4Stage;
	}
	in.close();

	for (unsigned int ct = 0; ct < lstClassificationT_SIZE; ++ct)
	{
		if (!configPPS2->classification_Active[ct])
			continue;
		fstream* stage_frag_stat_ins[FragmentationCriterionT_SIZE];
		fstream* stage_frag_sizes_ins[FragmentationCriterionT_SIZE];
		fstream  in_integral, *in_frag_stat, *in_frag_sizes;

		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];

		specificName = "_Summary";
		GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
			loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT_none);

		in_integral.open(fileName, ios::in);
		if (!in_integral.is_open())
			return false;

		for (unsigned int ri = 0; ri < 4; ++ri)
			getline(in_integral, buf);

		specificName = "Brittleness";
		GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
			loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT_none);
		fstream brittlein(fileName, ios::in);
		if (!brittlein.is_open())
			return false;
		diffClassificationPtr->BrittlenessIndicators_Read(brittlein);

		if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_Stages)
		{
			for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
			{
				stage_frag_stat_ins[fc] = NULL;
				stage_frag_sizes_ins[fc] = NULL;

				if (configPPS2->factors2DecideFragmented[fc] > 0.0)
				{
					in_frag_stat = new fstream();
					stage_frag_stat_ins[fc] = in_frag_stat;

					specificName = "StatFragmentation";
					GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
						loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT(fc));

					in_frag_stat->open(fileName.c_str(), ios::in);
					if (!in_frag_stat->is_open())
						return false;
					OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Read_Header(*in_frag_stat);

					in_frag_sizes = new fstream();
					stage_frag_sizes_ins[fc] = in_frag_sizes;
					specificName = "SizeFragmentation";
					GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
						loadingStagesT_none, lstClassificationT(ct), FragmentationCriterionT(fc));

					in_frag_sizes->open(fileName.c_str(), ios::in);
					if (!in_frag_sizes->is_open())
						return false;
				}
			}
		}

		///// now reading data
		for (unsigned int ls = 0; ls < loadingStagesT_SIZE; ++ls)
		{
			loadingStagesT lst = (loadingStagesT)ls;
			OneTimeValuePPS2Data* ptvpd = &diffClassificationPtr->data4Stages[lst];
			int timeIndex = ptvpd->timeIndex_4Stage;
			vector<double>* vecPtr = &ptvpd->space_spacetime_integrals;
			vecPtr->resize(timeSequenceSummaryStat.numFlds);
			for (unsigned int j = 0; j < vecPtr->size(); ++j)
				in_integral >> (*vecPtr)[j];

//			double timeVal = ptvpd->timeValue_4Stage;
			if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_Stages)
			{
				bool readFromRunOutputFolder = false, copyRunOutputFile2PPS2Folder = false;
				OneTimeInterfaceFlds_FragmentationPPS2* fieldFragInfoPtr = &ptvpd->fragmentation4Stages;
				if (timeIndex > 0)
					segmentInfo.Read_OneTime_Interface_DSU_Fragment(subdomainNo, ptvpd->timeValue_4Stage, *fieldFragInfoPtr, readFromRunOutputFolder, copyRunOutputFile2PPS2Folder, lst, (lstClassificationT)ct);
				else
					fieldFragInfoPtr->Empty_DUSigma(segmentInfo.numInterfaces);

				for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
				{
					if (configPPS2->factors2DecideFragmented[fc] <= 0.0)
						continue;
					OneTimeOneCriterionFragmentationRawData* fragmentationDatPtr = &ptvpd->fragmentation4Stages.fragmentationDat[fc];
					in_frag_stat = stage_frag_stat_ins[fc];
					*in_frag_stat >> buf >> buf;
					fragmentationDatPtr->OneTimeOneCriterionFragmentationRawData_Read_Data(*in_frag_stat);
					in_frag_sizes = stage_frag_sizes_ins[fc];
					fragmentationDatPtr->Fragmentation_Sizes_Read(*in_frag_sizes);
				}
			}
		}

		for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
		{
			if (stage_frag_stat_ins[fc] != NULL)
				delete stage_frag_stat_ins[fc];
			if (stage_frag_sizes_ins[fc] != NULL)
				delete stage_frag_sizes_ins[fc];
		}
	}
	return true;
}
bool OneSubdomainPPS2Data::Compute_OneTimeValuePPS2Data_4ActualTime(PPS2_TimeStamp timeStamp, OneTimeValuePPS2Data &oneTimeSlns, bool extract_space_spacetimeIntegrals, bool extract_DSUFields, bool compute_FragmentationStats)
{
	//! A. Calculating the time value and various indices
	if (oneTimeSlns.timeIndex_4Stage < 0)
	{
		double tol = segmentInfo.maxTime * 1e-6;
		if ((timeStamp.actualTime_index == timeIndexType_none) || (timeStamp.actualTime_val > 0))
			oneTimeSlns.timeValue_4Stage = timeStamp.actualTime_val;
		else
		{
			if (timeStamp.actualTime_indexType == ti_AllTimeSteps)
				oneTimeSlns.timeValue_4Stage = segmentInfo.timeStep * timeStamp.actualTime_index;
			else if (timeStamp.actualTime_indexType == ti_DSU_Fragment)
				oneTimeSlns.timeValue_4Stage = segmentInfo.timeStep_Interface_DSU_Fragment_Print * timeStamp.actualTime_index;
			else if (timeStamp.actualTime_indexType == ti_BulkInterfacePoints)
				oneTimeSlns.timeValue_4Stage = segmentInfo.timeStep_BulkInterfacePoints_Print * timeStamp.actualTime_index;
			else if (timeStamp.actualTime_indexType == ti_RawFinalSln_Scalars)
				oneTimeSlns.timeValue_4Stage = (segmentInfo.maxTime / segmentInfo.numTimeStep_InterfaceRawFinalSlnScalars_AllSpace_Print_4PP) * timeStamp.actualTime_index;
			else
			{
				cout << "timeStamp.actualTime_index\t" << timeStamp.actualTime_index << '\n';
				THROW("Invalid timeStamp.actualTime_index\n");
			}
		}
		//	if (segmentInfo.maxTime - oneTimeSlns.timeValue_4Stage > tol)
		//		return false;
		if (oneTimeSlns.timeValue_4Stage < -tol)
			return false;
		oneTimeSlns.timeValue_4Stage = MIN(MAX(oneTimeSlns.timeValue_4Stage, 0.0), segmentInfo.maxTime);
		oneTimeSlns.timeIndex_4Stage = (int)floor(oneTimeSlns.timeValue_4Stage / segmentInfo.timeStep + 1e-7);
		oneTimeSlns.timeIndex_4Stage = MIN(MAX(oneTimeSlns.timeIndex_4Stage, 0), (int)segmentInfo.totalTimeSteps);
		oneTimeSlns.SetAllTimeValsIndices_FromActualTime(segmentInfo);
	}
	if (extract_space_spacetimeIntegrals && (oneTimeSlns.space_spacetime_integrals.size() == 0))
	{
		Ensure_openning_complete_time_history_space_spacetime_integrals_summary_file();
		oneTimeSlns.space_spacetime_integrals = timeSequenceSummary.data_vals[oneTimeSlns.timeIndex_4Stage];
	}
	if (compute_FragmentationStats)
		extract_DSUFields = true;
	if ((extract_DSUFields) && (oneTimeSlns.fragmentation4Stages.sigma->size() == 0))
	{
		bool readFromRunOutputFolder = true, copyRunOutputFile2PPS2Folder = false;
		loadingStagesT lsi = loadingStagesT_none;
		lstClassificationT ci = lstClassificationT_none;
		if (segmentInfo.Read_OneTime_Interface_DSU_Fragment(subdomainNo, oneTimeSlns.timeValue_4Stage, oneTimeSlns.fragmentation4Stages, readFromRunOutputFolder, copyRunOutputFile2PPS2Folder, lsi, ci) == false)
			return false;

	}
	if (compute_FragmentationStats)
		Do_AllFragmentationAnalysis(oneTimeSlns.fragmentation4Stages);
	return true;
}

bool OneSubdomainPPS2Data::Get_Scalar_Or_Vector_Output(PPS2_dataPointer & datPointer, double & scalarVal, vector<double>& vecVal, bool & outputIsScalar, OneTimeValuePPS2Data & modifiableOneTimeSlns, double temporalFieldTimeStepValOrNumSteps)
{
	outputIsScalar = true;
	/// first check if brittleness indicators are needed
	if (datPointer.brittlenessRatioT != BrittlenessIndicatorRatioType_none)
	{
		BrittlnessIndicatorFieldType bf = datPointer.brittlenessFieldT;
		if (bf == BrittlnessIndicatorFieldType_none)
			bf = birt_eps;
		lstClassificationT lct = datPointer.timeStamp.lct;
		if (lct == lstClassificationT_none)
			lct = configPPS2->default_ClassificationType;
		scalarVal = diffClassifications[lct].brittlenessIndicators[bf][datPointer.brittlenessRatioT];
		return (scalarVal >= 0.0);
	}

	// first check if this is a time history data
	if (datPointer.stat_sv_type == sso_all_field_time_val)
	{
		outputIsScalar = false;
		Ensure_openning_complete_time_history_space_spacetime_integrals_summary_file();
		int fldIndex = configPPS2->get_space_spacetime_integral_index_from_name(datPointer.fldName);
		if (fldIndex < 0)
			return false;
		double delT = temporalFieldTimeStepValOrNumSteps;
		if (delT < 0)
			delT = -segmentInfo.maxTime / delT;
		int timeStep_Factor, max_timeStep;
		if (delT < segmentInfo.timeStep)
		{
			cout << "delT\t" << delT << '\n';
			cout << "segmentInfo.timeStep\t" << segmentInfo.timeStep << '\n';
			THROW("Increase delT = temporalFieldTimeStepValOrNumSteps\n");
		}
		timeStep_Factor = (int)floor(delT / segmentInfo.timeStep + 1e-7);
		max_timeStep = segmentInfo.totalTimeSteps / timeStep_Factor;
		unsigned int sz = max_timeStep + 1;
		if ((max_timeStep *  timeStep_Factor) >= (int)timeSequenceSummary.numDataPoints)
		{
			cout << "max_timeStep\t" << max_timeStep << '\n';
			cout << "timeStep_Factor\t" << timeStep_Factor << '\n';
			cout << "timeSequenceSummary.numDataPoints\t" << timeSequenceSummary.numDataPoints << '\n';
			THROW("invalid size\n");
		}
		vecVal.resize(sz);
		for (unsigned int ti = 0; ti < sz; ++ti)
			vecVal[ti] = timeSequenceSummary.data_vals[ti * timeStep_Factor][fldIndex];
		return true;
	}
	scalarVal = invalidNum;
	// so no brittleness indicator from this point-on 
	// next check if the data corresponds to the stats from space/spacetime integrals (i.e. timeSequenceSummaryStat)
	bool ptNeedsFragmentation = (datPointer.fragmentationCriterion != FragmentationCriterionT_none);
	outputIsScalar = (datPointer.stat_sv_type != sso_all_field_val);
	bool isScalarStatistics = ((datPointer.stat_sv_type != sso_none) && outputIsScalar);

	if (isScalarStatistics && !ptNeedsFragmentation) // stat type, non-fragmentation, and scalar -> timeSequenceSummaryStat
	{
		int fldIndex = configPPS2->get_space_spacetime_integral_index_from_name(datPointer.fldName);
		if (fldIndex < 0)
			return false;
		scalarVal = timeSequenceSummaryStat.data_vals[(int)datPointer.stat_sv_type][fldIndex];
		return true;
	}

	/////////////////////////////////////////////////////
	//// from this point on, we are dealing with a specific OneTimeValuePPS2Data
	// checking to see if we are dealing with classification + stage pair OR a true given time value
	OneTimeValuePPS2Data* oneTimeAllSolutions;
//	bool classification_stage_type = true;
	// classification / stage type
	if ((int)datPointer.timeStamp.timeStampType >= (int)lst_lin)
	{
		lstClassificationT lct = datPointer.timeStamp.lct;
		if (lct == lstClassificationT_none)
			lct = configPPS2->default_ClassificationType;
		oneTimeAllSolutions = &diffClassifications[lct].data4Stages[datPointer.timeStamp.timeStampType];
	}
	else
	{
//		bool compute_FragmentationStats = (configPPS2->bActualTime_compute_FragmentationStats || ptNeedsFragmentation);
//		bool extract_DSUFields = (configPPS2->bActualTime_extract_DSUFields || !outputIsScalar || compute_FragmentationStats);
//		bool extract_space_spacetimeIntegrals = (configPPS2->bActualTime_extract_space_spacetimeIntegrals || (outputIsScalar && !isScalarStatistics));

		bool compute_FragmentationStats = ptNeedsFragmentation;
		bool extract_DSUFields = (!outputIsScalar || compute_FragmentationStats);
		bool extract_space_spacetimeIntegrals = (outputIsScalar && !isScalarStatistics);
		if (Compute_OneTimeValuePPS2Data_4ActualTime(datPointer.timeStamp, modifiableOneTimeSlns, extract_space_spacetimeIntegrals, extract_DSUFields, compute_FragmentationStats) == false)
			return false;
		oneTimeAllSolutions = &modifiableOneTimeSlns;
//		classification_stage_type = false;
	}
	/// now there are 3 possibilitis for this one time solution
	// A. Fragmentation related
	if (ptNeedsFragmentation)
	{
		OneTimeOneCriterionFragmentationRawData* fragDag = &oneTimeAllSolutions->fragmentation4Stages.fragmentationDat[datPointer.fragmentationCriterion];
		return fragDag->Get_Scalar_Or_Vector_Output(datPointer.stat_sv_type, scalarVal, vecVal, outputIsScalar);
	}
	// B. field related values (D, U, DelU, ...)
	if (!outputIsScalar) // dealing with D, U, DelU, ... fields
		return oneTimeAllSolutions->fragmentation4Stages.Get_Vector_SolutionField(segmentInfo, datPointer.fldName, configPPS2->defaultDir, vecVal);

	// C. dealing with space/spacetime interal values
	int fldIndex = configPPS2->get_space_spacetime_integral_index_from_name(datPointer.fldName);
	if (fldIndex < 0)
		return false;
	scalarVal = oneTimeAllSolutions->space_spacetime_integrals[fldIndex];
	return true;
}

void OneSubdomainPPS2Data::Set_All_Stage_Indices_Times()
{
	/// strain, stresses
	Set_Time_Stage_Indices_Times();
	Set_EpsSigma_Stage_Indices_Times(ct_epssig_bc);
	Set_EpsSigma_Stage_Indices_Times(ct_epssig_bulk);
	Set_Damage_Stage_Indices_Times();
	Set_PhysInterfaceDiss_Stage_Indices_Times();
	Set_U_Stage_Indices_Times();

	for (int ct = 0; ct < lstClassificationT_SIZE; ++ct)
	{
		OneClassificationPPS2Data* diffClassificationPtr = &diffClassifications[ct];
#if 0
		diffClassificationPtr->subdomainNo = subdomainNo;
		diffClassificationPtr->ct = (lstClassificationT)ct;
#endif
		for (int lst = 0; lst < loadingStagesT_SIZE; ++lst)
		{
#if 0
			diffClassificationPtr->data4Stages[lst].subdomainNo = subdomainNo;
			diffClassificationPtr->data4Stages[lst].ct = (lstClassificationT)ct;
			diffClassificationPtr->data4Stages[lst].lst = (loadingStagesT)lst;
#endif
			diffClassificationPtr->data4Stages[lst].SetAllTimeValsIndices_FromActualTime(segmentInfo);
			int timeIndex = diffClassificationPtr->data4Stages[lst].timeIndex_4Stage;
			if (timeIndex >= 0)
				diffClassificationPtr->data4Stages[lst].space_spacetime_integrals = timeSequenceSummary.data_vals[timeIndex];
			else
			{
				diffClassificationPtr->data4Stages[lst].space_spacetime_integrals.resize(timeSequenceSummary.numFlds);
				fill(diffClassificationPtr->data4Stages[lst].space_spacetime_integrals.begin(), diffClassificationPtr->data4Stages[lst].space_spacetime_integrals.end(), invalidNum);
			}
		}
	}
}

void OneSubdomainPPS2Data::Compute_BrittlenessIndices(lstClassificationT ci, int fldStrn, int fldEne)
{
	if (fldStrn < 0)
	{
		if (configPPS2->classification_Active[ct_epssig_bulk] == true)
		{
			fldStrn = configPPS2->sind_eps_bulk;
			fldEne = configPPS2->sind_phi_tot_bulk;
		}
		else
		{
			fldStrn = configPPS2->sind_eps_bc;
			fldEne = configPPS2->sind_phi_tot_bc;
		}
	}

	OneClassificationPPS2Data* diffClassification = &diffClassifications[ci];
	int lin_index = diffClassification->data4Stages[lst_lin].timeIndex_4Stage;
	int max_index = diffClassification->data4Stages[lst_max].timeIndex_4Stage;
	int final_index = diffClassification->data4Stages[lst_final].timeIndex_4Stage;

	double max_time, max_strn, max_ene;
	bool hasMax = (max_index > 0); // equal is not acceptable
	if (hasMax)
	{
		max_time = diffClassification->data4Stages[lst_max].timeValue_4Stage;
		max_strn = timeSequenceSummary.data_vals[max_index][fldStrn];
		max_ene = timeSequenceSummary.data_vals[max_index][fldEne];
	}
	double final_time, final_strn, final_ene;
	bool hasFinal = (final_index > 0); // equal is not acceptable
	if (hasFinal)
	{
		final_time = diffClassification->data4Stages[lst_final].timeValue_4Stage;
		final_strn = timeSequenceSummary.data_vals[final_index][fldStrn];
		final_ene = timeSequenceSummary.data_vals[final_index][fldEne];
	}

	BrittlenessIndicatorRatioType br;
	if (lin_index >= 0) // i/m & i/f
	{
		double lin_time = diffClassification->data4Stages[lst_lin].timeValue_4Stage;
		double lin_strn = timeSequenceSummary.data_vals[lin_index][fldStrn];
		double lin_ene = timeSequenceSummary.data_vals[lin_index][fldEne];

		if (hasMax)
		{
			br = birt_i2m;
			diffClassification->brittlenessIndicators[birt_time][br] = lin_time / max_time;
			diffClassification->brittlenessIndicators[birt_eps][br] = lin_strn / max_strn;
			diffClassification->brittlenessIndicators[birt_psi][br] = lin_ene / max_ene;
		}
		if (hasFinal)
		{
			br = birt_i2f;
			diffClassification->brittlenessIndicators[birt_time][br] = lin_time / final_time;
			diffClassification->brittlenessIndicators[birt_eps][br] = lin_strn / final_strn;
			diffClassification->brittlenessIndicators[birt_psi][br] = lin_ene / final_ene;
		}
	}
	if (hasMax && hasFinal)
	{
		br = birt_m2f;
		diffClassification->brittlenessIndicators[birt_time][br] = max_time / final_time;
		diffClassification->brittlenessIndicators[birt_eps][br] = max_strn / final_strn;
		diffClassification->brittlenessIndicators[birt_psi][br] = max_ene / final_ene;
	}
}

void OneSubdomainPPS2Data::Set_Time_Stage_Indices_Times()
{
	lstClassificationT ci = ct_time;
	double actualTimesProvided_4TimeStage[loadingStagesT_SIZE];
	for (unsigned int i = 0; i < loadingStagesT_SIZE; ++i)
	{
		loadingStagesT lsi = (loadingStagesT)i;
		actualTimesProvided_4TimeStage[i] = configPPS2->actualTimesProvided_4TimeStage[i];
		if ((lsi == lst_final) && (actualTimesProvided_4TimeStage[i] < 0))
			actualTimesProvided_4TimeStage[i] = timeSequenceSummary.data_vals[timeSequenceSummary.numDataPoints - 1][configPPS2->sind_time] * (1.0 - 1e-9);
		if (actualTimesProvided_4TimeStage[i] > 0)
		{
			unsigned int closestSortingIndex;
			double closestSortingValOut;
			timeSequenceSummary.GetIndexClosestSortingVal(configPPS2->sind_time, actualTimesProvided_4TimeStage[i], closestSortingIndex, closestSortingValOut);

			diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = closestSortingIndex;
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = closestSortingValOut;
		}
	}

	loadingStagesT lsi = lst_lin;
	if (actualTimesProvided_4TimeStage[lsi] <= 0)
	{
		// no time is provided for lin stage ->
		// for lin position can use the time when input energy does not increase any more
		// computing final stage for input energy

//		double classification_tol4FinalValCheck = configPPS2->tol4FinalValCheck_energyInput;
		double classification_tol4FinalValCheck = configPPS2->classification_tol4LinValChecks[ci];
		bool classification_forwardFinalValueSearch = configPPS2->classification_forwardFinalValueSearches[ci];
		bool classification_usePositive4MaxValue = configPPS2->classification_usePositive4MaxValues[ci];
		unsigned int fldNo_sind_EneInp = configPPS2->sind_EneInp;

		bool forceTerminalvalue2ProvidedValue = false;
		double providedTerminalValue = 0.0;
		bool  bValuesPositive;
		int index_maxAbsVal, index_finalVal;
		double maxAbsVal, finalVal;
		timeSequenceSummary.GetMaxAbs_FinalValue(timeSequenceSummaryStat, fldNo_sind_EneInp, classification_usePositive4MaxValue, forceTerminalvalue2ProvidedValue, providedTerminalValue, classification_forwardFinalValueSearch, classification_tol4FinalValCheck,
			bValuesPositive, index_maxAbsVal, maxAbsVal, index_finalVal, finalVal);

		diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_finalVal;
		if (index_finalVal >= 0)
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_finalVal][configPPS2->sind_time];
		else
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
	}
	Compute_BrittlenessIndices(ct_time);
}

void OneSubdomainPPS2Data::Set_EpsSigma_Stage_Indices_Times(lstClassificationT ci)
{
	if (configPPS2->classification_Active[ci] == false)
	{
		for (unsigned int lsi = 0; lsi < loadingStagesT_SIZE; ++lsi)
		{
			diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = -1;
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
		}
		return;
	}
	double classification_tol4FinalValCheck = configPPS2->classification_tol4FinalValChecks[ci];
	bool classification_forwardFinalValueSearch = configPPS2->classification_forwardFinalValueSearches[ci];
	bool classification_usePositive4MaxValue = configPPS2->classification_usePositive4MaxValues[ci];
	unsigned int fldNo_sigma = configPPS2->sind_sig_bc, fldno_psi_lost = configPPS2->sind_phi_diss_lost_bc;
	
	if (ci == ct_epssig_bulk)
	{
		fldNo_sigma = configPPS2->sind_sig_bulk;
		fldno_psi_lost = configPPS2->sind_phi_diss_lost_bulk;
	}
	// computing max and final stages
	bool forceTerminalvalue2ProvidedValue = true;
	bool  bValuesPositive;
	int index_maxAbsVal, index_finalVal;
	double maxAbsVal, finalVal;
	timeSequenceSummary.GetMaxAbs_FinalValue(timeSequenceSummaryStat, fldNo_sigma, classification_usePositive4MaxValue, forceTerminalvalue2ProvidedValue, 0.0, classification_forwardFinalValueSearch, classification_tol4FinalValCheck,
		bValuesPositive, index_maxAbsVal, maxAbsVal, index_finalVal, finalVal);

	loadingStagesT lsi = lst_max;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage= index_maxAbsVal;
	if (index_maxAbsVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_maxAbsVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	lsi = lst_final;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_finalVal;
	if (index_finalVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_finalVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	// lin limit
	double classification_tol4LinValCheck = configPPS2->classification_tol4LinValChecks[ci];
	double tol = classification_tol4LinValCheck;
	if (tol > 0)
	{
		// absolute tolerance for stress, need to turn this to a relative one based on the maximum stress observed
		if (index_maxAbsVal < 0)
		{
			THROW("Cannot find max/min stress");
		}
		tol = -classification_tol4LinValCheck / maxAbsVal;
	}

	int index_zero_psi_lost;
	timeSequenceSummary.GetLinValue(timeSequenceSummaryStat, fldno_psi_lost, tol, index_zero_psi_lost);
	lsi = lst_lin;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_zero_psi_lost;
	if (index_zero_psi_lost >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_zero_psi_lost][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
	int fldStrn = configPPS2->sind_eps_bc, fldEne = configPPS2->sind_phi_tot_bc;
	if (ci == ct_epssig_bulk)
	{
		fldStrn = configPPS2->sind_eps_bulk;
		fldEne = configPPS2->sind_phi_tot_bulk;
	}
	Compute_BrittlenessIndices(ci, fldStrn, fldEne);
}

void OneSubdomainPPS2Data::Set_Damage_Stage_Indices_Times()
{
	lstClassificationT ci = ct_Damage;
	if (configPPS2->classification_Active[ci] == false)
	{
		for (unsigned int lsi = 0; lsi < loadingStagesT_SIZE; ++lsi)
		{
			diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = -1;
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
		}
		return;
	}
	double classification_tol4FinalValCheck = configPPS2->classification_tol4FinalValChecks[ci];
	bool classification_forwardFinalValueSearch = configPPS2->classification_forwardFinalValueSearches[ci];
	unsigned int fldNo_Dmax = configPPS2->sind_Dmax, fldno_DsrcMax = configPPS2->sind_DsrcMax;

	// computing final stage
	bool forceTerminalvalue2ProvidedValue = true;
	double providedTerminalValue = 1.0;
	bool  bValuesPositive;
	int index_maxAbsVal, index_finalVal;
	double maxAbsVal, finalVal;
	timeSequenceSummary.GetFirstLocationGreaterThanValue(timeSequenceSummaryStat, fldNo_Dmax, providedTerminalValue, classification_tol4FinalValCheck, index_finalVal);

	loadingStagesT lsi = lst_final;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_finalVal;
	if (index_finalVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_finalVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	// max is not computed from this, but lin limit is computed here
	int index_zero_Dmax;
	double classification_tol4LinValCheck = configPPS2->classification_tol4LinValChecks[ci];
	timeSequenceSummary.GetLinValue(timeSequenceSummaryStat, fldNo_Dmax, classification_tol4LinValCheck, index_zero_Dmax);
	lsi = lst_lin;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_zero_Dmax;
	if (index_zero_Dmax >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_zero_Dmax][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;


	///// for max, damage source is used
	forceTerminalvalue2ProvidedValue = true;
	providedTerminalValue = 0.0;
	if (classification_tol4FinalValCheck > 0)
		classification_tol4FinalValCheck = -classification_tol4FinalValCheck; // want to use a relative tolerance for damage source (although we don't use the final value anyways)
	timeSequenceSummary.GetMaxAbs_FinalValue(timeSequenceSummaryStat, fldno_DsrcMax, true, forceTerminalvalue2ProvidedValue, providedTerminalValue, classification_forwardFinalValueSearch, classification_tol4FinalValCheck,
		bValuesPositive, index_maxAbsVal, maxAbsVal, index_finalVal, finalVal);

	lsi = lst_max;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_maxAbsVal;
	if (index_maxAbsVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_maxAbsVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	Compute_BrittlenessIndices(ci);
}

void OneSubdomainPPS2Data::Set_PhysInterfaceDiss_Stage_Indices_Times()
{
	lstClassificationT ci = ct_interfaceDissEne;
	if (configPPS2->classification_Active[ci] == false)
	{
		for (unsigned int lsi = 0; lsi < loadingStagesT_SIZE; ++lsi)
		{
			diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = -1;
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
		}
		return;
	}
	double classification_tol4FinalValCheck = configPPS2->classification_tol4FinalValChecks[ci];
	bool classification_forwardFinalValueSearch = configPPS2->classification_forwardFinalValueSearches[ci];
	unsigned int fldNo_diss_lost = configPPS2->sind_diss_interface_lost, fldno_dissPower_lost = configPPS2->sind_dissPower_interface_lost;

	// computing final stage
	bool forceTerminalvalue2ProvidedValue = false;
	double providedTerminalValue = 0.0;
	bool  bValuesPositive;
	int index_maxAbsVal, index_finalVal;
	double maxAbsVal, finalVal;
	timeSequenceSummary.GetMaxAbs_FinalValue(timeSequenceSummaryStat, fldNo_diss_lost, true, forceTerminalvalue2ProvidedValue, providedTerminalValue, classification_forwardFinalValueSearch, classification_tol4FinalValCheck,
		bValuesPositive, index_maxAbsVal, maxAbsVal, index_finalVal, finalVal);

	loadingStagesT lsi = lst_final;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_finalVal;
	if (index_finalVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_finalVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
	bool classification_usePositive4MaxValue = configPPS2->classification_usePositive4MaxValues[ci];

	// max is not computed from this, but lin limit is computed here
	int index_zero;
	double classification_tol4LinValCheck = configPPS2->classification_tol4LinValChecks[ci];
	timeSequenceSummary.GetLinValue(timeSequenceSummaryStat, fldNo_diss_lost, classification_tol4LinValCheck, index_zero);
	lsi = lst_lin;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_zero;
	if (index_zero >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_zero][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	///// for max, interface dissippation power is used
	forceTerminalvalue2ProvidedValue = true;
	providedTerminalValue = 0.0;
	if (classification_tol4FinalValCheck > 0)
		classification_tol4FinalValCheck = -0.01; // a relative value is used for terminal, even though terminal value is not needed
	timeSequenceSummary.GetMaxAbs_FinalValue(timeSequenceSummaryStat, fldno_dissPower_lost, classification_usePositive4MaxValue, forceTerminalvalue2ProvidedValue, providedTerminalValue, classification_forwardFinalValueSearch, classification_tol4FinalValCheck,
		bValuesPositive, index_maxAbsVal, maxAbsVal, index_finalVal, finalVal);

	lsi = lst_max;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_maxAbsVal;
	if (index_maxAbsVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_maxAbsVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	Compute_BrittlenessIndices(ci);
}

void OneSubdomainPPS2Data::Set_U_Stage_Indices_Times()
{
	lstClassificationT ci = ct_U;
	if (configPPS2->classification_Active[ci] == false)
	{
		for (unsigned int lsi = 0; lsi < loadingStagesT_SIZE; ++lsi)
		{
			diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = -1;
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
		}
		return;
	}
	double classification_tol4FinalValCheck = configPPS2->classification_tol4FinalValChecks[ci];
	bool classification_forwardFinalValueSearch = configPPS2->classification_forwardFinalValueSearches[ci];
	unsigned int fldNo_U = configPPS2->sind_U;

	// computing max and final stage
	bool forceTerminalvalue2ProvidedValue = true;
	double providedTerminalValue = 0.0;
	bool  bValuesPositive;
	int index_maxAbsVal, index_finalVal;
	double maxAbsVal, finalVal;
	timeSequenceSummary.GetMaxAbs_FinalValue(timeSequenceSummaryStat, fldNo_U, true, forceTerminalvalue2ProvidedValue, providedTerminalValue, classification_forwardFinalValueSearch, classification_tol4FinalValCheck,
		bValuesPositive, index_maxAbsVal, maxAbsVal, index_finalVal, finalVal);

	loadingStagesT lsi = lst_max;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_maxAbsVal;
	if (index_maxAbsVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_maxAbsVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	lsi = lst_final;
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = index_finalVal;
	if (index_finalVal >= 0)
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = timeSequenceSummary.data_vals[index_finalVal][configPPS2->sind_time];
	else
		diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;

	// lin is taken from bc/bulk strain
	lsi = lst_lin;
	lstClassificationT ci_eps = ct_epssig_bulk;
	if (configPPS2->classification_Active[ci_eps] == false)
	{
		ci_eps = ct_epssig_bc;
		if (configPPS2->classification_Active[ci_eps] == false)
		{
			diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = -1;
			diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = -1.0;
			return;
		}
	}
	diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage = diffClassifications[ci_eps].data4Stages[lsi].timeIndex_4Stage;
	diffClassifications[ci].data4Stages[lsi].timeValue_4Stage = diffClassifications[ci_eps].data4Stages[lsi].timeValue_4Stage;

	Compute_BrittlenessIndices(ci);
}

void OneSubdomainPPS2Data::ComputeAllFramentationAlaysis()
{
	bool readFromRunOutputFolder = true, copyRunOutputFile2PPS2Folder = true;
	if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_Stages)
	{
		for (unsigned int ci = 0; ci < lstClassificationT_SIZE; ++ci)
		{
			if (configPPS2->classification_Active[ci] == false)
				continue;
			OneClassificationPPS2Data* diffClassification = &diffClassifications[ci];
			for (unsigned int lsi = 0; lsi < loadingStagesT_SIZE; ++lsi)
			{
				int timeIndex = diffClassifications[ci].data4Stages[lsi].timeIndex_4Stage;
				if (timeIndex < 0)
					continue;
				double timeValue = diffClassifications[ci].data4Stages[lsi].timeValue_4Stage;
				OneTimeInterfaceFlds_FragmentationPPS2* fieldFragInfoPtr = &diffClassification->data4Stages[lsi].fragmentation4Stages;
				segmentInfo.Read_OneTime_Interface_DSU_Fragment(subdomainNo, timeValue, *fieldFragInfoPtr, readFromRunOutputFolder, copyRunOutputFile2PPS2Folder, (loadingStagesT)lsi, (lstClassificationT)ci);
				Do_AllFragmentationAnalysis(*fieldFragInfoPtr);
			}
		}
	}
	if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_AllTimes)
	{
		readFromRunOutputFolder = true, copyRunOutputFile2PPS2Folder = false;
		fstream* outPtrs[FragmentationCriterionT_SIZE];
		fstream* outFragSizesPtrs[FragmentationCriterionT_SIZE];

		for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
		{
			outPtrs[fc] = NULL;
			outFragSizesPtrs[fc] = NULL;
			if (configPPS2->factors2DecideFragmented[fc] > 0.0)
			{
				FragmentationCriterionT fct = (FragmentationCriterionT)fc;
				string fileName;
				string specificName = "StatFragmentation";
				GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
					loadingStagesT_none, lstClassificationT_none, fct, "txt");

				outPtrs[fc] = new fstream();
				outPtrs[fc]->open(fileName.c_str(), ios::out);
				OneTimeOneCriterionFragmentationRawData::OneTimeOneCriterionFragmentationRawData_Write_Header(*(outPtrs[fc]));

				if (configPPS2->doFragmentation_DSigmaFieldExtraction_4_AllTimes_IncludedDetailedFragmentSizes)
				{
					string specificName = "SizeFragmentation";
					GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomainNo, -1, specificName,
						loadingStagesT_none, lstClassificationT_none, fct, "txt");
					outFragSizesPtrs[fc] = new fstream();
					outFragSizesPtrs[fc]->open(fileName.c_str(), ios::out);
				}
			}
		}
		/// now the printing time
		for (unsigned int ti = 0; ti < segmentInfo.numTimeStep_Interface_DSU_Fragment_Print_4PP; ++ti)
		{
			double timeValue = ti * segmentInfo.timeStep_Interface_DSU_Fragment_Print;
			unsigned int timeIndex = (unsigned int)floor(timeValue / segmentInfo.timeStep + 1e-7);
			OneTimeInterfaceFlds_FragmentationPPS2 fieldFragInfo;
			segmentInfo.Read_OneTime_Interface_DSU_Fragment(subdomainNo, timeValue, fieldFragInfo, readFromRunOutputFolder, copyRunOutputFile2PPS2Folder, loadingStagesT_none, lstClassificationT_none);
			Do_AllFragmentationAnalysis(fieldFragInfo);
			for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
			{
				if (outPtrs[fc] != NULL)
					fieldFragInfo.fragmentationDat[fc].OneTimeOneCriterionFragmentationRawData_Write_Data(*(outPtrs[fc]), timeIndex, timeValue);
				if (outFragSizesPtrs[fc] != NULL)
					fieldFragInfo.fragmentationDat[fc].Fragmentation_Sizes_Write(*(outFragSizesPtrs[fc]), timeIndex, timeValue);
			}
		}
		for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
		{
			if (outPtrs[fc] != NULL)
				delete outPtrs[fc];
			if (outFragSizesPtrs[fc] != NULL)
				delete outFragSizesPtrs[fc];
		}
	}
}

void OneSubdomainPPS2Data::Do_AllFragmentationAnalysis(OneTimeInterfaceFlds_FragmentationPPS2& fragDat)
{
	for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
	{
		FragmentationCriterionT fragCriterion = (FragmentationCriterionT)fc;
		double factor2DecideFragmented = configPPS2->factors2DecideFragmented[fc];
		if ((factor2DecideFragmented < 0) || (fragDat.fragmentationDat[fc].numFragments > 0))
			continue;
		Compute_FragmentationSizeData_Step1_DetermineBreakingPoints(fragCriterion, factor2DecideFragmented, fragDat);
		Compute_FragmentationSizeData_Step2_DetermineFragments(fragCriterion, fragDat);
		fragDat.fragmentationDat[fc].ComputeFragmentStatsFromFragmentSizes();
	}
}

void OneSubdomainPPS2Data::Compute_FragmentationSizeData_Step1_DetermineBreakingPoints(FragmentationCriterionT fragCriterion, double factor2DecideFragmented, OneTimeInterfaceFlds_FragmentationPPS2& fragDat)
{
	OneTimeOneCriterionFragmentationRawData* fragmentationDatPtr = &fragDat.fragmentationDat[fragCriterion];
	fragmentationDatPtr->fragmented_interface_indices.clear();

	unsigned int defaultDir = configPPS2->defaultDir;
	// this is not fully correct if we have symmetry BC, but ignired now
	unsigned int first_interior_interface = segmentInfo.interface_offset, end_interior_interface = segmentInfo.numInterfaces- 1 - segmentInfo.interface_offset;
	if (fragCriterion == fct_Damage)
	{
		for (unsigned int ii = first_interior_interface; ii <= end_interior_interface; ++ii)
			if (fragDat.D[ii] >= factor2DecideFragmented)
				fragmentationDatPtr->fragmented_interface_indices.push_back(ii);
	}
	else if (fragCriterion == fct_maxEffDelU)
	{
		for (unsigned int ii = first_interior_interface; ii <= end_interior_interface; ++ii)
			if (fragDat.maxEffDelU[ii] >= factor2DecideFragmented * segmentInfo.deltaCs[ii])
				fragmentationDatPtr->fragmented_interface_indices.push_back(ii);
	}
	else if (fragCriterion == fct_DelU)
	{
		double delu;
		double factor2DecideFragmentedMaxDelU = configPPS2->factors2DecideFragmented[fct_maxEffDelU];

		for (unsigned int ii = first_interior_interface; ii <= end_interior_interface; ++ii)
		{
			delu = fragDat.uR[defaultDir][ii] - fragDat.uL[defaultDir][ii];
			if (defaultDir > 0)
				delu = fabs(delu);
			if ((fragDat.maxEffDelU[ii] >= factor2DecideFragmentedMaxDelU * segmentInfo.deltaCs[ii]) &&
				(delu >= factor2DecideFragmented * segmentInfo.deltaCs[ii]))
				fragmentationDatPtr->fragmented_interface_indices.push_back(ii);
		}
	}
}

void OneSubdomainPPS2Data::Compute_FragmentationSizeData_Step2_DetermineFragments(FragmentationCriterionT fragCriterion, OneTimeInterfaceFlds_FragmentationPPS2& fragDat)
{
	OneTimeOneCriterionFragmentationRawData* fragmentationDatPtr = &fragDat.fragmentationDat[fragCriterion];
	fragmentationDatPtr->fragment_lengths.clear();
	vector<unsigned int> fragmented_interface_indices = fragmentationDatPtr->fragmented_interface_indices;
	if (!segmentInfo.isPeriodic)
	{
		// numbering convention 
		// Interface 0 - bulk 0 - interface 1 - bulk 1 - ..... - interface(n-1) bulk(n - 1) interface(n)
		// n = number of segments
		fragmented_interface_indices.push_back(segmentInfo.numInterfaces - 1);
		fragmentationDatPtr->numFragments = fragmented_interface_indices.size();
		if (fragmentationDatPtr->numFragments == 1)
			fragmentationDatPtr->fragment_lengths.push_back(segmentInfo.length);
		else
		{
			unsigned int st = 0, en;
			double segLength;
			for (unsigned int i = 0; i < fragmentationDatPtr->numFragments; ++i)
			{
				segLength = 0;
				en = fragmented_interface_indices[i];
				for (unsigned int j = st; j < en; ++j)
					segLength += segmentInfo.lengths[j];
				fragmentationDatPtr->fragment_lengths.push_back(segLength);
				st = en;
			}
		}
	}
	else
	{
		// numbering convention 
		// [Interface n-1 = periodic] - bulk 0 - interface 0 - bulk 1 - interface 1..... - interface(n-2) bulk(n - 1) interface(n-1)
		// n = number of segments
		// so interfce numbers are subtracted by one relative to previous case
		// periodic, this time 
		fragmentationDatPtr->numFragments = fragmented_interface_indices.size();
		if (fragmentationDatPtr->numFragments == 0)
		{
			fragmentationDatPtr->numFragments = 1;
			fragmentationDatPtr->fragment_lengths.push_back(segmentInfo.length);
		}
		else
		{
			unsigned int lastInterfaceIndex = segmentInfo.numInterfaces - 1;
			bool periodInterfaceBroken = Find(fragmented_interface_indices, lastInterfaceIndex);
			if (!periodInterfaceBroken)
				fragmented_interface_indices.push_back(lastInterfaceIndex);

			int st = -1, en;
			double segLength;
			for (unsigned int i = 0; i < fragmented_interface_indices.size(); ++i)
			{
				segLength = 0;
				en = fragmented_interface_indices[i];
				for (int j = st + 1; j <= en; ++j)
					segLength += segmentInfo.lengths[j];
				fragmentationDatPtr->fragment_lengths.push_back(segLength);
				st = en;
			}
			if (!periodInterfaceBroken)
			{
				// the last segment length needs to be added to the first one and its value removed from the link as this link is not broken
				double lastPartialSegLength = fragmentationDatPtr->fragment_lengths[fragmentationDatPtr->fragment_lengths.size() - 1];
				fragmentationDatPtr->fragment_lengths[0] += lastPartialSegLength;
				// removing that value
				fragmentationDatPtr->fragment_lengths.pop_back();
			}
		}
	}
	// some basic checks to ensure values are correct
	if (fragmentationDatPtr->fragment_lengths.size() != fragmentationDatPtr->numFragments)
	{
		cout << "fragmentationDatPtr->numFragments\t" << fragmentationDatPtr->numFragments << '\n';
		cout << "fragmentationDatPtr->fragment_lengths.size()\t" << fragmentationDatPtr->fragment_lengths.size() << '\n';
		THROW("Invalid size\n");
	}
	double totLength = 0;
	for (unsigned int i = 0; i < fragmentationDatPtr->numFragments; ++i)
		totLength += fragmentationDatPtr->fragment_lengths[i];
	double rat = computeRatio(totLength, segmentInfo.length);
	if (fabs(rat - 1.0) > 1e-4)
	{
		cout << "totLength\t" << totLength << '\n';
		cout << "segmentInfo.length\t" << segmentInfo.length << '\n';
		THROW("sizes don't match\n");
	}
}

void OneSubdomainPPS2Data::Ensure_openning_complete_time_history_space_spacetime_integrals_summary_file()
{
	if (timeSequenceSummary.numDataPoints > 0)
		return;
	string fileNameIn;
	string specificName = "_Summary";
	GetSubdomainIndexed_TimeIndexed_FileName(fileNameIn, subdomainNo, -1, specificName);
	timeSequenceSummary.DataW4LineHeader_ReadFileName(fileNameIn);
}

DomainPostProcessS2::DomainPostProcessS2()
{
	mainSubdomainNo = -1;
	numSubdomains = 0;
	defaultDir = 0;
#if DiM1
	defaultDirStr = "";
#else
	defaultDirStr = "0";
#endif
	for (unsigned int lsi = 0; lsi < loadingStagesT_SIZE; ++lsi)
		actualTimesProvided_4TimeStage[lsi] = -1.0;

	double tolLin = 0.00001, tolFinal = 0.001;
	for (unsigned int ci = 0; ci < lstClassificationT_SIZE; ++ci)
	{
		classification_Active[ci] = true;
		classification_tol4LinValChecks[ci] = -tolLin;
		classification_tol4FinalValChecks[ci] = -tolFinal;
		classification_forwardFinalValueSearches[ci] = true;
		classification_usePositive4MaxValues[ci] = true;
	}
	classification_tol4LinValChecks[ct_Damage] = tolLin;
	classification_tol4FinalValChecks[ct_Damage] = tolFinal;

//	tol4FinalValCheck_energyInput = -tolFinal;

	for (unsigned int fc = 0; fc < FragmentationCriterionT_SIZE; ++fc)
		factors2DecideFragmented[fc] = 0.40; // 0.90;

	doFragmentation_DSigmaFieldExtraction_4_Stages = true;
	doFragmentation_DSigmaFieldExtraction_4_AllTimes = true;
	doFragmentation_DSigmaFieldExtraction_4_AllTimes_IncludedDetailedFragmentSizes = true;

	doCalculateStages = true;
//	doCalculateActualTimeIntervals = false;
	stageSolutionsExist = false;
	redoPPS2Calculations = false;

	/// for time stage
//	bActualTime_extract_space_spacetimeIntegrals = true;
//	bActualTime_extract_DSUFields = true;
//	bActualTime_compute_FragmentationStats = true;
}

void DomainPostProcessS2::DomainPostProcessS2_Read_WO_Initialization(istream& in)
{
	double tmpd;
	string buf;
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
		if (buf == "doCalculateStages")
		{
			READ_NBOOL(in, buf, doCalculateStages);
		}
		else if (buf == "redoPPS2Calculations")
		{
			READ_NBOOL(in, buf, redoPPS2Calculations);
		}
		else if (buf == "doFragmentation_DSigmaFieldExtraction_4_Stages")
		{
			READ_NBOOL(in, buf, doFragmentation_DSigmaFieldExtraction_4_Stages);
		}
		else if (buf == "doFragmentation_DSigmaFieldExtraction_4_Stages")
		{
			READ_NBOOL(in, buf, doFragmentation_DSigmaFieldExtraction_4_Stages);
		}
		else if (buf == "doFragmentation_DSigmaFieldExtraction_4_AllTimes")
		{
			READ_NBOOL(in, buf, doFragmentation_DSigmaFieldExtraction_4_AllTimes);
		}
		else if (buf == "doFragmentation_DSigmaFieldExtraction_4_AllTimes_IncludedDetailedFragmentSizes")
		{
			READ_NBOOL(in, buf, doFragmentation_DSigmaFieldExtraction_4_AllTimes_IncludedDetailedFragmentSizes);
		}
#if 0
		else if (buf == "doCalculateActualTimeIntervals")
		{
			READ_NBOOL(in, buf, doCalculateActualTimeIntervals);
		}
		else if (buf == "bActualTime_extract_space_spacetimeIntegrals")
		{
			READ_NBOOL(in, buf, bActualTime_extract_space_spacetimeIntegrals);
		}
		else if (buf == "bActualTime_extract_DSUFields")
		{
			READ_NBOOL(in, buf, bActualTime_extract_DSUFields);
		}
		else if (buf == "bActualTime_compute_FragmentationStats")
		{
			READ_NBOOL(in, buf, bActualTime_compute_FragmentationStats);
		}
#endif		
		else if (buf == "factors2DecideFragmented")
		{
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("istream should start with {");
			}
			READ_NSTRING(in, buf, buf);
			while (buf != "}")
			{
				FragmentationCriterionT fc;
				name2Type(buf, fc);
				READ_NDOUBLE(in, buf, tmpd);
				factors2DecideFragmented[fc] = tmpd;
				READ_NSTRING(in, buf, buf);
			}
		}
		else if (buf == "actualTimesProvided_4TimeStage")
		{
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("istream should start with {");
			}
			READ_NSTRING(in, buf, buf);
			while (buf != "}")
			{
				loadingStagesT ls;
				name2Type(buf, ls);
				READ_NDOUBLE(in, buf, tmpd);
				actualTimesProvided_4TimeStage[ls] = tmpd;
				READ_NSTRING(in, buf, buf);
			}
		}
#if 0
		else if (buf == "actualTimesProvided_4TimeStage")
		{
			READ_NDOUBLE(in, buf, tol4FinalValCheck_energyInput);
		}
#endif
		else if (buf == "classification")
		{
			READ_NSTRING(in, buf, buf);
			if (buf == "classificationType")
				READ_NSTRING(in, buf, buf);
			lstClassificationT lct;
			name2Type(buf, lct);
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("istream should start with {");
			}
			READ_NSTRING(in, buf, buf);
			while (buf != "}")
			{
				if (buf == "Active")
				{
					READ_NBOOL(in, buf, classification_Active[lct]);
				}
				else if (buf == "tol4LinValChecks")
				{
					READ_NDOUBLE(in, buf, classification_tol4LinValChecks[lct]);
				}
				else if (buf == "tol4FinalValChecks")
				{
					READ_NDOUBLE(in, buf, classification_tol4FinalValChecks[lct]);
				}
				else if (buf == "forwardFinalValueSearches")
				{
					READ_NBOOL(in, buf, classification_forwardFinalValueSearches[lct]);
				}
				else if (buf == "usePositive4MaxValues")
				{
					READ_NBOOL(in, buf, classification_usePositive4MaxValues[lct]);
				}
				else
				{
					cout << "Invalid buf\t" << buf << '\n';
					THROW("exit\n");
				}
				READ_NSTRING(in, buf, buf);
			}
		}
		READ_NSTRING(in, buf, buf);
	}
}

void DomainPostProcessS2::Main_DomainPostProcessS2()
{
	g_SL_desc_data.Read_tdLoad();
	Initialize_DomainPostProcessS2();
	Read_SummaryData();
	if (!stageSolutionsExist)
		Compute_AllSubdomain_PPS2();

	string fn;
	string specificName = "_PPS2_success";
	GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fn, -1, -1, specificName);
	fstream outcheck(fn.c_str(), ios::out);
	outcheck << "1";
}

bool DomainPostProcessS2::Get_Scalar_Or_Vector_Output(PPS2_dataPointer& datPointer, double& scalarVal, vector<double>& vecVal, bool& outputIsScalar, OneTimeValuePPS2Data& modifiableOneTimeSlns, double temporalFieldTimeStepValOrNumSteps, int spatialFieldResolutionCorrector)
{
	vecVal.clear();
	scalarVal = invalidNum;

	bool retVal = true;
	int subDomainNo = datPointer.subdomainNo;
	if (subDomainNo < 0)
		subDomainNo = mainSubdomainNo;
	bool isInifield = false;
	if (datPointer.fldName == "field_finalResolution")
	{
		isInifield = true;
		int index = 0;
		if (datPointer.timeStamp.actualTime_index >= 0)
			index = MIN(datPointer.timeStamp.actualTime_index, (int)fields_finalResolution.size() - 1);
		outputIsScalar = false;
		vecVal = fields_finalResolution[index];
	}
	else
	{
		retVal = onesubdomainPPS2[subDomainNo].Get_Scalar_Or_Vector_Output(datPointer, scalarVal, vecVal, outputIsScalar, modifiableOneTimeSlns, temporalFieldTimeStepValOrNumSteps);
		if (retVal == false)
		{
			// if dealing with a spatial field, can size it correctly and set values to invalidNum
			if ((datPointer.stat_sv_type == sso_all_field_val) && (datPointer.fragmentationCriterion == FragmentationCriterionT_none)) 
			{
				vecVal.resize(onesubdomainPPS2[subDomainNo].segmentInfo.numInterfaces);
				fill(vecVal.begin(), vecVal.end(), invalidNum);
			}
			else
				return false;
		}
		else
		{
			if (outputIsScalar)
			{
				if (!Is_Valid_Num(scalarVal))
					retVal = false;
			}
			else
			{
				unsigned int sz = vecVal.size();
				for (unsigned int i = 0; i < sz; ++i)
				{
					if (!Is_Valid_Num(vecVal[i]))
					{
						retVal = false;
						break;
					}
				}
			}
		}
	}
	/////  Now post-processings:
	// Step A: normalizing, Step B: Operation (log, ...), Step C: (changing the resolution)
	// C is only for fields
	if (outputIsScalar)
	{
		if (retVal)
		{
			// Step A
			if (datPointer.normalizationMode != RunNormalizationQuantT_none)
				scalarVal /= normalizations[(int)datPointer.normalizationMode];
			// Step B
			if (datPointer.oprType != vpo_none)
				scalarVal = get_valPOperT(scalarVal, datPointer.oprType);
		}
		return retVal;
	}
	// dealing with field values
	unsigned int szField = vecVal.size();

	// Step A
	if (retVal)
	{
		if (datPointer.normalizationMode != RunNormalizationQuantT_none)
		{
			for (unsigned int i = 0; i < szField; ++i)
				vecVal[i] /= normalizations[(int)datPointer.normalizationMode];
		}
		// Step B
		if (datPointer.oprType != vpo_none)
		{
			for (unsigned int i = 0; i < szField; ++i)
				vecVal[i] = get_valPOperT(vecVal[i], datPointer.oprType);
		}
	}
	// Step C: Changing the resolution, only applies for spatial fields (not fragmentation sizes)
	bool ptNeedsFragmentation = (datPointer.fragmentationCriterion != FragmentationCriterionT_none);
	if ((spatialFieldResolutionCorrector == 0) || (ptNeedsFragmentation) || (datPointer.stat_sv_type == sso_all_field_time_val))
		return retVal;
	// dealing with a spatial field whose resolution needs to be changed
	// only simple cases are considered
	OneSubdomainPPS2Data_runInfo* segInfo = &onesubdomainPPS2[subDomainNo].segmentInfo;
	unsigned int numSegments = segInfo->numSegments;
	int numSegmentsNew = spatialFieldResolutionCorrector;
	unsigned int ratio; //numSegments_over_numSegmentsNew;
	if (numSegmentsNew < 0) // i.e. it's given as a reduction factor
	{
		ratio = -numSegmentsNew;
		//	ratio = numSegments / numSegmentsNew;	-> numSegmentsNew = numSegments / ratio 
		if (numSegments % ratio != 0)
		{
			cout << "numSegments\t" << numSegments << '\n';
			cout << "ratio = numSegments / numSegmentsNew\t" << ratio << '\n';
			cout << "numSegments % ratio\t" << numSegments % ratio << '\n';
			THROW("cannot do the simple resolution reduction\n");
		}
		numSegmentsNew = numSegments / ratio;
	}
	else
	{
		if (numSegments % numSegmentsNew != 0)
		{
			cout << "numSegments\t" << numSegments << '\n';
			cout << "numSegmentsNew\t" << numSegmentsNew << '\n';
			cout << "numSegments % numSegmentsNew\t" << numSegments % numSegmentsNew << '\n';
			THROW("cannot do the simple resultion reduction\n");
		}
		ratio = numSegments / numSegmentsNew;
	}
	// offset is used as the values are given at vertices for D,S,U, etc.
	unsigned int offset = segInfo->interface_offset;
	if (isInifield)
	{
		if (vecVal.size() == segInfo->numSegments)
			offset = 0;
		else if (vecVal.size() != segInfo->numInterfaces)
		{
			cout << "vecVal is not of the right size\t" << vecVal.size();
			cout << "segInfo->numSegments\t" << segInfo->numSegments << '\n';
			cout << "segInfo->numInterfaces\t" << segInfo->numInterfaces << '\n';
			cout << "segInfo->interface_offset" << segInfo->interface_offset << '\n';
			THROW("Invalid size\n");
		}
	}
	vector<double> vecBK = vecVal;
	vecVal.resize(numSegmentsNew + offset);
	for (int i = 0; i < numSegmentsNew; ++i)
	{
		int ii = i * ratio;
		vecVal[i] = vecBK[ii];
	}
	if (offset == 1)
		vecVal[numSegmentsNew] = vecBK[numSegments];
	return retVal;
}

int DomainPostProcessS2::get_space_spacetime_integral_index_from_name(const string& fldName)
{
	if (fldName == "time")
		return sind_time;
	if (fldName == "phi")
		return sind_phi;
	if (fldName == "K")
		return sind_K;
	if (fldName == "U")
		return sind_U;
	if (fldName == "phys_diss_tot")
		return sind_diss_tot;
	if (fldName == "phys_diss_lost")
		return sind_diss_lost;
	if (fldName == "phys_diss_interface_lost")
		return sind_diss_interface_lost;

	if (fldName == "eps_bc")
		return sind_eps_bc;
	if (fldName == "sigma_bc")
		return sind_sig_bc;
	if (fldName == "eps_bulk")
		return sind_eps_bulk;
	if (fldName == "sigma_bulk")
		return sind_sig_bulk;
	if (fldName == "eps_bulk_intfc")
		return sind_eps_bulk_intfc;
	if (fldName == "linMomentum")
		return sind_P;

	if (fldName == "phys_diss_recoverable")
		return sind_diss_recov;
	if (fldName == "energy_phys_diss_2_input")
		return sind_energy_phys_diss_2_input;
	if (fldName == "energy_phys_diss_lost_2_input")
		return sind_energy_phys_diss_lost_2_input;

	if (fldName == "input_energy")
		return sind_EneInp;
	if (fldName == "eneIDiss")
		return sind_eneIDiss;
	if (fldName == "energyIDiss_Recoverable_2_input")
		return sind_energyIDiss_Recoverable_2_input;
	if (fldName == "energy_numerical_diss_2_input")
		return sind_energy_numerical_diss_2_input;
	if (fldName == "powerIDiss")
		return sind_dissPower_interface_lost;

	if (fldName == "energy_eps_total_bc")
		return sind_phi_tot_bc;
	if (fldName == "energy_eps_diss_bc")
		return sind_phi_diss_lost_bc;
	if (fldName == "energy_eps_recoverable_bc")
		return sind_phi_recov_bc;

	if (fldName == "energy_eps_total_bulk")
		return sind_phi_tot_bulk;
	if (fldName == "energy_eps_diss_bulk")
		return sind_phi_diss_lost_bulk;
	if (fldName == "energy_eps_recoverable_bulk")
		return sind_phi_recov_bulk;

	if (fldName == "mean_interface_damage")
		return sind_Dbar;
	if (fldName == "max_interface_damage")
		return sind_Dmax;
	if (fldName == "min_interface_damage")
		return sind_Dmin;
	if (fldName == "max_interface_damage_source")
		return sind_DsrcMax;

	if (fldName == "energy_BC")
		return sind_EneBC;
	if (fldName == "energy_L")
		return sind_EneL;
	if (fldName == "energy_R")
		return sind_EneR;
	if (fldName == "phi0")
		return sind_phi0;

	DataW4LineHeader* summaryStat = &onesubdomainPPS2[mainSubdomainNo].timeSequenceSummaryStat;
	vector<string> namesNonLatex = summaryStat->headerStrs3;
	return Find(namesNonLatex, fldName);
}

void DomainPostProcessS2::Initialize_DomainPostProcessS2()
{
	stageSolutionsExist = false;
//	if ((doCalculateActualTimeIntervals == false) && (redoPPS2Calculations == false))// no need to load full data
	if (redoPPS2Calculations == false)// no need to load full data
	{
		if (numSubdomains == 0)
		{
			subDomainActive4PPS2.clear();
			bool exists = true;
			unsigned int subdomain_number = 0;
			while (exists)
			{
				string fileName;
				string specificName = "_Summary_stat";
				GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileName, subdomain_number, -1, specificName);
				fstream in(fileName.c_str(), ios::in);
				exists = (in.is_open());
				if (exists)
				{
					subDomainActive4PPS2.push_back(exists);
					++subdomain_number;
				}
			}
			numSubdomains = subDomainActive4PPS2.size();
		}
		stageSolutionsExist = (numSubdomains > 0);
		if (stageSolutionsExist)
		{
			string fn;
			string specificName = "_PPS2_success";
			GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fn, -1, -1, specificName);
			fstream incheck(fn.c_str(), ios::in);
			if (!incheck.is_open())
				stageSolutionsExist = false;
			int success = -1;
			incheck >> success;
			if (success != 1)
				stageSolutionsExist = false;
		}
	}
	if (numSubdomains == 0)
	{
		subDomainActive4PPS2.clear();
		bool exists = true;
		unsigned int subdomain_number = 0;
		while (exists)
		{
			string fileName;
			string specificName = "_Summary";
			GetSubdomainIndexed_TimeIndexed_FileName(fileName, subdomain_number, -1, specificName);
			fstream in(fileName.c_str(), ios::in);
			exists = (in.is_open());
			if (exists)
			{
				subDomainActive4PPS2.push_back(exists);
				++subdomain_number;
			}
		}
		numSubdomains = subDomainActive4PPS2.size();
	}
}

bool DomainPostProcessS2::Read_SummaryData()
{
	string fileName;
	string specificName = "mainSubdomainNo";
	bool addUnderline = true;
	GetSubdomainIndexed_TimeIndexed_FileName(fileName, -1, -1, specificName, "txt", addUnderline);
	fstream in(fileName.c_str(), ios::in);
	if (in.is_open())
	{
		in >> mainSubdomainNo;
		in.close();
	}
	onesubdomainPPS2.resize(numSubdomains);
	if (!stageSolutionsExist)
	{
		for (unsigned int si = 0; si < numSubdomains; ++si)
		{
			if (!subDomainActive4PPS2[si])
				continue;
			if (mainSubdomainNo < 0)
				mainSubdomainNo = si;

			string fileNameIn;
			string specificName = "_Summary";
			GetSubdomainIndexed_TimeIndexed_FileName(fileNameIn, si, -1, specificName);

			specificName = "_Summary_stat";
			string fileNameOut;
			GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileNameOut, si, -1, specificName);
			//			loadingStagesT_none, lstClassificationT_none, FragmentationCriterionT_none, "txt", false);

			onesubdomainPPS2[si].configPPS2 = this;
			onesubdomainPPS2[si].subdomainNo = si;
			onesubdomainPPS2[si].timeSequenceSummary.Initialize_DataW4LineHeader(fileNameIn, onesubdomainPPS2[si].timeSequenceSummaryStat, fileNameOut);

			bool addUnderline = true;
			specificName = "keyParameters";
			GetSubdomainIndexed_TimeIndexed_FileName(fileNameIn, si, -1, specificName, "txt", addUnderline);
			fstream inkey(fileNameIn.c_str(), ios::in);
			if (!inkey.is_open())
			{
				cout << fileNameIn << '\n';
				THROW("Cannot open file\n");
			}
			inkey >> onesubdomainPPS2[si].segmentInfo;

			GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(fileNameOut, si, -1, specificName,
				loadingStagesT_none, lstClassificationT_none, FragmentationCriterionT_none, "txt", true);
			fileOperation(copyF, fileNameIn, fileNameOut);
		}
	}
	else
	{
		for (unsigned int si = 0; si < numSubdomains; ++si)
		{
			if (!subDomainActive4PPS2[si])
				continue;
			if (mainSubdomainNo < 0)
				mainSubdomainNo = si;
			onesubdomainPPS2[si].StageRelatedData_PPS2_Read(si, this);
		}
	}
	if (mainSubdomainNo < 0)
		return false;

	default_ClassificationType = lstClassificationT_none;
	for (unsigned int ct = (int)ct_epssig_bulk; ct < lstClassificationT_SIZE; ++ct)
	{
		if (classification_Active[ct])
		{
			default_ClassificationType = (lstClassificationT)ct;
			break;
		}
	}
	if ((default_ClassificationType == lstClassificationT_none) && (classification_Active[ct_time] == true))
		default_ClassificationType = ct_time;

	/// setting the infices
	DataW4LineHeader* summaryStat = &onesubdomainPPS2[mainSubdomainNo].timeSequenceSummaryStat;
	vector<string> namesNonLatex = summaryStat->headerStrs3;
	
	string tmps = "time";								sind_time = Find(namesNonLatex, tmps);
	tmps = "phi";										sind_phi = Find(namesNonLatex, tmps);
	tmps = "K";											sind_K = Find(namesNonLatex, tmps);
	tmps = "U";											sind_U = Find(namesNonLatex, tmps);
	tmps = "input_energy";								sind_EneInp = Find(namesNonLatex, tmps);
	tmps = "energy_BC";									sind_EneBC = Find(namesNonLatex, tmps);
	tmps = "energy_L";									sind_EneL = Find(namesNonLatex, tmps);
	tmps = "energy_R";									sind_EneR = Find(namesNonLatex, tmps);
	tmps = "phi0";										sind_phi0 = Find(namesNonLatex, tmps);
	tmps = "phys_diss_tot";								sind_diss_tot = Find(namesNonLatex, tmps);
	tmps = "phys_diss_lost";							sind_diss_lost = Find(namesNonLatex, tmps);
	tmps = "eneIDiss";									sind_eneIDiss = Find(namesNonLatex, tmps);
	tmps = "phys_diss_interface_lost";					sind_diss_interface_lost = Find(namesNonLatex, tmps);
	tmps = "phys_diss_recoverable";						sind_diss_recov = Find(namesNonLatex, tmps);
	tmps = "energy_phys_diss_2_input";					sind_energy_phys_diss_2_input = Find(namesNonLatex, tmps);
	tmps = "energy_phys_diss_lost_2_input";				sind_energy_phys_diss_lost_2_input = Find(namesNonLatex, tmps);
	tmps = "energyIDiss_Recoverable_2_input";			sind_energyIDiss_Recoverable_2_input = Find(namesNonLatex, tmps);
	tmps = "energy_numerical_diss_2_input";				sind_energy_numerical_diss_2_input = Find(namesNonLatex, tmps);

	tmps = "powerIDiss";								sind_dissPower_interface_lost = Find(namesNonLatex, tmps);


	tmps = "eps_bc" + defaultDirStr;					sind_eps_bc = Find(namesNonLatex, tmps);
	tmps = "eps_bulk_intfc" + defaultDirStr;			sind_eps_bulk_intfc = Find(namesNonLatex, tmps);

	tmps = "sigma_bc" + defaultDirStr;					sind_sig_bc = Find(namesNonLatex, tmps);
	tmps = "energy_eps_total_bc";						sind_phi_tot_bc = Find(namesNonLatex, tmps);
	tmps = "energy_eps_diss_bc";						sind_phi_diss_lost_bc = Find(namesNonLatex, tmps);
	tmps = "energy_eps_recoverable_bc";					sind_phi_recov_bc = Find(namesNonLatex, tmps);

	tmps = "eps_bulk" + defaultDirStr;					sind_eps_bulk = Find(namesNonLatex, tmps);
	tmps = "sigma_bulk" + defaultDirStr;				sind_sig_bulk = Find(namesNonLatex, tmps);
	tmps = "energy_eps_total_bulk";						sind_phi_tot_bulk = Find(namesNonLatex, tmps);
	tmps = "energy_eps_diss_bulk";						sind_phi_diss_lost_bulk = Find(namesNonLatex, tmps);
	tmps = "energy_eps_recoverable_bulk";				sind_phi_recov_bulk = Find(namesNonLatex, tmps);

	tmps = "mean_interface_damage";						sind_Dbar = Find(namesNonLatex, tmps);
	tmps = "max_interface_damage";						sind_Dmax = Find(namesNonLatex, tmps);
	tmps = "min_interface_damage";						sind_Dmin = Find(namesNonLatex, tmps);
	tmps = "max_interface_damage_source";				sind_DsrcMax = Find(namesNonLatex, tmps);
	tmps = "linMomentum";								sind_P = Find(namesNonLatex, tmps);

	/// setting normalizations
	for (unsigned int nt = 0; nt < RunNormalizationQuantT_SIZE; ++nt)
		normalizations[nt] = 1.0;

	OneSubdomainPPS2Data_runInfo* segmentInfoPtr = &onesubdomainPPS2[mainSubdomainNo].segmentInfo;
	if (summaryStat->numDataPoints > 0)
	{
		for (unsigned int nt = 0; nt < RunNormalizationQuantT_SIZE; ++nt)
			normalizations[nt] = 0.0;
		for (unsigned int si = 0; si < numSubdomains; ++si)
		{
			if (!subDomainActive4PPS2[si])
				continue;
			DataW4LineHeader* summaryStati = &onesubdomainPPS2[si].timeSequenceSummaryStat;
			normalizations[ene_KInitial] += summaryStati->data_vals[sso_valStart][sind_K];
			normalizations[ene_UInitial] += summaryStati->data_vals[sso_valStart][sind_U];
			normalizations[ene_PhiInitial] += summaryStati->data_vals[sso_valStart][sind_phi];

			normalizations[ene_KFinal] += summaryStati->data_vals[sso_valEnd][sind_K];
			normalizations[ene_UFinal] += summaryStati->data_vals[sso_valEnd][sind_U];
			normalizations[ene_PhiFinal] += summaryStati->data_vals[sso_valEnd][sind_phi];

			normalizations[eneBCFinal] += summaryStati->data_vals[sso_valEnd][sind_EneBC];
			normalizations[eneInputFinal] += summaryStati->data_vals[sso_valEnd][sind_EneInp];
		}
	}
	if (g_SL_desc_data.tdLoad_inputEnergy > 0.0)
	{
		normalizations[ene_PhiInitial] = g_SL_desc_data.tdLoad_inputEnergy;
		normalizations[eneBCFinal] = g_SL_desc_data.tdLoad_inputEnergy;
	}

	normalizations[rnq_loadTimeScale] = segmentInfoPtr->loadTimeScale;
	normalizations[rnq_EScale] = segmentInfoPtr->EScale;
	normalizations[rnq_rhoScale] = segmentInfoPtr->rhoScale;
	normalizations[rnq_dampingScale] = segmentInfoPtr->dampingScale;
	normalizations[rnq_sigmaCScale] = segmentInfoPtr->sigmaCScale;
	normalizations[rnq_deltaCScale] = segmentInfoPtr->deltaCScale;
	normalizations[rnq_cScale] = segmentInfoPtr->cScale;
	normalizations[rnq_tauCScale] = segmentInfoPtr->tauCScale;
	normalizations[rnq_strainCScale] = segmentInfoPtr->strainCScale;
	normalizations[rnq_velCScale] = segmentInfoPtr->velCScale;
	normalizations[rnq_psiCScale] = segmentInfoPtr->psiCScale;

	/// reading random fields
	int cntr = 0;
	bool contLoop = true;
	while (contLoop)
	{
		string str, irfileSource, irfileDestination, frfileDestination;
		toString(cntr, str);
		irfileDestination = g_prefileNamePPS2 + "field" + str + "_initialResolution.txt";
		frfileDestination = g_prefileNamePPS2 + "field" + str + "_finalResolution.txt";
		fstream in(irfileDestination.c_str(), ios::in);
		if (!in.is_open())
		{
			irfileSource = g_prefileName + "field" + str + "_initialResolution.txt";
			in.open(irfileSource.c_str(), ios::in);
			if (!in.is_open())
				contLoop = false;
			else
			{
				fileOperation(copyF, irfileSource, irfileDestination);
				string tmp;
				tmp = g_prefileName + "field" + str + "_finalResolution.txt";
				fileOperation(copyF, tmp, frfileDestination);
			}
			if (contLoop)
				in.open(irfileDestination.c_str(), ios::in);
		}
		if (contLoop)
		{
			vector<double> tmpV;
			double tmp;
			string tmps;
			in >> tmps;
			bool valid = (fromString(tmps, tmp) && (!in.eof()));
			while (valid)
			{ 
				tmpV.push_back(tmp);			
				in >> tmps;
				valid = (fromString(tmps, tmp) && (!in.eof()));
			}
			fields_initialResolution.push_back(tmpV);
			in.close(); tmpV.clear();
			in.open(frfileDestination.c_str(), ios::in);
			in >> tmps;
			valid = (fromString(tmps, tmp) && (!in.eof()));
			while (valid)
			{
				tmpV.push_back(tmp);
				in >> tmps;
				valid = (fromString(tmps, tmp) && (!in.eof()));
			}
			fields_finalResolution.push_back(tmpV);
			++cntr;
		}
	}
	return true;
}

void DomainPostProcessS2::Compute_AllSubdomain_PPS2()
{
	for (unsigned int si = 0; si < numSubdomains; ++si)
	{
		if (!subDomainActive4PPS2[si])
			continue;
		onesubdomainPPS2[si].Compute_OneSubdomain_PPS2();
	}
}

void GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(string& fileName, unsigned int subdomainNo, int timeIndex, const string& specificName,
	loadingStagesT lsi, lstClassificationT ci, FragmentationCriterionT fc, string ext, bool addUnderline)
{
	string tmps, posStr = "", nameBase, subDomainStr;
	int tmpi;
	if (lsi != loadingStagesT_none)
	{
		posStr = "__tStage_";
		tmpi = (int)lsi;
		toString(tmpi, tmps);
		posStr += tmps;
		tmps = getName(lsi);
		posStr += "_";
		posStr += tmps;
	}
	else
	{
		if (timeIndex >= 0)
		{
			toString(timeIndex, posStr);
			posStr = "_tPos_" + posStr;
		}
		else
			posStr = "__tAll_";
	}
	toString(subdomainNo, subDomainStr);
	string root, sd;
	root = g_prefileNamePPS2;
	if (!addUnderline)
		sd = "sd_";
	else
		sd = "_sd_";
	if (subdomainNo == -1)
	{

		subDomainStr = "";
		sd = "_";

	}

	string clsStr = "";
	if (ci != lstClassificationT_none)
	{
		clsStr = "_Class_";
		tmpi = (int)ci;
		toString(tmpi, tmps);
		clsStr += tmps;
		tmps = getName(ci);
		clsStr += "_";
		clsStr += tmps;
	}
	string fcStr = "";
	if (fc != FragmentationCriterionT_none)
	{
		fcStr = "_FragCrn_";
		tmpi = (int)fc;
		toString(tmpi, tmps);
		fcStr += tmps;
		tmps = getName(fc);
		fcStr += "_";
		fcStr += tmps;
	}
	nameBase = root + sd + subDomainStr + posStr + clsStr + fcStr + + "_" + specificName;
//	nameBase = root + sd + subDomainStr + "_" + specificName + posStr;
	fileName = nameBase + "." + ext;
}

void MAIN_DomainPostProcessS2(string fileName)
{
	fstream in(fileName.c_str(), ios::in);
	DomainPostProcessS2 dpps2;
	dpps2.DomainPostProcessS2_Read_WO_Initialization(in);
	dpps2.Main_DomainPostProcessS2();
}
