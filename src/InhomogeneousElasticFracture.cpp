#include "InhomogeneousElasticFracture.h"
#include "globalMacros.h"

string getName(elasticInhomogT dat)
{
	if (dat == eih_undecided)
		return "undecided";
	if (dat == eih_Crho_const)
		return "Crho_const";
	if (dat == eih_rho_const)
		return "rho_const";
	if (dat == eih_C_const)
		return "C_const";
	if (dat == eih_c_const)
		return "c_const";
	if (dat == eih_Z_const)
		return "Z_const";
	if (dat == eih_nonconst)
		return "nonconst";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, elasticInhomogT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= elasticInhomogT_SIZE)
			THROW("too large of a number\n");
		typeVal = (elasticInhomogT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < elasticInhomogT_SIZE; ++i)
	{
		typeVal = (elasticInhomogT)i; // casting integer to elasticInhomogT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading elasticInhomogT\n");
}

//operator for output
ostream& operator<<(ostream& out, elasticInhomogT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, elasticInhomogT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

////////////////////////////////////
string getName(fractureInhomogT dat)
{
	if (dat == fih_undecided)
		return "undecided";
	if (dat == fih_sigmadelta_const)
		return "sigmadelta_const";
	if (dat == fih_delta_const)
		return "delta_const";
	if (dat == fih_sigma_const)
		return "sigma_const";
	if (dat == fih_G_const)
		return "G_const";
	if (dat == fih_nonconst)
		return "nonconst";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, fractureInhomogT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= fractureInhomogT_SIZE)
			THROW("too large of a number\n");
		typeVal = (fractureInhomogT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < fractureInhomogT_SIZE; ++i)
	{
		typeVal = (fractureInhomogT)i; // casting integer to fractureInhomogT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading fractureInhomogT\n");
}

//operator for output
ostream& operator<<(ostream& out, fractureInhomogT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, fractureInhomogT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

ostream& operator<<(ostream& out, const ElasticFractureInhomogFactors& dat)
{
	out << dat.CFactor << '\t';
	out << dat.rhoFactor << '\t';
	out << dat.dampingFactor << '\t';
	out << dat.sigmaFactor << '\t';
	out << dat.deltaFactor << '\t';
	out << dat.iniDamage;
	return out;
}

ElasticFractureInhomogFactors::ElasticFractureInhomogFactors()
{
	CFactor = 1.0;
	rhoFactor = 1.0;
	dampingFactor = 1.0;
	sigmaFactor = 1.0;
	deltaFactor = 1.0;
	iniDamage = 0.0;
}

ElasticFractureInhomogField::ElasticFractureInhomogField()
{
	serialNumberPtr = NULL;
	isPeriodicPtr = NULL;
	xMPtr = NULL;
	xmPtr = NULL;
	baseName_WOExt_oihf_elastic1 = "none";
	baseName_WOExt_oihf_elastic2 = "none";
	baseName_WOExt_oihf_elastic_damping = "none";

	eit = eih_undecided;

	baseName_WOExt_oihf_fracture1 = "none";
	baseName_WOExt_oihf_fracture2 = "none";
	fit = fih_undecided;
	baseName_WOExt_oihf_fracture_iniDamage = "none";

	//////////////////////////////////////
	// computed
	serialNumber_str = "";
	oihf_elastic_C = NULL;
	oihf_elastic_rho = NULL;
	oihf_elastic_damping = NULL;
	oihf_fracture_sigma = NULL;
	oihf_fracture_delta = NULL;
	oihf_fracture_iniDamage = NULL;

	resolutionFactor = 1;
	sso_elasticC = sso_mean_harmonic;
	sso_elasticRho = sso_mean_arithmetic;
	sso_elasticDamping = sso_mean_arithmetic;
	sso_sigma = sso_min;
	sso_delta = sso_mean_arithmetic;
	sso_iniDamage = sso_max;
}

ElasticFractureInhomogField::~ElasticFractureInhomogField()
{
	if (serialNumberPtr != NULL)
		delete serialNumberPtr;
	if (isPeriodicPtr != NULL)
		delete isPeriodicPtr;
	if (xMPtr != NULL)
		delete xMPtr;
	if (xmPtr != NULL)
		delete xmPtr;

	if (oihf_elastic_C != NULL)
		delete oihf_elastic_C;
	if (oihf_elastic_rho != NULL)
		delete oihf_elastic_rho;
	if (oihf_elastic_damping != NULL)
		delete oihf_elastic_damping;
	if (oihf_fracture_sigma != NULL)
		delete oihf_fracture_sigma;
	if (oihf_fracture_delta != NULL)
		delete oihf_fracture_delta;
	if (oihf_fracture_iniDamage != NULL)
		delete oihf_fracture_iniDamage;
}

ElasticFractureInhomogField::ElasticFractureInhomogField(const ElasticFractureInhomogField& other)
{
	serialNumberPtr = NULL;
	isPeriodicPtr = NULL;
	xMPtr = NULL;
	xmPtr = NULL;
	oihf_elastic_C = NULL;
	oihf_elastic_rho = NULL;
	oihf_elastic_damping = NULL;
	oihf_fracture_sigma = NULL;
	oihf_fracture_delta = NULL;
	oihf_fracture_iniDamage = NULL;
	(*this) = other;
}

ElasticFractureInhomogField& ElasticFractureInhomogField::operator=(const ElasticFractureInhomogField& other)
{
	if (other.serialNumberPtr != NULL)
	{
		if (serialNumberPtr == NULL)
			serialNumberPtr = new int();
		*serialNumberPtr = *other.serialNumberPtr;
	}
	if (other.isPeriodicPtr != NULL)
	{
		if (isPeriodicPtr == NULL)
			isPeriodicPtr = new bool();
		*isPeriodicPtr = *other.isPeriodicPtr;
	}
	if (other.xMPtr != NULL)
	{
		if (xMPtr == NULL)
			xMPtr = new double();
		*xMPtr = *other.xMPtr;
	}
	if (other.xmPtr != NULL)
	{
		if (xmPtr == NULL)
			xmPtr = new double();
		*xmPtr = *other.xmPtr;
	}
	baseName_WOExt_oihf_elastic1 = other.baseName_WOExt_oihf_elastic1;
	baseName_WOExt_oihf_elastic2 = other.baseName_WOExt_oihf_elastic2;
	baseName_WOExt_oihf_elastic_damping = other.baseName_WOExt_oihf_elastic_damping;
	eit = other.eit;
	baseName_WOExt_oihf_fracture1 = other.baseName_WOExt_oihf_fracture1;
	baseName_WOExt_oihf_fracture2 = other.baseName_WOExt_oihf_fracture2;
	fit = other.fit;
	baseName_WOExt_oihf_fracture_iniDamage = other.baseName_WOExt_oihf_fracture_iniDamage;
	serialNumber_str = other.serialNumber_str;

	if (other.oihf_elastic_C != NULL)
	{
		if (oihf_elastic_C == NULL)
			oihf_elastic_C = new OneIHField();
		*oihf_elastic_C = *other.oihf_elastic_C;
	}
	if (other.oihf_elastic_rho != NULL)
	{
		if (oihf_elastic_rho == NULL)
			oihf_elastic_rho = new OneIHField();
		*oihf_elastic_rho = *other.oihf_elastic_rho;
	}
	if (other.oihf_elastic_damping != NULL)
	{
		if (oihf_elastic_damping == NULL)
			oihf_elastic_damping = new OneIHField();
		*oihf_elastic_damping = *other.oihf_elastic_damping;
	}
	if (other.oihf_elastic_damping != NULL)
	{
		if (oihf_fracture_sigma == NULL)
			oihf_fracture_sigma = new OneIHField();
		*oihf_fracture_sigma = *other.oihf_fracture_sigma;
	}
	if (other.oihf_fracture_delta != NULL)
	{
		if (oihf_fracture_delta == NULL)
			oihf_fracture_delta = new OneIHField();
		*oihf_fracture_delta = *other.oihf_fracture_delta;
	}
	if (other.oihf_fracture_iniDamage != NULL)
	{
		if (oihf_fracture_iniDamage == NULL)
			oihf_fracture_iniDamage = new OneIHField();
		*oihf_fracture_iniDamage = *other.oihf_fracture_iniDamage;
	}
	xs = other.xs;
	numSegments = other.numSegments;
	numVertices = other.numVertices;

	baseName_WOExt_oihf_elastic1_inst = other.baseName_WOExt_oihf_elastic1_inst;
	baseName_WOExt_oihf_elastic2_inst = other.baseName_WOExt_oihf_elastic2_inst;
	baseName_WOExt_oihf_elastic_damping_inst = other.baseName_WOExt_oihf_elastic_damping_inst;
	baseName_WOExt_oihf_fracture1_inst = other.baseName_WOExt_oihf_fracture1_inst;
	baseName_WOExt_oihf_fracture2_inst = other.baseName_WOExt_oihf_fracture2_inst;
	baseName_WOExt_oihf_fracture_iniDamage_inst = other.baseName_WOExt_oihf_fracture_iniDamage_inst;

	resolutionFactor = other.resolutionFactor;
	sso_elasticC = other.sso_elasticC;
	sso_elasticRho = other.sso_elasticRho;
	sso_elasticDamping = other.sso_elasticDamping;
	sso_sigma = other.sso_sigma;
	sso_delta = other.sso_delta;
	sso_iniDamage = other.sso_iniDamage;

	return *this;
}
void ElasticFractureInhomogField::Read_ElasticFractureInhomogField(unsigned int subdomainNumber, string configFileName, int* serialNumberPtrIn, bool* isPeriodicPtrIn, double* xMPtrIn, double* xmPtrIn)
{
	fstream in(configFileName.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "configFileName\t" << configFileName << '\n';
		THROW("Cannot open config file\n");
	}
	Read_ElasticFractureInhomogField(in, serialNumberPtrIn, isPeriodicPtrIn, xMPtrIn, xmPtrIn, subdomainNumber);
}

void ElasticFractureInhomogField::Read_ElasticFractureInhomogField(istream & in, int* serialNumberPtrIn, bool *isPeriodicPtrIn, double *xMPtrIn, double *xmPtrIn, unsigned int subdomainNumber)
{
	string key, str_val; map<string, string>* mpPtr;
	setStatOp_type sso;
	double value = -1;
	subdomainNo = subdomainNumber;

	if (serialNumberPtrIn != NULL)
	{
		serialNumberPtr = new int();
		*serialNumberPtr = *serialNumberPtrIn;
	}
	if (isPeriodicPtrIn != NULL)
	{
		isPeriodicPtr = new bool();
		*isPeriodicPtr = *isPeriodicPtrIn;
	}
	if (xMPtrIn != NULL)
	{
		xMPtr = new double();
		*xMPtr = *xMPtrIn;
	}
	if (xmPtrIn != NULL)
	{
		xmPtr = new double();
		*xmPtr = *xmPtrIn;
	}
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
		else
		{
			string toCheck = "infile_inhomogeneous_elastic_frac";
			string subdomainNumber_str;
			toString(subdomainNumber, subdomainNumber_str);
			toCheck += subdomainNumber_str;
			while ((buf != toCheck) && (!in.eof()))
				READ_NSTRING(in, buf, buf);
			if (in.eof())
				THROW("Reached end of file looking for infile_icbc\n");
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
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "serialNumber")
		{
			int serialNumber;
			READ_NINTEGER(in, buf, serialNumber);
			if (serialNumberPtr == NULL)
			{
				serialNumberPtr = new int();
				*serialNumberPtr = serialNumber;
			}
		}
		else if (buf == "isPeriodic")
		{
			bool tmpb;
			READ_NBOOL(in, buf, tmpb);
			if (isPeriodicPtr == NULL)
			{
				isPeriodicPtr = new bool();
				*isPeriodicPtr = tmpb;
			}
		}
		else if (buf == "xM")
		{
			double tmp;
			READ_NDOUBLE(in, buf, tmp);
			if (xMPtr == NULL)
			{
				xMPtr = new double();
				*xMPtr = tmp;
			}
		}
		else if (buf == "xm")
		{
			double tmp;
			READ_NDOUBLE(in, buf, tmp);
			if (xmPtr == NULL)
			{
				xmPtr = new double();
				*xmPtr = tmp;
			}
		}
		else if (buf == "resolutionFactor")
		{
			READ_NINTEGER(in, buf, resolutionFactor);
			if (sfcm.direct_resolutionFactor != 0)
				resolutionFactor = sfcm.direct_resolutionFactor;
			key = "resFact";
			bool found = Find_Version_Value(key, value, mpPtr);
			if (!found)
			{
				key = "resolutionFactor";
				found = Find_Version_Value(key, value, mpPtr);
			}
			if (found)
				resolutionFactor = (int)round(value);
		}
		else if (buf == "sso_elasticC")
		{
			in >> sso_elasticC;
			key = "ssoEE";
			bool found = Find_Version_String(key, str_val, mpPtr);
			if (found && (name2Type(str_val, sso) == true))
				sso_elasticC = sso;
		}
		else if (buf == "sso_elasticRho")
		{
			in >> sso_elasticRho;
			key = "ssoER";
			bool found = Find_Version_String(key, str_val, mpPtr);
			if (found && (name2Type(str_val, sso) == true))
				sso_elasticRho = sso;
		}
		else if (buf == "sso_elasticDamping")
		{
			in >> sso_elasticDamping;
			key = "ssoED";	str_val = "none";
			bool found = Find_Version_String(key, str_val, mpPtr);
			if (found && (name2Type(str_val, sso) == true))
				sso_elasticDamping = sso;
		}
		else if (buf == "sso_sigma")
		{
			in >> sso_sigma;
			key = "ssoFS";	str_val = "none";
			bool found = Find_Version_String(key, str_val, mpPtr);
			if (found && (name2Type(str_val, sso) == true))
				sso_sigma = sso;
		}
		else if (buf == "sso_delta")
		{
			in >> sso_delta;
			key = "ssoFD";	str_val = "none";
			bool found = Find_Version_String(key, str_val, mpPtr);
			if (found && (name2Type(str_val, sso) == true))
				sso_delta = sso;
		}
		else if (buf == "sso_iniDamage")
		{
			in >> sso_iniDamage;
			key = "ssoFI";	str_val = "none";
			bool found = Find_Version_String(key, str_val, mpPtr);
			if (found && (name2Type(str_val, sso) == true))
				sso_iniDamage = sso;
		}
		else if (buf == "baseName_WOExt_oihf_elastic1")
		{
			READ_NSTRING(in, buf, baseName_WOExt_oihf_elastic1);
			if (baseName_WOExt_oihf_elastic1 == "inline")
			{
				READ_NSTRING(in, buf, baseName_WOExt_oihf_elastic1);
				baseName_WOExt_oihf_elastic1_inst = g_prefileName + "/_baseName_WOExt_oihf_elastic1.inst";
				WriteContent_BetweenCurlyBrack2OtherFile(in, baseName_WOExt_oihf_elastic1_inst);
			}
			else
				baseName_WOExt_oihf_elastic1_inst = baseName_WOExt_oihf_elastic1 + ".inst";
		}
		else if (buf == "baseName_WOExt_oihf_elastic2")
		{
			READ_NSTRING(in, buf, baseName_WOExt_oihf_elastic2);
			if (baseName_WOExt_oihf_elastic2 == "inline")
			{
				READ_NSTRING(in, buf, baseName_WOExt_oihf_elastic2);
				baseName_WOExt_oihf_elastic2_inst = g_prefileName + "/_baseName_WOExt_oihf_elastic2.inst";
				WriteContent_BetweenCurlyBrack2OtherFile(in, baseName_WOExt_oihf_elastic2_inst);
			}
			else
				baseName_WOExt_oihf_elastic2_inst = baseName_WOExt_oihf_elastic2 + ".inst";
		}
		else if (buf == "baseName_WOExt_oihf_elastic_damping")
		{
			READ_NSTRING(in, buf, baseName_WOExt_oihf_elastic_damping);
			if (baseName_WOExt_oihf_elastic_damping == "inline")
			{
				READ_NSTRING(in, buf, baseName_WOExt_oihf_elastic_damping);
				baseName_WOExt_oihf_elastic_damping_inst = g_prefileName + "/_baseName_WOExt_oihf_elastic_damping.inst";
				WriteContent_BetweenCurlyBrack2OtherFile(in, baseName_WOExt_oihf_elastic_damping_inst);
			}
			else
				baseName_WOExt_oihf_elastic_damping_inst = baseName_WOExt_oihf_elastic_damping + ".inst";
		}
		else if (buf == "eit")
		{
			in >> eit;
		}
		else if (buf == "baseName_WOExt_oihf_fracture1")
		{
			READ_NSTRING(in, buf, baseName_WOExt_oihf_fracture1);
			if (baseName_WOExt_oihf_fracture1 == "inline")
			{
				READ_NSTRING(in, buf, baseName_WOExt_oihf_fracture1);
				baseName_WOExt_oihf_fracture1_inst = g_prefileName + "/_baseName_WOExt_oihf_fracture1.inst";
				WriteContent_BetweenCurlyBrack2OtherFile(in, baseName_WOExt_oihf_fracture1_inst);
			}
			else
				baseName_WOExt_oihf_fracture1_inst = baseName_WOExt_oihf_fracture1 + ".inst";
		}
		else if (buf == "baseName_WOExt_oihf_fracture2")
		{
			READ_NSTRING(in, buf, baseName_WOExt_oihf_fracture2);
			if (baseName_WOExt_oihf_fracture2 == "inline")
			{
				READ_NSTRING(in, buf, baseName_WOExt_oihf_fracture2);
				baseName_WOExt_oihf_fracture2_inst = g_prefileName + "/_baseName_WOExt_oihf_fracture2.inst";
				WriteContent_BetweenCurlyBrack2OtherFile(in, baseName_WOExt_oihf_fracture2_inst);
			}
			else
				baseName_WOExt_oihf_fracture2_inst = baseName_WOExt_oihf_fracture2 + ".inst";
		}
		else if (buf == "baseName_WOExt_oihf_fracture_iniDamage")
		{
			READ_NSTRING(in, buf, baseName_WOExt_oihf_fracture_iniDamage);
			if (baseName_WOExt_oihf_fracture_iniDamage == "inline")
			{
				READ_NSTRING(in, buf, baseName_WOExt_oihf_fracture_iniDamage);
				baseName_WOExt_oihf_fracture_iniDamage_inst = g_prefileName + "/_baseName_WOExt_oihf_fracture_iniDamage.inst";
				WriteContent_BetweenCurlyBrack2OtherFile(in, baseName_WOExt_oihf_fracture_iniDamage_inst);
			}
			else
				baseName_WOExt_oihf_fracture_iniDamage_inst = baseName_WOExt_oihf_fracture_iniDamage + ".inst";
		}
		else if (buf == "fit")
		{
			in >> fit;
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Initialize_ElasticFractureInhomogField();
}

void ElasticFractureInhomogField::getFactorsIniDamage_ByIndex(unsigned int index, ElasticFractureInhomogFactors & efif) const
{
	if (oihf_elastic_C != NULL)
	{
		if (eit != eih_Crho_const)
			efif.CFactor = oihf_elastic_C->getValueByIndex(index);
		if (eit == eih_c_const)
			efif.rhoFactor = efif.CFactor;
		else if (eit == eih_Z_const)
			efif.rhoFactor = 1.0 / efif.CFactor;
		else if (eit == eih_nonconst)
			efif.rhoFactor = oihf_elastic_rho->getValueByIndex(index);
	}
	else if (eit == eih_C_const)
		efif.rhoFactor = oihf_elastic_rho->getValueByIndex(index);
	if (oihf_elastic_damping != NULL)
		efif.dampingFactor = oihf_elastic_damping->getValueByIndex(index);
	else
		efif.dampingFactor = 1.0;
	
	// fracture ones
	if (oihf_fracture_iniDamage != NULL)
		efif.iniDamage = MAX(oihf_fracture_iniDamage->getValueByIndex(index), 0.0);

	if (oihf_fracture_sigma != NULL)
	{
		efif.sigmaFactor = oihf_fracture_sigma->getValueByIndex(index);
		if (fit == fih_G_const)
			efif.deltaFactor = 1.0 / efif.sigmaFactor;
		else if (fit == fih_nonconst)
			efif.deltaFactor = oihf_fracture_delta->getValueByIndex(index);
	}
	else if (fit == fih_sigma_const)
		efif.deltaFactor = oihf_fracture_delta->getValueByIndex(index);
}

void ElasticFractureInhomogField::getFactorsIniDamage_By_x(double x, ElasticFractureInhomogFactors & efif) const
{
	if (oihf_elastic_C != NULL)
	{
		if (eit != eih_Crho_const)
			efif.CFactor = oihf_elastic_C->getValueByCoord(x);
		if (eit == eih_c_const)
			efif.rhoFactor = efif.CFactor;
		else if (eit == eih_Z_const)
			efif.rhoFactor = 1.0 / efif.CFactor;
		else if (eit == eih_nonconst)
			efif.rhoFactor = oihf_elastic_rho->getValueByCoord(x);
	}
	else if (eit == eih_C_const)
		efif.rhoFactor = oihf_elastic_rho->getValueByCoord(x);
	if (oihf_elastic_damping != NULL)
		efif.dampingFactor = oihf_elastic_damping->getValueByCoord(x);
	else
		efif.dampingFactor = 1.0;

	// fracture ones
	if (oihf_fracture_iniDamage != NULL)
		efif.iniDamage = MAX(oihf_fracture_iniDamage->getValueByCoord(x), 0.0);

	if (oihf_fracture_sigma != NULL)
	{
		efif.sigmaFactor = oihf_fracture_sigma->getValueByCoord(x);
		if (fit == fih_G_const)
			efif.deltaFactor = 1.0 / efif.sigmaFactor;
		else if (fit == fih_nonconst)
			efif.deltaFactor = oihf_fracture_delta->getValueByCoord(x);
	}
	else if (fit == fih_sigma_const)
		efif.deltaFactor = oihf_fracture_delta->getValueByCoord(x);
}

unsigned int ElasticFractureInhomogField::getFactorsIniDamage_4AllIndices(vector<ElasticFractureInhomogFactors>& efifs) const
{
	efifs.clear();
	unsigned int sz = numVertices;
	efifs.resize(numVertices);
	for (unsigned int i = 0; i < sz; ++i)
		getFactorsIniDamage_ByIndex(i, efifs[i]);
	return sz;
}

unsigned int ElasticFractureInhomogField::getFactorsIniDamage_By_Equai_distance_xs(vector<ElasticFractureInhomogFactors>& efifs, vector<double>& xss, int numSegmentsIn, double* xMPtrIn, double* xmPtrIn) const
{
	efifs.clear();
	if (numSegmentsIn < 0)
		numSegmentsIn *= -numSegments;
	double xm, xM;
	if (xMPtrIn != NULL)
		xM = *xMPtrIn;
	else
		xM = *xMPtr;
	if (xmPtrIn != NULL)
		xm = *xmPtrIn;
	else
		xm = *xmPtr;
	unsigned int numVert = numSegmentsIn + 1;
	double del = (xM - xm) / numSegmentsIn;

	efifs.resize(numVert);
	xss.resize(numVert);
	for (unsigned int i = 0; i < numVert; ++i)
	{
		double x = xm + del * i;
		xss[i] = x;
		getFactorsIniDamage_By_x(x, efifs[i]);
	}
	return numVert;
}

void ElasticFractureInhomogField::Initialize_ElasticFractureInhomogField()
{
#if 0
	string key = "dd2";	map<string, string>* mpPtr;
	double value;
	if ((Find_Version_Value(key, value, mpPtr)) && (value < 1e-7))
	{
		baseName_WOExt_oihf_elastic1 = "none";
		baseName_WOExt_oihf_elastic2 = "none";
		baseName_WOExt_oihf_elastic_damping = "none";

		eit = eih_Crho_const;

		baseName_WOExt_oihf_fracture1 = "none";
		baseName_WOExt_oihf_fracture2 = "none";
		fit = fih_sigmadelta_const;
		baseName_WOExt_oihf_fracture_iniDamage = "none";
	}
#endif
	if (serialNumberPtr != NULL)
	{
		toString(*serialNumberPtr, serialNumber_str);
		serialNumber_str = "_" + serialNumber_str;
	}
	else
		serialNumber_str = "";

	bool b_baseName_WOExt_oihf_elastic1, b_baseName_WOExt_oihf_elastic2;
	b_baseName_WOExt_oihf_elastic1 = (baseName_WOExt_oihf_elastic1 != "none");
	b_baseName_WOExt_oihf_elastic2 = (baseName_WOExt_oihf_elastic2 != "none");
	if (b_baseName_WOExt_oihf_elastic1)
	{
		if (b_baseName_WOExt_oihf_elastic2)
		{
			// both files are given
			if (eit == eih_undecided)
				eit = eih_nonconst;
		}
		else // first files is given, there is no second file
		{
			if (eit == eih_undecided)
				eit = eih_rho_const;
		}
	}
	else // none of the files are given
	{
		eit = eih_Crho_const;
	}
	if (oihf_elastic_C != NULL)
	{
		delete oihf_elastic_C;
		oihf_elastic_C = NULL;
	}
	if (oihf_elastic_rho != NULL)
	{
		delete oihf_elastic_rho;
		oihf_elastic_rho = NULL;
	}
	if (oihf_elastic_damping != NULL)
	{
		delete oihf_elastic_damping;
		oihf_elastic_damping = NULL;
	}
	if ((eit == eih_rho_const) || (eit == eih_c_const) || (eit == eih_Z_const) || (eit == eih_nonconst))
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_elastic1, 1, 0);
		string fileNameDat = baseName_WOExt_oihf_elastic1 + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_elastic1_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_elastic_C = new OneIHField();
		oihf_elastic_C->Read_Initialize_OneIHField(inDat, inInstruction, 0, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_elasticC);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	else if (eit == eih_C_const)
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_elastic1, 1, 0);
		string fileNameDat = baseName_WOExt_oihf_elastic1 + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_elastic1_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_elastic_rho = new OneIHField();
		oihf_elastic_rho->Read_Initialize_OneIHField(inDat, inInstruction, 0, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_elasticRho);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	if (eit == eih_nonconst)
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_elastic2, 2, 0);
		string fileNameDat = baseName_WOExt_oihf_elastic2 + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_elastic2_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_elastic_rho = new OneIHField();
		oihf_elastic_rho->Read_Initialize_OneIHField(inDat, inInstruction, 0, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_elasticRho);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	///////// iniDamage
	if (baseName_WOExt_oihf_elastic_damping != "none")
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_elastic_damping, 3, 0);
		string fileNameDat = baseName_WOExt_oihf_elastic_damping + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_elastic_damping_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_elastic_damping = new OneIHField();
		oihf_elastic_damping->Read_Initialize_OneIHField(inDat, inInstruction, 0, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_elasticDamping);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	/////////////////// elastic files read above

	/////////////////// fracture files read below
	bool b_baseName_WOExt_oihf_fracture1, b_baseName_WOExt_oihf_fracture2;
	b_baseName_WOExt_oihf_fracture1 = (baseName_WOExt_oihf_fracture1 != "none");
	b_baseName_WOExt_oihf_fracture2 = (baseName_WOExt_oihf_fracture2 != "none");
	if (b_baseName_WOExt_oihf_fracture1)
	{
		if (b_baseName_WOExt_oihf_fracture2)
		{
			// both files are given
			if (fit == fih_undecided)
				fit = fih_nonconst;
		}
		else // first files is given, there is no second file
		{
			if (fit == fih_undecided)
				fit = fih_delta_const;
		}
	}
	else // none of the files are given
	{
		fit = fih_sigmadelta_const;
	}
	if (oihf_fracture_sigma != NULL)
	{
		delete oihf_fracture_sigma;
		oihf_fracture_sigma = NULL;
	}
	if (oihf_fracture_delta != NULL)
	{
		delete oihf_fracture_delta;
		oihf_fracture_delta = NULL;
	}
	if ((fit == fih_delta_const) || (fit == fih_G_const) || (fit == fih_nonconst))
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_fracture1, 0, 1);
		string fileNameDat = baseName_WOExt_oihf_fracture1 + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_fracture1_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_fracture_sigma = new OneIHField();
		oihf_fracture_sigma->Read_Initialize_OneIHField(inDat, inInstruction, 1, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_sigma);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	else if (fit == fih_sigma_const)
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_fracture1, 0, 1);
		string fileNameDat = baseName_WOExt_oihf_fracture1 + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_fracture1_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_fracture_delta = new OneIHField();
		oihf_fracture_delta->Read_Initialize_OneIHField(inDat, inInstruction, 1, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_delta);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	if (fit == fih_nonconst)
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_fracture2, 0, 2);
		string fileNameDat = baseName_WOExt_oihf_fracture2 + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_fracture2_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_fracture_delta = new OneIHField();
		oihf_fracture_delta->Read_Initialize_OneIHField(inDat, inInstruction, 1, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_delta);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	///////// iniDamage
	if (baseName_WOExt_oihf_fracture_iniDamage != "none")
	{
		ModifyElasticFracture_RandomFieldInput(baseName_WOExt_oihf_fracture_iniDamage, 0, 3);
		string fileNameDat = baseName_WOExt_oihf_fracture_iniDamage + serialNumber_str + ".txt";
		string fileNameInstruction = baseName_WOExt_oihf_fracture_iniDamage_inst;
		fstream inDat(fileNameDat.c_str(), ios::in);
		if (!inDat.is_open())
		{
			cout << fileNameDat << '\n';
			THROW("Cannot open file\n");
		}
		fstream *inInstruction = new fstream();
		inInstruction->open(fileNameInstruction.c_str(), ios::in);
		if (!inInstruction->is_open())
		{
			delete inInstruction;
			inInstruction = NULL;
		}
		oihf_fracture_iniDamage = new OneIHField();
		oihf_fracture_iniDamage->Read_Initialize_OneIHField(inDat, inInstruction, 1, isPeriodicPtr, xMPtr, xmPtr, resolutionFactor, sso_iniDamage);
		if (inInstruction != NULL)
			delete inInstruction;
	}
	//////////////////////////////////////////////
	// finding the number of vertices, segments, xs ...
	OneIHField *oihf = NULL;
	if (oihf_elastic_C != NULL)
		oihf = oihf_elastic_C;
	else if (oihf_fracture_sigma != NULL)
		oihf = oihf_fracture_sigma;
	else if (oihf_fracture_iniDamage != NULL)
		oihf = oihf_fracture_iniDamage;
	else if (oihf_fracture_delta != NULL)
		oihf = oihf_fracture_delta;
	else if (oihf_elastic_rho != NULL)
		oihf = oihf_elastic_rho;

	numSegments = 0;
	numVertices = 0;
	if (oihf != NULL)
	{
		numSegments = oihf->getNumSegments();
		numVertices = oihf->getNumVertices();
		oihf->get_xs(xs);
		if (xmPtr == NULL)
		{
			xmPtr = new double();
			*xmPtr = oihf->get_xm();
		}
		if (xMPtr == NULL)
		{
			xMPtr = new double();
			*xMPtr = oihf->get_xM();
		}
	}
}

void ElasticFractureInhomogField::ModifyElasticFracture_RandomFieldInput(string& baseName_WOExt, int elasticNo, int factureNo)
{
	bool versionChange = ((sfcm.success) && (sfcm_gen.subdomainNo == subdomainNo));
	if (!versionChange)
		return;
	string llc_str = "llcnone";
	string key;	map<string, string>* mpPtr;
	key = "llc";
	if (!Find_Version_String(key, llc_str, mpPtr))
		return;
	g_logout << "\tllc_str\t" << llc_str;

//	string np_str = "1025";
	string np_str = "16385";
/*
	if (mpPtr != NULL)
	{
		map<string, string>::iterator it = mpPtr->find("np");
		if (it != mpPtr->end())
			np_str = it->second;
	}
*/
	if (sfcm.success && (sfcm.number_of_elements > 0))
	{
		int np = sfcm.number_of_elements;
//		if (factureNo > 0)
			++np;
		toString(np, np_str);
	}
	double llc = -1;
	double d_np;
	double limt = -4.9999;
	if (fromString(np_str, d_np) == true)
	{
		limt = -(log(d_np) / log(10.0));
	}
	fromString(llc_str, llc);
	if (llc >= limt)//(llc <= -5.9999)
		baseName_WOExt = "InhomogeneousFiles/cl" + llc_str + "_np" + np_str + "/initial_values";
	else
		baseName_WOExt = "InhomogeneousFiles/clz_np" + np_str + "/initial_values";
}

void WriteContent_BetweenCurlyBrack2OtherFile(istream& in, const string& otherFileName)
{
	string buf;
	fstream out(otherFileName.c_str(), ios::out);
	READ_NSTRING(in, buf, buf);
	while (buf != "endinline")
	{
		out << buf << '\n';
		READ_NSTRING(in, buf, buf);
	}
}

void TestInhomogeneousElasticFractorField(string configName, int* serialNumberPtrIn, bool * isPeriodicPtrIn, double * xMPtrIn, double * xmPtrIn)
{
	/// Giang: how to use this
	// your input files are:
	// TestInhomogeneousElasticFractorField_config.txt -> this gives the name of initial Damage file
	// Data
	// TestInhomogeneousElasticFractorField_fracture_iniDamage_2.txt: the 101 values here are NOT initial damage. They are supposedly from a standard normal distribution
	// Instruction
	// TestInhomogeneousElasticFractorField_fracture_iniDamage.inst

	// NOTE: you see that there is a _2 for the data file. You can have from _0 to say _1000 (1001 runs) data files for initial damage
	// serialNumber 2 is the serial number of the run. 
	// serial number is not used in instruction file	

	ElasticFractureInhomogField efif;
	efif.Read_ElasticFractureInhomogField(0, configName, serialNumberPtrIn, isPeriodicPtrIn, xMPtrIn, xmPtrIn);

	// how to use it in your code
	ElasticFractureInhomogFactors efiFactor;
	double x = 0.65;
	efif.getFactorsIniDamage_By_x(x, efiFactor);
	double iniDamage_at_x = efiFactor.iniDamage;

	// End for Giang

	vector< ElasticFractureInhomogFactors> efifs;
	fstream out("TestFiles/TestInhomogeneousElasticFractorField_index.txt", ios::out);
	unsigned int sz = efif.getFactorsIniDamage_4AllIndices(efifs);
	for (unsigned int i = 0; i < sz; ++i)
		out << i << '\t' << efifs[i] << '\n';
	out.close();

	vector<double> xss;
	out.open("TestFiles/TestInhomogeneousElasticFractorField_x.txt", ios::out);
	int numSegmentsIn = -10;
	sz = efif.getFactorsIniDamage_By_Equai_distance_xs(efifs, xss, numSegmentsIn);
	for (unsigned int i = 0; i < sz; ++i)
		out << xss[i] << '\t' << efifs[i] << '\n';
}
