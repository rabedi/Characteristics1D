#ifndef INHOMOGENEOUS_Elastic_Fracture__H
#define INHOMOGENEOUS_Elastic_Fracture__H

#include "InhomogeousField.h"

// C stands for stiffness
// rho mass density
// c: wave speed
// Z: impedance
//		CPtr: Inhomogeneous generator for C				rhoPtr: Inhomogeneous generator for rho
//															//	CPtr			//	rhoPtr
// eih_Crho_const:	C & rho constant							NULL				NULL
// eih_rho_const:	rho constant								CPtr				NULL
// eih_C_const:		C constant									NULL				rhoPtr
// eih_c_const:		c constant (rhoFactor = CFractor)			CPtr				NULL
// eih_Z_const:		Z constant (rhoFactor = 1/CFractor)			CPtr				NULL
// eih_nonconst:	nothing is constant							CPtr				rhoPtr
typedef enum { eih_undecided, eih_Crho_const, eih_rho_const, eih_C_const, eih_c_const, eih_Z_const, eih_nonconst, elasticInhomogT_SIZE } elasticInhomogT;

string getName(elasticInhomogT dat);
void name2Type(string& name, elasticInhomogT& typeVal);
ostream& operator<<(ostream& out, elasticInhomogT dat);
istream& operator>>(istream& in, elasticInhomogT& dat);


// sigma stands for strength scale
// delta displacement scale
// G: fracture energy scale
//		sigmaPtr: Inhomogeneous generator for sigma				deltaPtr: Inhomogeneous generator for delta
//																//	sigmaPtr			//	deltaPtr
// fih_sigmadelta_const:	sigma & delta constant					NULL					NULL
// fih_delta_const:			delta constant							sigmaPtr				NULL
// fih_sigma_const:			sigma constant							NULL					deltaPtr
// fih_G_const:		G constant (deltaFactor = 1/sigmaFactor)		sigmaPtr				NULL
// fih_nonconst:			nothing is constant						sigmaPtr				deltaPtr
typedef enum { fih_undecided, fih_sigmadelta_const, fih_delta_const, fih_sigma_const, fih_G_const, fih_nonconst, fractureInhomogT_SIZE } fractureInhomogT;

string getName(fractureInhomogT dat);
void name2Type(string& name, fractureInhomogT& typeVal);
ostream& operator<<(ostream& out, fractureInhomogT dat);
istream& operator>>(istream& in, fractureInhomogT& dat);

class ElasticFractureInhomogFactors
{
	friend ostream& operator<<(ostream& out, const ElasticFractureInhomogFactors& dat);
public:
	ElasticFractureInhomogFactors();
	// factors for elastic and fracture properties and initial damage
	// elastic
	double CFactor, rhoFactor, dampingFactor;
	// fracture
	double sigmaFactor, deltaFactor;
	double iniDamage;
};

class ElasticFractureInhomogField
{
public:
	ElasticFractureInhomogField();
	~ElasticFractureInhomogField();
	void Read_ElasticFractureInhomogField(unsigned int subdomainNumber, string configFileName, int* serialNumberPtrIn = NULL, bool* isPeriodicPtrIn = NULL, double* xMPtrIn = NULL, double* xmPtrIn = NULL);
public:	
	// access functions
	////////////////////////////////////////////////
	// index is the index of bulk (for elastic) or interface (for fracture)
	void getFactorsIniDamage_ByIndex(unsigned int index, ElasticFractureInhomogFactors& efif) const;
	// x is spatial position
	void getFactorsIniDamage_By_x(double x, ElasticFractureInhomogFactors & efif) const;
	/// useful functions to access all at once
	// return value is size
	unsigned int getFactorsIniDamage_4AllIndices(vector< ElasticFractureInhomogFactors>& efifs) const;
	// numSegmentsIn < 0 -> it's a factor of numSegments (for example -3 -> 3 * numSegments is used)
	// if xMPtrIn or xmPtrIn are NULL, inside xM and xm pointers are used to decide x bounds
	// return value is the size of efifs
	unsigned int getFactorsIniDamage_By_Equai_distance_xs(vector< ElasticFractureInhomogFactors>& efifs, vector<double>& xss, int numSegmentsIn = -1, double* xMPtrIn = NULL, double* xmPtrIn = NULL) const;

	// serial number for random fields
	int* serialNumberPtr;

	// isPeriodic and domain bounds (if !NULL they overwrite individual field values)
	bool* isPeriodicPtr;
	double* xMPtr; double* xmPtr;

	// resolutionFactor 
	//					== 0 or +/-1 nothing happens
	//					>  1  -> number of segments is DECREASED by this factor (e.g. if resolutionFactor ==  10 and numSegments = 1000 -> numSegments becoms 100)
	//					<  -1 -> number of segments is INCREASED by this factor (e.g. if resolutionFactor == -10 and numSegments = 1000 -> numSegments becoms 10000)
	int resolutionFactor;

	// elastic fields
	////////////////////////////////////////
	// 1 and 2 correspond to random fields for C factor and rhoFactor except only for eih_C_const where rho is given from the first field
	// these are base names for elastic 1 and 2 data files (will get serial number if serialNumber >= 0) and instruction file
	string baseName_WOExt_oihf_elastic1, baseName_WOExt_oihf_elastic2;
	string baseName_WOExt_oihf_elastic_damping;
	// for certain cases, the type of elastic field must be provided, otherwise the default will be chosen
	elasticInhomogT eit;
	// if the read elastic field is going to be coarsened, we need to know how to do the coarsening (averaging) the values
	setStatOp_type sso_elasticC, sso_elasticRho, sso_elasticDamping;

	////////////////////////////////////////
	// 1 and 2 correspond to random fields for sigma and delta fields except only for fih_sigma_const where delta is given from the first one
	string baseName_WOExt_oihf_fracture1, baseName_WOExt_oihf_fracture2;
	// for certain cases, the type of fracture field must be provided, otherwise the default will be chosen
	fractureInhomogT fit;
	// initial damage base file name
	string baseName_WOExt_oihf_fracture_iniDamage;
	setStatOp_type sso_sigma, sso_delta, sso_iniDamage;

	//////////////////////////////////////
	// computed
	string serialNumber_str;
	OneIHField *oihf_elastic_C, *oihf_elastic_rho, *oihf_elastic_damping;
	OneIHField *oihf_fracture_sigma, *oihf_fracture_delta, *oihf_fracture_iniDamage;

	vector<double> xs;
	int numSegments, numVertices;

	ElasticFractureInhomogField(const ElasticFractureInhomogField& other);
	ElasticFractureInhomogField& operator=(const ElasticFractureInhomogField& other);
private:
	string baseName_WOExt_oihf_elastic1_inst, baseName_WOExt_oihf_elastic2_inst;
	string baseName_WOExt_oihf_elastic_damping_inst;
	string baseName_WOExt_oihf_fracture1_inst, baseName_WOExt_oihf_fracture2_inst;
	string baseName_WOExt_oihf_fracture_iniDamage_inst;
	int subdomainNo;

	void Read_ElasticFractureInhomogField(istream& in, int* serialNumberPtrIn = NULL, bool* isPeriodicPtrIn = NULL, double* xMPtrIn = NULL, double* xmPtrIn = NULL, unsigned int subdomainNumber = 0);
	void Initialize_ElasticFractureInhomogField();


	// elasticNo and factureNo refer to files below. for example elasticNo = 0, fractureNo = 1 refers to fracture1 file
//	Elastic  numbers 1, 2, 3	baseName_WOExt_oihf_elastic1, baseName_WOExt_oihf_elastic2, baseName_WOExt_oihf_elastic_damping
//  Fracture numbers 1, 2, 3 baseName_WOExt_oihf_fracture1, baseName_WOExt_oihf_fracture2, baseName_WOExt_oihf_fracture_iniDamage;
	void ModifyElasticFracture_RandomFieldInput(string& baseName_WOExt, int elasticNo, int factureNo);
};

void WriteContent_BetweenCurlyBrack2OtherFile(istream& in, const string& otherFileName);
void TestInhomogeneousElasticFractorField(string configName = "TestFiles/TestInhomogeneousElasticFractorField_config.txt", int* serialNumberPtrIn = NULL, bool* isPeriodicPtrIn = NULL, double* xMPtrIn = NULL, double* xmPtrIn = NULL);

#endif