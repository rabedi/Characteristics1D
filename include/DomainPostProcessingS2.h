#ifndef DOMAIN_POST_PROCESSING_S2__H
#define DOMAIN_POST_PROCESSING_S2__H


#include "DomainPostProcessing.h"
#include "RandomVariable.h"
#include "globalTypesClasses.h"

//////////////////////////////////////////
// This file is Level 2
// 3 Levels of Post-Processing in C++
// Level 1.				Space + Spacetime integrals are computed for EVERY time step in "RUNNAME/sd_[subdomainNo]__Summary.txt
//						plus extraction of D, S, U for ALL INTERFACES and ALL BULKS -> see files sd_[subdomainNo]_BulkInterface_tPos_[timeIndex].txt		
//						Done during the solution (updated after each time step)
//						Files: DomainPostProcessing.*
// Level 2				Gets outputs from Level 1 and post-process new results including
//							Stages: lin -> response starts to deviate from linear, max (max stress, U, ... atained), an final (full damage) are characterized based on different classification criteria
//							Fragmentation data: Sizes of fragments and their statistics (min, max, number, ... of framents)
//						Process can be time multiple times with different parameters and flags after solution is completed
//						Files: DomainPostProcessingS2.*
// Level 3				Loads results from Level2 and generates CSV files for ML and other simple plotting/processing
//						Files: DomainPostProcessingS3.*


// loadingStagesT_none: given data is not time related (e.g. initial profile of strength, etc) 

// lst_actualTime: A real time is given either through indices (for one of the index types in timeIndexType) or actual time value
////// those related to fracture simulation
// lin -> deviate from linear response
// max -> failure is mostly active and/or max load/energy experienced
// final -> zero load/energy or final time of simulation
typedef enum { loadingStagesT_none = -2, lst_actualTime, lst_lin, lst_max, lst_final, loadingStagesT_SIZE} loadingStagesT;

string getName(loadingStagesT dat);
void name2Type(string& name, loadingStagesT& typeVal);
ostream& operator<<(ostream& out, loadingStagesT dat);
istream& operator>>(istream& in, loadingStagesT& dat);

// various classification criteria are used to determin lin, max, and final stages above.
// ct_time: solution time is used
// ct_epssig_bc:				eps & sig homogenized over bc are used
// ct_epssig_bulk:				eps & sig homogenized over bulk domain are used
// ct_interfaceDissEnee:		energy dissipation is used
// ct_Damage:					damage field is used
// ct_U:						Internal energy is used
typedef enum { lstClassificationT_none = -1, ct_time, ct_epssig_bulk, ct_epssig_bc, ct_interfaceDissEne, ct_Damage, ct_U, 
ct_frag_Damage, ct_frag_maxEffDelU, ct_frag_DelU, lstClassificationT_SIZE } lstClassificationT;
#define lstClassificationT_no_frag_SIZE	ct_frag_Damage

string getName(lstClassificationT dat);
void name2Type(string& name, lstClassificationT& typeVal);
ostream& operator<<(ostream& out, lstClassificationT dat);
istream& operator>>(istream& in, lstClassificationT& dat);

typedef enum { BrittlenessIndicatorRatioType_none = -1, birt_i2m, birt_i2f, birt_m2f, BrittlenessIndicatorRatioType_SIZE} BrittlenessIndicatorRatioType;
string getLatexName(BrittlenessIndicatorRatioType dat);
string getName(BrittlenessIndicatorRatioType dat);
void name2Type(string& name, BrittlenessIndicatorRatioType& typeVal);
ostream& operator<<(ostream& out, BrittlenessIndicatorRatioType dat);
istream& operator>>(istream& in, BrittlenessIndicatorRatioType& dat);

typedef enum { BrittlnessIndicatorFieldType_none = -1, birt_time, birt_eps, birt_psi, BrittlnessIndicatorFieldType_SIZE } BrittlnessIndicatorFieldType;
string getLatexName(BrittlnessIndicatorFieldType dat);
string getName(BrittlnessIndicatorFieldType dat);
void name2Type(string& name, BrittlnessIndicatorFieldType& typeVal);
ostream& operator<<(ostream& out, BrittlnessIndicatorFieldType dat);
istream& operator>>(istream& in, BrittlnessIndicatorFieldType& dat);

typedef enum { FragmentationCriterionT_none = -1, fct_Damage, fct_maxEffDelU, fct_DelU, FragmentationCriterionT_SIZE } FragmentationCriterionT;
string getLatexName(FragmentationCriterionT dat);
string getName(FragmentationCriterionT dat);
void name2Type(string& name, FragmentationCriterionT& typeVal);
ostream& operator<<(ostream& out, FragmentationCriterionT dat);
istream& operator>>(istream& in, FragmentationCriterionT& dat);

typedef enum { RunNormalizationQuantT_none = -1, RunNormalizationQuantT_default, 
	ene_KInitial, ene_UInitial, ene_PhiInitial, ene_KFinal, ene_UFinal, ene_PhiFinal, enePhiDissLossFinal, eneBCFinal, eneInputFinal,
	rnq_loadTimeScale, rnq_EScale, rnq_rhoScale, rnq_dampingScale, rnq_sigmaCScale, rnq_deltaCScale, rnq_energyCScale, rnq_cScale, rnq_ZScale, rnq_tauCScale, rnq_strainCScale, rnq_velCScale, rnq_psiCScale,
	RunNormalizationQuantT_SIZE} RunNormalizationQuantT;

string getLatexName(RunNormalizationQuantT dat);
string getName(RunNormalizationQuantT dat);
void name2Type(string& name, RunNormalizationQuantT& typeVal);
ostream& operator<<(ostream& out, RunNormalizationQuantT dat);
istream& operator>>(istream& in, RunNormalizationQuantT& dat);

typedef enum {sfo_scalar, sfo_field_x, sfo_field_t, sfo_field_other, ScalarFieldOutputT_SIZE} ScalarFieldOutputT;

string getLatexName(ScalarFieldOutputT dat);
string getName(ScalarFieldOutputT dat);
bool name2Type(string& name, ScalarFieldOutputT& typeVal);
ostream& operator<<(ostream& out, ScalarFieldOutputT dat);
istream& operator>>(istream& in, ScalarFieldOutputT& dat);

// There are multiple outputs of the code output at different resolutions (each certain number of time steps).
// The corresponding time steps are given in the main config files and set values are printed in _sd_[subdomainNo]_keyParameters.txt
// These time steps are listed below

//--------------------
// timeIndexType_none:						no time index
//--------------------
// ti_AllTimeSteps:							index for the solution time step -> corresponds to the rows of the file
//
//			example file					sd_[subdomainNo]__Summary.txt  (e.g. sd_0__Summary.txt)
//											solution time step index is
//																			timeStepIndex
//--------------------
// ti_DSU_Fragment:							Files that have D (damage), S (stress), U (displacement -> maxDelU, uL, uR) 
//											for all INTERFACEs
//			example file					sd_[subdomainNo]_InterfaceRawFinalSln_tPos_[timeStepIndex], (e.g. sd_0_Interface_DSU_Fragment_tPos_400.txt)
//--------------------
// ti_BulkInterfacePoints:					Files that have displacement, velocity, displacements across in the subdomain
//											for all INTERFACES and BULKS (certain number of points)
//											purpose: 
//													1. These files are used to compute spatial and spacetime integrals in sd_[subdomainNo]__Summary.txt 
//													2. They are used to plot high resolution solutions using the MATLAB function plot_space_spacetime_results.m
//			example file					sd_[subdomainNo]_BulkInterface_tPos_[timeStepIndex].txt (e.g. sd_0_BulkInterface_tPos_4000.txt)
//--------------------
// ti_RawFinalSln_Scalars:					Stores Raw solution at INTERFACES for
//													a. Final solutions (e.g. uL,R  vL,R, sigma, D, Dsrc) in files sd_0_InterfaceRawFinalSln_tPos_400.txt
//													b. Interface scalars (k, D, contact fraction, slip fraction), ... 
//													from the solution stage (not PP)
//			example file					a. sd_[subdomainNo]_InterfaceRawFinalSln_tPos_[timeStepIndex].txt (e.g. sd_0_InterfaceRawFinalSln_tPos_400.txt)
//			example file					b. sd_[subdomainNo]_InterfaceRawScalars_tPos_[timeStepIndex].txt (e.g. sd_0_InterfaceRawScalars_tPos_400.txt)

typedef enum { timeIndexType_none = -1, ti_AllTimeSteps, ti_DSU_Fragment, ti_BulkInterfacePoints, ti_RawFinalSln_Scalars, timeIndexType_SIZE } timeIndexType;
string getName(timeIndexType dat);
void name2Type(string& name, timeIndexType& typeVal);
ostream& operator<<(ostream& out, timeIndexType dat);
istream& operator>>(istream& in, timeIndexType& dat);


// the enumeration RMode_IOT is similar to RiemannMode_StorageTbut mainly for input / output purposes
#if DiM2a3_F
typedef enum { rm_allcs, rm_sp, rm_cst, rm_csl, RMode_IOT_SIZE } RMode_IOT;
#else
typedef enum { rm_allcs, rm_sp, rm_cst, RMode_IOT_SIZE } RMode_IOT;
#endif

string getLatexName(RMode_IOT dat);
string getName(RMode_IOT dat);
bool name2Type(string& name, RMode_IOT& typeVal);
ostream& operator<<(ostream& out, RMode_IOT dat);
istream& operator>>(istream& in, RMode_IOT& dat);

// DamageState_IOT is similar to SLFF_InterfacialDamageModeType but for I/O purposes
//				Hijacked for other uses
//		ds_inactive				->			not referring to statistics of fragmented points ...
//									OR
//											for denominator it means no denominator exists
//					alpha and beta from Contact_Damage_State_Config::Finalize_OneFragmentOneClassificationStatSet
//		ds_alpha
//		ds_beta 
//		ds_tz					->	dealing with things at time zero (e.g. min, max, ...) of the initial strength field
////			Real physical meanings
// nz	: nonzero
// z	: zero
// gzn1 : greater than zero less than 1
// 1	: equal to 1
typedef enum {ds_inactive = -4, ds_alpha = -3, ds_beta = -2, ds_t0 = -1, ds_Dall, ds_Dnz, ds_Dz, ds_Dgzn1, ds_D1, DamageState_IOT_SIZE} DamageState_IOT;

string getLatexName(RMode_IOT dat);
string getName(RMode_IOT dat);
bool name2Type(string& name, RMode_IOT& typeVal);
ostream& operator<<(ostream& out, RMode_IOT dat);
istream& operator>>(istream& in, RMode_IOT& dat);

typedef enum {fio_strength, fio_stress, Field_IOT_SIZE} Field_IOT;

string getLatexName(Field_IOT dat);
string getName(Field_IOT dat);
bool name2Type(string& name, Field_IOT& typeVal);
ostream& operator<<(ostream& out, Field_IOT dat);
istream& operator>>(istream& in, Field_IOT& dat);

// this is a class that goes with things stored in Contact_Damage_State_Config
class FragmentationPtsStamp
{
public:
	FragmentationPtsStamp();
	void MakeVoid();
	// these 4 values refer to what we get from Contact_Damage_State_Config storage
	//	d_numerator_or_flag ==
	//1.	inactive					nothing is obtained from Contact_Damage_State_Config): A trick to turn this off
	//2.	alpha, beta					alpha and beta values computed in Contact_Damage_Stat_Process_All_Interfaces are returned
	//3.	t0							values of the initial mesh are returned
	//4.	Dall, Dnz, Dz, Dgzn1, D1	choosing numerator interface point damage state (allD, 0<D, D=0,0<D<0,D=1)

	// ONLY for case 4 above other entries will have a meaning
	//i.		cs_numerator	STS, ST, S		choosing numerator interface point contact/separation state 
	//											STS=stick + separaton (all points in terms of contact / separation mode)
	//											ST=stick (zero displacement jump)
	//											S=separation (> 0 displacement jump)

	//ii. d_denominator_or_flag:
	//ii.a)						inactive: no denominator exists
	//ii.b)						Dall, Dnz, Dz, Dgzn1, D1 : denominator exists AND cs_denominator is meaningful

	//iii. cs_denominator		ONLY meaningful under ii.b) same values of STS, ST, S

	DamageState_IOT	d_numerator_or_flag;
	DamageState_IOT	d_denominator_or_flag;
	RMode_IOT		cs_numerator;
	RMode_IOT		cs_denominator;
};

// this is the class that stores what time we are referring to in the post-process stage (or actually no time through loadingStagesT_none  options
class OneSubdomainPPS2Data_runInfo;
class DomainPostProcessS2;
class OneTimeValuePPS2Data;
class OneTimeInterfaceFlds_FragmentationPPS2;
class Contact_Damage_State_Config;

class PPS2_TimeStamp
{
	friend istream& operator>>(istream& in, PPS2_TimeStamp& dat);
public:
	PPS2_TimeStamp();
	void MakeVoid_PPS2_TimeStamp();
	void PPS2_TimeStamp_Read(istream& in);

	loadingStagesT timeStampType;

	// for timeStampType == lst_actualTime
	// if actualTime_indexType == timeIndexType_none or actualTime_index < 0 actualTimeval is used, else index + index time provides actual time
	timeIndexType actualTime_indexType;
	double actualTime_val;
	int actualTime_index;

	// for timeStampType == lst_lin OR lst_max OR lst_final
	lstClassificationT lct;
	//	classification type nails down the exact time
};


// with the content of this class we can access various things stored in 
// DomainPostProcessS2::vector<OneSubdomainPPS2Data> onesubdomainPPS2 + fields_initialResolution, fields_finalResolution

class PPS2_dataPointer
{
public:
	PPS2_dataPointer();
	void MakeVoid_PPS2_dataPointer();
	// returns false if have reached the end of the list
	bool PPS2_dataPointer_Read(istream& in, string& buf);

	// if this name matches the key of the map
	//	DomainPostProcessS3::timeStampOverwriters 
	// the timestamp / active / name_CVS / name_Latex will be modified
	string overwriterName;

	// this provides a simple means to turn this on or off when reading from the text file. if off, the entry is ignored
	bool isActive;
	// these two names are what will be used in the plots and in CVS file eventually
	string name_In_CSV_file;
	string name_Latex;

	int subdomainNo; // -1 -> it takes the value of mainSubdomainNo below

	// stat_sv_type specifies the following:
	//	== sso_all_field_val					-> dealing with a vector output:
	//				FA. Initial inhomogeneous field
	//				FB. D, sigma, U at interfaces at later time
	//				FC. Fragment sizes 
	//	else dealing with scalar values wherein stats are taken over
	//				SA. over time for space/spacetime integrals -> it includes time0, timeFinal, mean over time, min over time, ....
	//				SB. over fragment segments when fc != FragmentationCriterionT_none 
	setStatOp_type stat_sv_type;

	// whether referring to a particular time
	PPS2_TimeStamp timeStamp;

	// fldName:
	//			FNA. == field_finalResolution -> returns field_finalResolution[subDomainNo]
	//			FNB. (else), refers to the field name for space/spacetime fields
	
	string fldName;

	// takes the value other than FragmentationCriterionT_none if dealing with a statistics
	FragmentationCriterionT fragmentationCriterion;

	// both of these take _none mode unless brittleness indicators are sought
	BrittlenessIndicatorRatioType brittlenessRatioT;
	BrittlnessIndicatorFieldType brittlenessFieldT;

	//////////////////////////////////////////////////////////

	//////// the following change of the data is done subsequently on the data obtained from variables above
	// step 1. the value may be normalized (e.g. by input energy)
	RunNormalizationQuantT normalizationMode;
	// step 2. log, .... of the value obtained from above may be taken
	valPOperT oprType;

	///////////////////////////////////////////////////////
	// accessing values pertained to Contact_Damage_State_Config
	FragmentationPtsStamp	fragPtStat;
};

class PPS2_dataPointerVec
{
public:
	PPS2_dataPointerVec();
	bool PPS2_dataPointerVec_Read(istream& in);
	inline PPS2_dataPointer& operator[](unsigned int i) { return pointers[i]; };
	inline const PPS2_dataPointer& operator[](unsigned int i) const { return pointers[i]; };
	bool Get_Scalar_Or_Vector_Output(DomainPostProcessS2* configPPS2Ptr, double& scalarVal, vector<double>& vecVal, ScalarFieldOutputT& sfot, OneTimeValuePPS2Data* modifiableOneTimeSlnsPtr, PPS2_TimeStamp* timeStamp4FixedTimePtr = NULL, double temporalFieldTimeStepValOrNumSteps = -400, int spatialFieldResolutionCorrector = 0);

	bool isvActive;
	unsigned int sz_pointers;
	vector<PPS2_dataPointer> pointers;
	vector<string> twoOperands;

	// temporary members
	string name_In_CSV_file;
	string name_Latex;
};

///////////////
// this is a fast way to overwrite the time stamp of PPS2_dataPointer quickly and if needed active / deactivate all its instances
class PPS2_dataPointer_StageOverwritter
{
public:
	PPS2_dataPointer_StageOverwritter();
	void MakeVoid_PPS2_dataPointer_StageOverwritter();
	// returns false if have reached the end of the list
	bool PPS2_dataPointer_StageOverwritter_Read(istream& in);
	void Overwrite_PPS2_dataPointer(PPS2_dataPointer& datPointer);

	int isActive;
	string post_name_In_CSV_file;
	string post_name_Latex;
	PPS2_TimeStamp timeStamp;
};

class DomainPostProcessS2;

class Contact_Damage_State_IO_Stat_1Field
{
public:
	void SetFieldName(const string& field_name, const string& field_nameLatex);
	// rt \in {rm_sp, rm_cst}, dt \in {ds_Dz, ds_Dgzn1, ds_D1}
	void UpdateValue(double value, double xPos, RMode_IOT rt, DamageState_IOT dt);
	// level 0: stat group, level 1: stat subgroup, level 2: name, level 3, name latex, level 4 values
	// for levels 0 to 3 prints, it prints time index and value, alpha, beta heading
	// for level 4 only prints alpha and beta
	void PrintLineLevel(ostream& out, unsigned int nodeCount, unsigned int level, bool print_sdiv, bool print_min_max, bool print_min_max_loc);

	string name_fld_fragT, nameLatex_fld_fragT;
	statHolder stats1Field[RMode_IOT_SIZE][DamageState_IOT_SIZE];
	// Levy_Molinari_2010_Dynamic fragmentation of ceramics, signature of defect sand scaling.pdf
	double alpha; // eqn (21), Nbroken/NodeCount
	double beta; // eqn (23) broken@time_t(mean - min)/@time_zero(mean - min) - division only occurs for strength values
	// I take a brute-force approach and do not merge the stats later, rather compute union ones directly
	// static members
	static unsigned int sz_r_ind_to_appl[RMode_IOT_SIZE];
	static vector<RMode_IOT> r_ind_to_appl[RMode_IOT_SIZE];

	static unsigned int sz_d_ind_to_appl[DamageState_IOT_SIZE];
	static vector<DamageState_IOT> d_ind_to_appl[DamageState_IOT_SIZE];
	static void SetStatics_Contact_Damage_State_IO_Stat_1Field();
};

class Contact_Damage_State_IO_Stat_AllsField
{
public:
	void Read_Contact_Damage_OneLineData(istream& in, Field_IOT fit, Contact_Damage_State_Config* contact_damage_confPtr);
	loadingStagesT lsi;
	double time;
	int timeIndex;

	Contact_Damage_State_IO_Stat_1Field stat_fields[Field_IOT_SIZE];
	vector<double> vals_strengthOfD1Bonds;
	vector<unsigned int> inds_strengthOfD1Bonds;
};

// for all fields and all times/stages
class Contact_Damage_State_IO_Stat_AllsField_Times
{
public:
	FragmentationCriterionT fc;
	lstClassificationT ct;
	Contact_Damage_State_IO_Stat_AllsField_Times();
	~Contact_Damage_State_IO_Stat_AllsField_Times();

	vector<Contact_Damage_State_IO_Stat_AllsField*> statsPtr;
	unsigned int statsPtr_curIndex;
	fstream* stat_outPtr[Field_IOT_SIZE];
	fstream strength_out;
};

class Contact_Damage_State_Config
{
public:
	Contact_Damage_State_Config();

	///////////////////// computation functions
	void Initialize_Config_ForWholeSubdomain(OneSubdomainPPS2Data_runInfo& segmentInfo, unsigned int subdomainNoIn);
	void Initialize_OneFragmentOneClassificationStatSet(Contact_Damage_State_IO_Stat_AllsField_Times& stats, FragmentationCriterionT fc, lstClassificationT ct);
	void Contact_Damage_Stat_Process_All_Interfaces(Contact_Damage_State_IO_Stat_AllsField_Times& stats, OneTimeInterfaceFlds_FragmentationPPS2* fragDatPtr, loadingStagesT lsi, double time, unsigned int timeIndex);

	///////////////////// retrieve value function
	// returns true if successful
	double GetValue(PPS2_dataPointer& datPointer, bool& success);
	///////////////////// reading data
	void Read_Contact_Damage_State_Config(DomainPostProcessS2* configPPS2);


	// configuration paras
	bool isActive;
	bool processField[Field_IOT_SIZE];
	bool process_ciriticalPoints, process_all_times;
	bool printStrengthD1_ciritcalPoints, printStrengthD1_all_times;
	bool print_sdiv, print_min_max, print_min_max_loc;
	
	/////////////////////////////////// set based on each individual run
	unsigned int nodeCount;
	/// calculated from initial strength stats
	// is: initial strength
	// min, max, mean, and sdiv are calculated from initial strength field
	// is_beta_denom = 1.0 / (is_mean - is_min); cf. eq (23) of Levy_Molinari_2010_Dynamic fragmentation of ceramics, signature of defect sand scaling.pdf
	double is_min, is_min_x, is_max, is_max_x, is_mean, is_sdiv, is_cov, is_beta_denom;
	unsigned int subdomainNo;
	// data computed
	Contact_Damage_State_IO_Stat_AllsField_Times stats_frag_class[FragmentationCriterionT_SIZE][lstClassificationT_SIZE];
	Contact_Damage_State_IO_Stat_AllsField_Times stats_frag_allTimes[FragmentationCriterionT_SIZE];

private:
	OneSubdomainPPS2Data_runInfo* segmentInfoPtr;
	// rt \in {rm_sp, rm_cst}, dt \in {ds_Dz, ds_Dgzn1, ds_D1}
	void UpdateValues(Contact_Damage_State_IO_Stat_AllsField_Times& stats, double strength, double stress, unsigned int xIndex, double xPos, RMode_IOT rt, DamageState_IOT dt);
	void StartOneTime_or_Stage_Computation(Contact_Damage_State_IO_Stat_AllsField_Times& stats, loadingStagesT lsi, double time, int timeIndex);
	void Update_Actual_Fragment_points(Contact_Damage_State_IO_Stat_AllsField_Times& stats, OneTimeInterfaceFlds_FragmentationPPS2* fragDatPtr);
	void Finalize_OneFragmentOneClassificationStatSet(Contact_Damage_State_IO_Stat_AllsField_Times& stats);
};


// this class can read files with 4 header lines from file nameWOExtIn.txt (e.g. _sd_0_summary.txt) and write 
// (min, max, mean, sdt, cov, var) 
// in rows of file nameWOExtInstat.txt (e.g. _sd_0_summarystat.txt)
class DataW4LineHeader
{
public:
	DataW4LineHeader();
	void Initialize_DataW4LineHeader(string& nameIn, DataW4LineHeader& statOut, string& nameOut);
	void Compute_Stat(DataW4LineHeader& statSummary);
	void DataW4LineHeader_Write(ostream& out);
	// true if the file is read
	bool DataW4LineHeader_ReadFileName(string& fileNameIn);

	//////////
	// Inputs
	// fldNo: column number for which max and terminal values are obtained
	// forceTerminalvalue2ProvidedValue:
	//						true (suitable for stress and internal energy for example)
	//							it looks when the value is close to providedTerminalValue (e.g. ZERO)
	//						false (suitable for total physical dissipation, e.g. values with nonzero terminal values)
	//							it looks when the value is close to final value (final value of the field in time)
	// forwardFinalValueSearch:			
	//						true:
	//							From the max position it goes forward and the very first point for which the value is close to the terminal or zero (based on the flag above) value it stops
	//						false
	//							It starts from the last value (e.g. corresponding to last time) and goes backward until the very first point is found not close enougth to terminal / zero value

	//	tol4FinalValCheck
	//						> 0 absTol = tol4FinalValCheck				
	//						< 0 absTol = -tolFinalvalCheck * maxAbsVal
	//		absTol is used to check closeness to terminal / zero value for finalVal
	//////////
	// Outputs: If indices < 0, the value is not found

	// maxAbs:
	//		bValuesPositive if maxAbsValue is positive
	//		index_maxAbsVal: index (point index) of this max value
	//		the actual maximum abs value
	// finalVal
	//		index_finalVal: index for this
	//		finalVal: value for this
	//		checkCrossing4FinalValue:			for sigma's it looks for the first time that the sign of sigma changes pass the peak value
	void GetMaxAbs_FinalValue(const DataW4LineHeader& timeSequenceSummaryStat, unsigned int fldNo, bool classification_usePositive4MaxValue, bool forceTerminalvalue2ProvidedValue, double providedTerminalValue, bool forwardFinalValueSearch, bool checkCrossing4FinalValue, double tol4FinalValCheck, bool &bValuesPositive, int& index_maxAbsVal, double& maxAbsVal, int& index_finalVal, double& finalVal);

	// to find the first time D becomes 1
	bool GetFirstLocationGreaterThanValue(const DataW4LineHeader& timeSequenceSummaryStat, unsigned int fldNo, double providedTerminalValue, double tol4FinalValCheck, int& index);


	// this is a data that is zero all early time and at a time that a nonvalue is experienced, nonlinear response is assumed
	// example energy dissipations, loss energy, damage parameters
	//	tol4FinalValCheck
	//						> 0 absTol = tol4FinalValCheck				
	//						< 0 absTol = -tolFinalvalCheck * maxAbsVal of the field
	//		absTol is used to check closeness to zero from the beginning
	void GetLinValue(const DataW4LineHeader& timeSequenceSummaryStat, unsigned int fldNo, double tol4ZeroValCheck, int& index_LinVal);

	// sortingFld is generally 0 (e.g. time)
	// for a given value (e.g. time), it finds the closes sorted value (e.g. times in this class) through its index closestSortingIndex and valueclosestSortingValOut
	void GetIndexClosestSortingVal(unsigned int sortingFldNo, double sortingValIn, unsigned int& closestSortingIndex, double& closestSortingValOut);
	// input
	string nameWOExt;

	// computed
	vector<vector<double> > data_vals;
	vector<string> headerStrs1, headerStrs2, headerStrs3, headerStrs4;
	string header1, header2, header3, header4;
	unsigned int numFlds, numDataPoints;
};

// this class stores fragmentation data based on one of the criteria (D, maxDelU, delU)
class OneTimeOneCriterionFragmentationRawData
{
public:

	OneTimeOneCriterionFragmentationRawData();
	void ComputeFragmentStatsFromFragmentSizes();
	// fldName	= (num, mean, min, max, sdiv, cov) scalarVal = scalar	- sfot = sfo_scalar
	//			= (sizes) -> vecVal = fragment_lengths					- sfot != sfo_scalar
	// returns true is data is valid and computed
	bool Get_Scalar_Or_Vector_Output(setStatOp_type statType, double& scalarVal, vector<double>& vecVal, ScalarFieldOutputT& sfot);

	void Fragmentation_Sizes_Write(ostream& out, int timeIndex, double timeVal);
	void Fragmentation_Sizes_Read(istream& in);
	void OneTimeOneCriterionFragmentationRawData_Write_Data(ostream& out, int timeIndex, double timeVal);
	void OneTimeOneCriterionFragmentationRawData_Read_Data(istream& in);
	static void OneTimeOneCriterionFragmentationRawData_Write_Header(ostream& out);
	static void OneTimeOneCriterionFragmentationRawData_Read_Header(istream& in);
	static void OneTimeOneCriterionFragmentationRawData_Write_Data_void(ostream& out, int timeIndex, double timeVal);

	unsigned int numFragments;
	double mean_fragmentSize, min_fragmentSize, max_fragmentSize, sdiv_fragmentSize, cov_fragmentSize;
	vector<unsigned int> fragmented_interface_indices;
	vector<double> fragment_lengths;
private:
	static void OneTimeOneCriterionFragmentationRawData_Get_Header_Labels(vector<HeaderLabels>& hLabels);
};

// this opens files sd_[subdomainNo]_tPos_[timeIndex]_Interface_DSU_Fragment.txt and reads it
class OneTimeInterfaceFlds_FragmentationPPS2
{
public:
	void Read_OneTimeInterfaceFlds_FragmentationPPS2(istream& in, unsigned int numInterfaces);
	void Empty_DUSigma(unsigned int numInterfaces);

	// names are "D", "maxEffDelU", "sigma", "DelU", "uL", "uR"
	// returns true if the vector is size > 0
	bool Get_Vector_SolutionField(const OneSubdomainPPS2Data_runInfo& segmentInfo, const string& fldName, unsigned int defaultDir, vector<double>& vecVal);
	vector<double> D;
	vector<double> maxEffDelU;
	vector<double> uL[DiM], uR[DiM], sigma[DiM];
#if DSU_PRINT_VS
	vector<double> vL[DiM], vR[DiM];
#endif

	OneTimeOneCriterionFragmentationRawData fragmentationDat[FragmentationCriterionT_SIZE];
};

// this class opens file _sd_[subdomainNumber]_keyParameters.txt and reads its content
class OneSubdomainPPS2Data_runInfo
{
	friend istream& operator>>(istream& in, OneSubdomainPPS2Data_runInfo& dat);
public:
	OneSubdomainPPS2Data_runInfo();
	bool isPeriodic;
	double maxTime, timeStep;
	unsigned int totalTimeSteps, numTimeStep_Interface_DSU_Fragment_Print_4PP, numTimeStep_BulkInterfacePoints_Print_4PP, numSpatialSubsegments_BulkInterfacePoints_Print_4PP, numTimeStep_InterfaceRawFinalSlnScalars_AllSpace_Print_4PP;

	vector<double> lengths;
	vector<double> deltaCs;
	vector<double> sigmaCs;
	double length, xm, xM;
	unsigned int numSegments, bulk_cn_st, bulk_cn_en;
	vector<double> xs;

	double loadTimeScale;
	double EScale;
	double rhoScale;
	double dampingScale;

	double sigmaCScale;
	double deltaCScale;
	double energyCScale;

	// computed (c: wave speed, Z impedance, tau time, psi: strain/stress energy)
	double cScale, ZScale, tauCScale, strainCScale, velCScale, psiCScale;

	unsigned int interface_offset;
	unsigned int numInterfaces;
	double timeStep_BulkInterfacePoints_Print;
	unsigned int maxIndex_BulkInterfacePoints_Print;
	double timeStep_Interface_DSU_Fragment_Print;
	unsigned int maxIndex_Interface_DSU_Fragment_Print;

	/// Auxiliary function
		// getting the index of "OneTimeInterfaceFlds_FragmentationPPS2_TimeIndex" from time value
	int get_Interface_DSU_Fragment_TimeIndexValue(double timeValue, double& mod_timeValue);
	int get_Interface_BulkInterfacePoints_TimeIndexValue(double timeValue, double& mod_timeValue);
	// return true if successful
	// for stage files, we copy these files from source (output of the run) to PPS2 folder so if the run folder is deleted, data is available
	// readFromRunOutputFolder:				true (for stage files, when first time the content is read) false (rerunning PPS2 -> it reads from the copied file to PPS2 folder)
	// copyRunOutputFile2PPS2Folder: copies the file from the source (output) to target (PPS2) folder
	bool Read_OneTime_Interface_DSU_Fragment(unsigned int subdomainNo, double timeValue, OneTimeInterfaceFlds_FragmentationPPS2& fieldFragInfo, bool readFromRunOutputFolder, bool copyRunOutputFile2PPS2Folder, loadingStagesT lsi = loadingStagesT_none, lstClassificationT ci = lstClassificationT_none);
	void ComputeAllScales();
};
// this includes data for one classification, e.g. bulk or bc strain/stress response, damage values, etc.. 
// OneSubdomainPPS2Data includes lstClassificationT_SIZE of them

string getBrittleness_IndicatorName(bool isLatex, BrittlnessIndicatorFieldType bf, BrittlenessIndicatorRatioType br);


class OneTimeValuePPS2Data
{
public:
	OneTimeValuePPS2Data();
	void SetAllTimeValsIndices_FromActualTime(OneSubdomainPPS2Data_runInfo segmentInfo);
	//	unsigned int subdomainNo;
//	lstClassificationT ct;
//	loadingStagesT lst;

	int timeIndex_4Stage;
	double timeValue_4Stage;
	// indices and time values with medium accuracy -> fields D, sigma, fragmentation data from sd_[subdomainNo]_tPos_[timeIndex]_Interface_DSU_Fragment.txt
	int timeIndex_DSU_Fragment4Stage;
	double timeValue_DSU_Fragment4Stage;
	// indices and time values with the coarsest accuracy -> fine (with in segment point) detail values of u, v, sigma, D, ... used in sd_[subdomainNo]_tPos_[timeIndex]plot_03_[fldNo]_[fldName].png / fig files
	int timeIndex_BulkInterfacePoints4Stage;
	double timeValue_BulkInterfacePoints4Stage;

	vector<double> space_spacetime_integrals;
	OneTimeInterfaceFlds_FragmentationPPS2 fragmentation4Stages;
};

class OneClassificationPPS2Data
{
public:
	OneClassificationPPS2Data();
	void BrittlenessIndicators_Write(ostream& out);
	void BrittlenessIndicators_Read(istream& in);
//	unsigned int subdomainNo;
//	lstClassificationT ct;
	double brittlenessIndicators[BrittlnessIndicatorFieldType_SIZE][BrittlenessIndicatorRatioType_SIZE];
	OneTimeValuePPS2Data data4Stages[loadingStagesT_SIZE];
};

// this class stores Stage2 (after simulation) post-process data for 
// ONE subdomain
class OneSubdomainPPS2Data
{
public:
	OneSubdomainPPS2Data();
	void Compute_OneSubdomain_PPS2();
	void StageRelatedData_PPS2_Write();
	// returns true if stats can be read
	bool StageRelatedData_PPS2_Read(int subdomainNoIn, DomainPostProcessS2* configPPS2In);

	// For the cases that the SPECIFIC TIME IS GIVEN (timeStamp.timeStampType == lst_actualTime)
	// this function calculates oneTimeSlns if needed (e.g. IsItAlreadySetUp4_ActualTimeType() returns false)
	// extract_space_spacetimeIntegrals: these are the entries from the file sd_[subdomainNo]__Summary.txt for the given time
	// extract_DSUFields: reads fields D, S, U for this time
	// compute_FragmentationStats: it calculates fragmentation statistics
	// return value = true, if calculation is successful
	bool Compute_OneTimeValuePPS2Data_4ActualTime(PPS2_TimeStamp  timeStamp, OneTimeValuePPS2Data& oneTimeSlns, bool extract_space_spacetimeIntegrals, bool extract_DSUFields, bool compute_FragmentationStats);

	// when timeStamp in dataPointer points to a actual time (as opposed to stage calculations) modifiableOneTimeSlns is created (if needed) and reused for next computations
	// return value is true if the computed value is set correctly (e.g. the stage exists, ...)
	bool Get_Scalar_Or_Vector_Output(PPS2_dataPointer& datPointer, double& scalarVal, vector<double>& vecVal, ScalarFieldOutputT& sfot, OneTimeValuePPS2Data& modifiableOneTimeSlns, double temporalFieldTimeStepValOrNumSteps = -400);

//	computed
	unsigned subdomainNo;
	OneSubdomainPPS2Data_runInfo segmentInfo;
	DataW4LineHeader timeSequenceSummary, timeSequenceSummaryStat;
	OneClassificationPPS2Data diffClassifications[lstClassificationT_SIZE];

	DomainPostProcessS2* configPPS2;

private:
	/////////////////// stage x [lin, max, final] computations
	// ci can only be ct_epssig_bc, ct_epssig_bulk
	void Set_All_Stage_Indices_Times();
		void Set_Time_Stage_Indices_Times();
		void Set_EpsSigma_Stage_Indices_Times(lstClassificationT ci);
		void Set_Damage_Stage_Indices_Times();
		void Set_PhysInterfaceDiss_Stage_Indices_Times();
		void Set_U_Stage_Indices_Times();
		void Set_lin_max_final_generic_data(lstClassificationT ci, int time_index_lin, int time_index_max, int time_index_final);
		void Finalize_OneStage(lstClassificationT ci);

		void Compute_BrittlenessIndices(lstClassificationT ci, int fldStrn = -1, int fldEne = -1);
	/////////////////// 
	// fragmentation analysis:
	void ComputeAllFramentationAlaysis();
		void Do_AllFragmentationAnalysis(OneTimeInterfaceFlds_FragmentationPPS2& fragDat);

	// fragCriterion: it's one of D, maxDelU, or DelU that will be used to decide framentation
	// criterion						factor2DecideFragmented meaning (when an interface is deemed to be broken)
	// fct_Damage						D >= factor2DecideFragmented
	// fct_maxEffDelU, fct_DelU			(max)delU >= deltaC * factor2DecideFragmented
		void Compute_FragmentationSizeData_Step1_DetermineBreakingPoints(FragmentationCriterionT fragCriterion, double factor2DecideFragmented, OneTimeInterfaceFlds_FragmentationPPS2& fragDat);
		void Compute_FragmentationSizeData_Step2_DetermineFragments(FragmentationCriterionT fragCriterion, OneTimeInterfaceFlds_FragmentationPPS2& fragDat);

	void Ensure_openning_complete_time_history_space_spacetime_integrals_summary_file();
};

// manager for second stage (i.e. after the solution) post-process
class DomainPostProcessS2
{
public:
	DomainPostProcessS2();
	void DomainPostProcessS2_Read_WO_Initialization(istream& in);

	void Main_DomainPostProcessS2();
	// for spatial fields over x, we may want to change the resolution. For example, if values given at 1001 points (resolution 1000), but we want the new resolution of 100 (101 points)
	// spatialFieldResolutionCorrector = 100 or spatialFieldResolutionCorrector = -10 (i.e. 10 fold reduction)
	// note that the solution can only decrease (from a fine mesh to a coarser representation)
	// spatialFieldResolutionCorrector == 0 -> does not change the field
	// for temporal data how many time steps should be taken from time 0 to final time
	bool Get_Scalar_Or_Vector_Output(PPS2_dataPointer& datPointer, double& scalarVal, vector<double>& vecVal, ScalarFieldOutputT& sfot, OneTimeValuePPS2Data& modifiableOneTimeSlns, double temporalFieldTimeStepValOrNumSteps = -400, int spatialFieldResolutionCorrector = 0);

	// return field index in summary classes (space + spacetime integrals)
	// -1 if the string is not found
	int get_space_spacetime_integral_index_from_name(const string& fldName);

	// INPUTS
	// if want to ensure redoing PPS2 calculation (e.g. if limit of D for fragmentation is changed) set the boolean below to true
	bool redoPPS2Calculations;
	///////////////////////////////////////////////////////////////
	////////		classification + stage calculation
	// the booleans below can mutually be true
	// if need to calculate lin, max, and final stages, set the boolean below to true
	bool doCalculateStages;
	/// classification flags

	bool classification_Active[lstClassificationT_SIZE];
	double classification_tol4LinValChecks[lstClassificationT_SIZE];
	double classification_tol4FinalValChecks[lstClassificationT_SIZE];
	bool classification_forwardFinalValueSearches[lstClassificationT_SIZE];
	// for maximum value if below is true -> value is used for computing maximum, if -value is used (i.e. minimum value is extracted)
	bool classification_usePositive4MaxValues[lstClassificationT_SIZE];
	/// for time stage
	// if values < 0, for final -> final time is used, for lin -> when energy input stbilizes, for mid nothing is done
	double actualTimesProvided_4TimeStage[loadingStagesT_SIZE];

	///////////////////////////////////////////////////////////////
	////////		fragmentation
	// if true: for all stages (lin, max, final) x (strain, ...., U) it also calculates fragmentation data
	bool doFragmentation_DSigmaFieldExtraction_4_Stages;
	// if true, calculates fragmentation analysis and D and sigma field extraction for all available times
	bool doFragmentation_DSigmaFieldExtraction_4_AllTimes;
	// if the above is true and below is true, the segment lengths are written. If only above is true, only the brief stat (mean, number, ...) of fragmentation sizes will be written
	bool doFragmentation_DSigmaFieldExtraction_4_AllTimes_IncludedDetailedFragmentSizes;
	// values < 0 mean that such options are invalid
	double factors2DecideFragmented[FragmentationCriterionT_SIZE];

	Contact_Damage_State_Config contact_damage_conf;
#if 0 // removed as they are not really necessary
	double tol4FinalValCheck_energyInput; // this is used for lin stage of classification time if actualTimesProvided_4TimeStage[lin] < 0

	// if need to output solution at actual time intervals (e.g. 20 intervals for the simulation) turn this to on
	bool doCalculateActualTimeIntervals;
	// bActualTime: coresponds to actual time output results (boolean above) or whenever a actual time is going to "this" to calculate results
	// under these conditions (as opposed to stage-based calculation) three booleans below specify the types of computation needed
	bool bActualTime_extract_space_spacetimeIntegrals;
	bool bActualTime_extract_DSUFields;
	bool bActualTime_compute_FragmentationStats;
#endif


	///// CALCULATED from ABOVE in 	Main_DomainPostProcessS2
	/////////////////////////////////////////////////////////////// indices for summary file
	vector<bool> subDomainActive4PPS2;
	unsigned int numSubdomains;
	vector<OneSubdomainPPS2Data> onesubdomainPPS2;

	unsigned int sind_time;
	unsigned int sind_phi, sind_K, sind_U, sind_U2phi, sind_EneInp, sind_EneL, sind_EneR, sind_EneBC, sind_phi0, sind_diss_tot,
		sind_diss_lost, sind_eneIDiss, sind_diss_interface_lost, sind_energy_phys_diss_2_input, 
		sind_energy_phys_diss_lost_2_input, sind_energyIDiss_Recoverable_2_input,
		sind_energy_numerical_diss_2_input, sind_diss_recov, sind_dissPower_interface_lost;
	unsigned int sind_eps_bc, sind_eps_bulk_intfc, sind_sig_bc, sind_phi_tot_bc, sind_phi_diss_lost_bc, sind_phi_recov_bc, sind_phi_diss2tot_bc, sind_phi_dissL2phid_bc;;
	unsigned int sind_eps_bulk, sind_sig_bulk, sind_phi_tot_bulk, sind_phi_diss_lost_bulk, sind_phi_recov_bulk, sind_phi_diss2tot_bulk, sind_phi_dissL2phid_bulk;
	unsigned int sind_Dbar, sind_Dmax, sind_Dmin, sind_DsrcMax, sind_P;
	unsigned int defaultDir;
	string defaultDirStr;

	double normalizations[RunNormalizationQuantT_SIZE];
	lstClassificationT default_ClassificationType;

	int mainSubdomainNo;
private:
	// index 0: Field position (can be for subdomain or the component that is random, e.g. (E, rho, ...) - these are field0_initialResolution , ..., in the run output folder
	// index 1: spatial position index, e.g. for a domain of size 1, resolution 100, index 50 is the middle of the domain
	vector< vector <double> > fields_initialResolution, fields_finalResolution;
	// the post-process S2 stage results are already calculated and stored in PPS2 folder, so they don't need to be recalculated
	bool stageSolutionsExist;
	// returns true if any subdomains are found
	void Initialize_DomainPostProcessS2();
	bool Read_SummaryData();
	void Compute_AllSubdomain_PPS2();
};

class SlnPP2FileMover
{
public:
	SlnPP2FileMover();
	bool SlnPP2FileMover_Read(istream& in, string& buf);
	bool SlnPP2FileMover_MoveFile();
	bool isActive;
	bool isPPS2;	// false file taken from solution file, otherwise from PPS2 file
	string fileName;
	string ext;
};

void GetSubdomainIndexed_TimeIndexed_PostProcess_FileName(string& fileName, unsigned int subdomainNo, int timeIndex, const string& specificName,
	loadingStagesT lsi = loadingStagesT_none, lstClassificationT ci = lstClassificationT_none, FragmentationCriterionT fc = FragmentationCriterionT_none, string ext = "txt", bool addUnderline = false);

void MAIN_DomainPostProcessS2(string fileName = "config/Domain/PPS2ConfigTest.txt");

#endif