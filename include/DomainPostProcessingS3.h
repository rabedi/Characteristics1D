#ifndef DOMAIN_POST_PROCESSING_S3__H
#define DOMAIN_POST_PROCESSING_S3__H

#include "DomainPostProcessingS2.h"

//////////////////////////////////////////
// This file is Level 3
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
typedef enum {pps3_om_scalars_vectors, pss3_om_scalars, pps3_om_vectors, pps3_om_T_SIZE} pps3_om_T;

string getName(pps3_om_T dat);
void name2Type(string& name, pps3_om_T& typeVal);
ostream& operator<<(ostream& out, pps3_om_T dat);
istream& operator>>(istream& in, pps3_om_T& dat);
string getAddedName4OutputFile(pps3_om_T dat);

// this break data to input and output groups (for example initial damage profile is considered input)
typedef enum {pps3_i, pps3_o, PPS3_IOT_SIZE} PPS3_IOT;

string getName(PPS3_IOT dat);
void name2Type(string& name, PPS3_IOT& typeVal);
ostream& operator<<(ostream& out, PPS3_IOT dat);
istream& operator>>(istream& in, PPS3_IOT& dat);

// output type for PPS3:
//		pps3_stageStat: uses post-process classification time x stages (lin, max, final) - it also includes stats such as min/max/etc of data over time
//		pps3_timeStep:	for each timeStep (generally with DSU output data) a separate CSV file is written
typedef enum {pps3_stageStat, pps3_timeStep, PPS3_TimeOT_SIZE} PPS3_TimeOT;

string getName(PPS3_TimeOT dat);
void name2Type(string& name, PPS3_TimeOT& typeVal);
ostream& operator<<(ostream& out, PPS3_TimeOT dat);
istream& operator>>(istream& in, PPS3_TimeOT& dat);

class DomainPostProcessS3
{
public:
	DomainPostProcessS3();
	void MAIN_DomainPostProcessS3(const string configFileName = "config/Domain/PPS3ConfigTest.txt");

	DomainPostProcessS2 configPPS2;

	/////////////////////////////////////////////////////////
	/// version controls (changes to delc, a, lc, .... through sfcm and sfcm_gen
	// true: for each version there is a separate folder for all random fields of One version
	// false: the output folder is always _PP3 and contains all versions and all serial numbers
	bool version_seperatePP3Folders;
	// true: values that are not computed (e.g. values when phid is stabilized) are printed as nan rather than 1e40
	bool print_uncomputed_vals_as_nan;
	// 0: version columns are not added to PP3 file
	// 1: they are added
	// 2: only columns of fields that have more than 1 values (e.g. more than 1 cor. length, ...) are written
	int version_print_version_columns;
	// true: writes version material before acceptance and serial number
	// false: these are written after accept flag and serial number
	bool version_print_before_accept_serial;
	// prints version Number from the list of versions in generator file
	bool version_print_version_No;
	// prints version Number for the output folder. These are version numbers from generator file, plus possible offset (SolveParameters::versionOffset)
	bool version_print_version_No_wOffset;
	// print indices of version columns (e.g. for index 5 correlation length -> 5 is printed)
	bool version_print_indices;
	// print values of version columns (e.g. for index 5 correlation length -> correlation[5] is printed)
	bool version_print_values;

	/////////////////////////////////////////////////////////
	// print, scalars + vectors, scalars, vectors, options are added in this file
	vector<pps3_om_T> output_modes;

	// out: output file
	string root;
	string out_baseName;
	string out_ext; // generally it will be csv for csv files, but xlsx and txt are some other options
	// 0 or 1 clear, -1 is an active only in large memory access case
	int outputTypeActive[PPS3_TimeOT_SIZE];
	vector<PPS2_dataPointer> dataPointers[PPS3_TimeOT_SIZE][PPS3_IOT_SIZE];
	// for spatial fields over x, we may want to change the resolution. For example, if values given at 1001 points (resolution 1000), but we want the new resolution of 100 (101 points)
	// spatialFieldResolutionCorrector = 100 or spatialFieldResolutionCorrector = -10 (i.e. 10 fold reduction)
	// note that the solution can only decrease (from a fine mesh to a coarser representation)
	// spatialFieldResolutionCorrector == 0 -> does not change the field
	int spatialFieldResolutionCorrector;

	// for temporal fields (e.g. energyDissipated, sigma_homogenized, etc.) there are two options:
	// temporalFieldTimeStepValOrNumSteps > 0 -> TimeStepVal = temporalFieldTimeStepValOrNumSteps
	//									  < 0 -> NumTimeSteps = -temporalFieldTimeStepValOrNumSteps (e.g. final time = 10, this value = -400 -> timeStep = 10/400)
	double temporalFieldTimeStepValOrNumSteps;
	// for stage type data, we can get invalid data (for example if the run does not reach the terminal stage)
	// in this case we generally don't print the data for this file. The values for such data are invalidNum (1e40)
	// If the boolean below is true, the values are written as invalidNum
	bool addInvalidData;
	// for time outputs, we can generate outputs for all DSU stems in each case timeStep4_DSU_outputs == 1, ever other of such time value for == 2, every third for this value == 3, ...
	int	 timeStep4_DSU_outputs;

	map<string, PPS2_dataPointer_StageOverwritter> timeStampOverwriters;
private:
	bool DomainPostProcessS3_Read_WO_Initialization(istream& in);
	void Get_Version_Header_Names(vector<string>& names, vector<string>& namesLatex);
	void Print_Version_Values(ostream& out);

	// returns true if data is added
	// time_outputIndex:	for tot == pps3_timeStep  -> file out_baseName_time_outputIndex.file_out_ext prints values for 
	//				timeIndex_DSU_Fragment4Stage = time_outputIndex * timeStep4_DSU_outputs
	// accebleOverAlVals stores whether this data has any invalid values (e.g. a stage that does not exist in the solution)
	// return value: true if the point is printed to the file (can be true even if accebleOverAlVals when addInvalidData == true 
	bool ComputePrint_Data(PPS3_TimeOT tot, bool& accebleOverAlVals, int time_outputIndex = -1, pps3_om_T output_mode = pps3_om_scalars_vectors);
	string sep;
	string scalarFieldAddOn[ScalarFieldOutputT_SIZE], scalarFieldAddOn_Latex[ScalarFieldOutputT_SIZE];
};

void MAIN_DomainPostProcessS3(const string configFileName = "config/Domain/PPS3ConfigTest.txt");

#endif