#ifndef GLOBAL_TYPE_CLASSES__H
#define GLOBAL_TYPE_CLASSES__H

#include "Dims.h"

// for 1D domains interfaces are left .... twosided ... right. 
// periodic BC is twoSided too
// ilt_rightWPeriodicShift is the right side of the very last interface which is basically the left side of the very first interface (that is not included in data)
typedef enum {ilt_noSided = -1, ilt_twoSided, ilt_left, ilt_right, interfaceLocation1DT_SIZE} InterfaceLocation1DT;
// characteristics is the case that incoming characteristic (zero for transmitting) from outside is given
// bct_Symmetric: for a given direction, the characteristics are equal from opposite side
// bct_AntiSymmetric: for a given direction, the characteristics are opposites from opposite side
//		For example in 2D for mode I , we need bct_Symmetric for dir 0 and bct_AntiSymmetric for dir 1 (plus dir 2 for 3D)
//						  for mode II, we need bct_AntiSymmetric for dir 0 and bct_Symmetric for dir 1
	typedef enum {bct_Undecided = -1, bct_Dirichlet, bct_Neumann, bct_Characteristics, bct_Symmetric, bct_AntiSymmetric, bct_PeriodicOrBloch, bct_Unspecified, BoundaryConditionT_SIZE} BoundaryConditionT;

string getName(BoundaryConditionT dat);
void name2Type(string& name, BoundaryConditionT& typeVal);
ostream& operator<<(ostream& out, BoundaryConditionT dat);
istream& operator>>(istream& in, BoundaryConditionT& dat);

typedef enum { allOffV, allOnV, partialV, onOffTSIZE } onOffT;
typedef enum { ect_notActive, ect_nonDimensional, ect_Dimensional, errorCheckT_SIZE } errorCheckT;

string getName(errorCheckT dat);
void name2Type(string& name, errorCheckT& typeVal);
ostream& operator<<(ostream& out, errorCheckT dat);
istream& operator>>(istream& in, errorCheckT& dat);

typedef enum { iof_none, iof_ascii, iof_binary, IOF_type_SIZE } IOF_type;

string getName(IOF_type dat);
void name2Type(string& name, IOF_type& typeVal);
ostream& operator<<(ostream& out, IOF_type dat);
istream& operator>>(istream& in, IOF_type& dat);
bool getExt(IOF_type dat, string& ext);

// the letters are the names of these fields
/// Domain ones
// sp		so_domain_sp: solves the domain and PPS3 it (PPS3 internally calls PPS2)
// s		so_domain_s: only solves the domain
// p		so_domain_p: only PPS3 the domain (and as needed does PPS2)
// p2		so_domain_p2: only does PPS2 on the domain
/// Interface ones
// i		so_interface_s: solves one interface for all times
/// One point
// o		so_onePoint_s: solves one point only
/// Random fields
// eff		so_elfrac_fields: tests 1 or a number of random fields that represent elastic/fracture problems
// f		so_one_field: tests only 1 random field
typedef enum { so_domain_sp, so_domain_s, so_domain_p, so_domain_p2, so_interface_s, so_onePoint_s, so_elfrac_fields, so_one_field, 
so_configGen, so_configRead, so_wnRandomFieldGen, solveOptions_SIZE } solveOptions;

string getName(solveOptions dat);
void name2Type(string& name, solveOptions& typeVal);
ostream& operator<<(ostream& out, solveOptions dat);
istream& operator>>(istream& in, solveOptions& dat);


// a_terminate_run_prematurely is added so when error is still too large but time step is becoming increasingly small, a_terminate_run_prematurely is issued to stop the simulation
// other terminate is when the run is ready to exit due to reaching final time, ...
enum AdaptivityF { a_unassigned, a_coarsen, a_none, a_refine, a_terminate_run_correctly, a_terminate_run_prematurely};

class AdaptivityS
{
	friend ostream& operator<<(ostream& out, const AdaptivityS& dat);
public:
	AdaptivityS();
	// if not worse, not updated
	void Update_AdaptFlag(AdaptivityF a_flagIn);
	AdaptivityF get_a_flag() const { return a_flag; }
	double a_delt;
private:
	AdaptivityF a_flag;
};

class HeaderLabels
{
public:
	// example phi = U + K = energy at time t:
	// group = "energy", subgroup = "tT" (t = T), textLabel = "energy", latexLabel = "\phi"
	string group;
	string subgroup;
	string textLabel;
	string latexLabel;
};

void PrintHeader(const vector<HeaderLabels>& hLabels, ostream& out, bool isCommaSeparated = false);

template <class T>
inline int Find(vector<T>& vec, const T &value)
{
	for (unsigned int i = 0; i < vec.size(); ++i)
		if (vec[i] == value)
			return i;
	return -1;
}

class SolveParameters
{
public:
	SolveParameters();
	void Read_SolveParameters(string SolveParametersConfigName);
	void InitializeAfterSetup();
	solveOptions sOpt;
	string configName;
	string configPPName;
	// what is the range of serial numbers -> looping over random fields
	int serialNumber_st, serialNumber_en;
	// version number (related to sfcm and sfcm_gen)
	int versionNumber_st, versionNumber_en;
	// the name of the config maker generator name that creates combinations for the run
	string version_configMakerGenName;
	// if true, the instruction file will be regenerated regardless of if it already exists or not
	bool version_configMaker_forceRewrite;
	// this ofset causes the versionStr used in output file names to be incremented by value below.
	// Use: We won't overwrite previous runs possibly from other instructions
	int versionOffset;
	// loop_version_first_serial_second: see the function void Solve_all_serialNumbers
	bool lv1s2;
	// to only write out necessary information
	int low_disk_space;

	// this is used to generate script for parallel runs (.sh files)
	// if the number < 0, the scripts are not generated
	int numParallelRuns;
	// the name for the main config that sometimes is used to read the content of this
	string solveParametersConfigName;
	// stores if the file above is read
	bool b_solveParametersConfigName;

	// the following are only used for random field tests - for that the configName should be provided without extension
	bool isPeriodic;
	double xm, xM;

	bool isDomain;
	bool PPS2_outside;
	// 0 does not delete, 1 deletes, 2 only deletes in low disk scenario
	int delete_runFolders;
	// 0 inside run folder, 1 outside run folders, 2 only outside when versionNo >= 0
	int vis_outside;
};

extern SolveParameters solvePara;
// timeIndex < 0, does not print time position
void GetSubdomainIndexed_TimeIndexed_FileName(string& fileName, int subdomainNo, int timeIndex, const string& specificName, string ext = "txt", bool addUnderline = false);
// spaceIndex < 0, does not print space position
void GetSubdomainIndexed_SpaceIndexed_FileName(string& fileName, int spaceIndexOr_domainIndexLeft, const string& specificName, string ext = "txt", int domainIndexRight = -1);


class extermumhold
{
	friend 	ostream& operator<<(ostream &output, extermumhold &exter);
	friend 	istream& operator>>(istream &input, extermumhold &exter);

public:
	extermumhold();
	double value;
	double x_loc;
	bool e_initialized;

	bool updateMin(double value_in, double x_loc_in);
	bool updateMax(double value_in, double x_loc_in);

	void MergeMin(const extermumhold& exte_in);
	void MergeMax(const extermumhold& exte_in);
};

class statHolder
{
	friend 	ostream& operator<<(ostream &output, const statHolder &stat);
	friend 	istream& operator>>(istream &input, statHolder &stat);

private:
	extermumhold max;
	extermumhold min;
	long counter;
	double sum;
	double sumSquares;
	double sumMeasure;
	bool useMeasure;

public:
	string name;
	string nameLatex;

	statHolder(bool useMeasureIn = false);
	statHolder(string& nameIn, string& nameLatexIn, bool useMeasureIn = false);
	void setEmpty();
	inline void set_useMeasure(bool useMeasureIn = true) { useMeasure = useMeasureIn; };

	//	void printBrief(ostream& output);	
	void setName(const string& nameIn, const string& nameLatexIn);
	int getCount()const { return counter; }
	void setCount(long counterIn) { counter = counterIn; }
	inline double getSum() const { return sum; };
	double getAverage()const;
	void setAverage(double aveIn);
	double getStandardDeviation() const;
	void setStandardDeviation(double sdiv) { sdiv_saved = sdiv; b_sdiv_saved = true;	};
	double getCOV() const;
	inline void set_rN(double rNIn) { rN = rNIn; }
	inline double get_rN() { return rN; }
	void update(double value_in, double x_loc_in, double weightIn = 1.0);
	void MergeStatHolder(statHolder& stat);
	inline double getMin() const { if (counter > 0) return min.value;  return 1e40; };
	inline double getMax() const { if (counter > 0) return max.value;  return 1e40; };
	inline double getMinLoc() const { if (counter > 0) return min.x_loc; return 1e40; };
	inline double getMaxLoc() const { if (counter > 0) return max.x_loc; return 1e40; };
	inline void setMin(double minVal, double minLoc) { min.updateMin(minVal, minLoc); }
	inline void setMax(double maxVal, double maxLoc) { max.updateMin(maxVal, maxLoc); }
	inline double get_measure() const { if (!useMeasure) return (double)get_counter(); if (counter > 0) return sumMeasure; return 1; }
	double getValue(int setStatOp_type_i);

private:
	void setSquareSumFromStandardDeviation(double sDeviation);
	inline long get_counter() const { if (counter > 0) return counter; return 1; }
	double rN; // normalized count
	double sdiv_saved;
	bool b_sdiv_saved;
};


#endif
