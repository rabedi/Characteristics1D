#ifndef MIN_CROSSING_SOLVER__H
#define MIN_CROSSING_SOLVER__H

#include "Dims.h"

typedef enum {sgnT_none = -2, sgn_m = -1, sgn_0, sgn_p, sgnT_SIZE} sgnT;
sgnT Sgn(double val, double tol);

// successful convergence:	ct1d_yes 
// else
//		ct1d_xmin_lim: reach xmin limit and cannot test smaller values
//		ct1d_xmax_lim: reach xmax limit and cannot test larger values
//		ct1d_notCrossing: closest point that gets to zero for finding crossing BUT the function always stays on one side
//						in this case if this holds true for the end points, the function returns cr1d_xmin_lim, ct1d_xmax_lim
//		ct1d_failedSln: an individual x -> y solution failed
//		ct1d_gettingSameXcrd: points are not distinguished well in x axis and the iteration keeps getting the same x
//		maximum number of iteration hit:
//				ct1d_large_delx:			delx too large
//				ct1d_large_dely:			dely too large
typedef enum {ct1d_yes, ct1d_notCrossing, ct1d_xmin_lim, ct1d_xmax_lim, ct1d_failedSln, ct1d_nanSln, ct1d_gettingSameXcrd, ct1d_large_delx, ct1d_large_dely, convergeT1D_SIZE} convergeT1D;
// -2, failed solution, -1 nan, 0 numerical solution but not converged, 1 numerical solution and converged
int Crossing_Sln_Mode(convergeT1D dat);
int Extremum_Sln_Mode(convergeT1D dat);
int Sln_Mode(convergeT1D dat, bool isExtremum);

string getName(convergeT1D dat);
bool name2Type(string& name, convergeT1D& typeVal);
ostream& operator<<(ostream& out, convergeT1D dat);
istream& operator>>(istream& in, convergeT1D& dat);

class gFx2y;
class Solver1D;

class ConvergenceLog
{
	friend ostream& operator<<(ostream& out, const ConvergenceLog& dat);
public:
	ConvergenceLog();
	// int isExtremumSln = -1, valid solution for crossing and extremum is printed,
	// 0 or 1 crossing or extremum os sought and the correct solution is printed
	void Print(ostream& out, bool printHeader, bool printData, bool printLong, bool expandConvType, int isExtremumSln = -1) const;
	convergeT1D convType;
	unsigned int iterCnt;
	// in extremum search, the extremum may stay the same point for many iterations, this keeps track on number of times this has not changed
	unsigned int iterExtremumNotChanged;
	vector<double> delxs, delys;
};

class genIndexVal
{
	friend istream& operator>>(istream& in, genIndexVal& dat);
	friend ostream& operator<<(ostream& out, const genIndexVal& dat);
public:
	genIndexVal(int index_main_in = 0, double x_in = 0.0);
	inline bool IsEmpty() const { return (index_sec < 0); };
	inline unsigned int getSize_y() const { return ys.size(); };
	bool operator<(const genIndexVal &other) const;
	bool operator==(const genIndexVal &other) const;
	void Write_genIndexVal(ostream& out) const;
	void Write_genIndexVal_JustVals(ostream& out) const;
	bool Read_genIndexVal(istream& in);
	void MakeNaN();

	int index_main;
	int index_sec;
	double x;
	vector<double> ys;
};
class genIndexVal_DBH
{
public:
	genIndexVal_DBH();
	bool found;
	unsigned int pos;
	unsigned int y_size;
//	genIndexVal* givPtr;
};

// stores a history of genIndexVal
class genIndexVal_DB
{
public:
	void Sort(bool equality_by_x = false, double tolx = -1, vector<genIndexVal>* dbIn = NULL);
	unsigned int size() const { return db.size(); };
	void CopySort(vector<genIndexVal>& dbOut, bool equality_by_x = false, double tolx = -1);

	// print_dbAsIs
	//		== true, prints db as is
	//		== fasle, creates a copy and prints based on the next two bools
	void Write_genIndexVal_DB(ostream& out, bool print_dbAsIs = true, bool equality_by_x = false, double tolx = -1);
	void Write_genIndexVal_DB_WithHeader(gFx2y* functionIn, ostream& out, bool print_dbAsIs = true, bool equality_by_x = false, double tolx = -1);
	// read_dbAsIs is similar to above but just for the reading
	void Read_genIndexVal_DB(istream& in, bool read_dbAsIs = true, bool equality_by_x = false, double tolx = -1);
	// adding point
	// checkIfPtExists checks if the point already exists
	//		if true: equality_by_x and tolx are used to check for existence
	// addPt may change if checkIfPtExists == true and the addPt is changed

	// return value (index of the added point, and size of it's y)
	genIndexVal_DBH Add_Pt(genIndexVal& addPt, bool checkIfPtExists = true, bool equality_by_x = false, double tolx = -1);

	// if true operator< sorts it by x, otherwise by main and sec indices
	static bool g_genIndexVal_DB_equality_by_x;
	static double g_genIndexVal_DB_tol_x;

	vector<genIndexVal> db;
};

// used to find best point closest to an x0 for crossing or optimal value

class cross_optimal_helper
{
public:
	unsigned int ii, i, inext;
	double x, xnext;
	sgnT signi, signinext;
};

// a general function that goes from x to y
// the parent class uses the test function
// y0 = (x - para[0])^2 - 9 for the range [para[0] - 2, para[0] + 4]. It has a minimum at -2 and zero at 5
// y1 = sin(2.0 x)
class gFx2y
{
public:
	// num_y is the number of y's the class returns
	// tol_ys need not to be set, if stays empty, they are set outside (e.g. from Solver1D_1posConf)
	virtual void InitializeValues(Solver1D* conf, const map<string, string>& str_mapIn, const vector<unsigned int>& indices4ParasIn, const vector<double>& parasIn, double& xMin, double& xMax, double& tol_x, vector<double>& primary_xs, vector<double>& secondary_xs, int& num_y, vector<double>& tol_ys);
	// x and indices of giv are set, y is sought
	// returns true if the solution is successful
	virtual bool ComputeValue(genIndexVal& giv);
	// calls the function below
	void Print_Header(ostream& out, int num_y = -1) const;
	virtual void Print_YHeader(ostream& out, int num_y = -1) const;
	vector<unsigned int> indices4Paras;
	vector<double> paras;
	// this map helps to read any other kind of data that gFx2y needs
	map<string, string> str_map;
};

class Solver1D_1posConf
{
	friend istream& operator>>(istream& in, Solver1D_1posConf& dat);
public:
	Solver1D_1posConf();
	unsigned int Initialize_Solver1D_1posConf();
	bool Read_Solver1D_1posConf(istream& in);

	string nameBase;
	unsigned int pos_y;
	double tol_y;
	bool isActive;
	bool isExtremum; // if true, extremum values are solved, else crossing values
	vector<bool> isMins;
	vector<double> crossing_ys;
	vector<double> x0s;

	///// set-up after initialization
	string namep0;
	unsigned int sz;
	vector <bool> vec_isExtremum;
	vector<string> vec_name;
	vector<double> tol_ys;
};

class Solver1D
{
public:
	Solver1D();
	~Solver1D();

	void Read_Solver1D(istream& in);

	// if configName == "none", posConfs2Solve should have been provided indirectly 
	void MAIN_ProcessConfigFile(gFx2y* functionIn, string confName = "none");

	// addAdditionalPoints: beyond initial points computed (e.g. mainIndex = 0, 1 provided by the function) no additional points are added and the best point in the data is chosen
	ConvergenceLog SolveYCrossingOrExtremum(genIndexVal& sln, gFx2y* functionIn, unsigned int y_pos, bool isExtremum, bool isMin, double crossing_y, bool addAdditionalPoints, double x0, double tol_y, bool doInitialization = false);

	ConvergenceLog SolveYCrossing(gFx2y* functionIn, genIndexVal& sln, double crossing_y, unsigned int y_pos, bool addAdditionalPoints, double x0, double tol_y, bool doInitialization = true);
	ConvergenceLog SolveYExtremum(gFx2y* functionIn, genIndexVal& sln, bool isMin, unsigned int y_pos, bool addAdditionalPoints, double tol_y, bool doInitialization = true);
	unsigned int SolveYCrossings(gFx2y* functionIn, vector<genIndexVal>& slns, vector<ConvergenceLog>& cls, double crossing_y, unsigned int y_pos, bool addAdditionalPoints, double tol_y, bool doInitialization);

	vector<Solver1D_1posConf> posConfs2Solve;
	// for each combination in "posConfs2Solve", there are 3 possible solve sequences:
	// possible cases are stored in sz_solve_mode
	// 1. Solve after main and secondary points are added. No solution for extremum or crossing points is done at this stage
	bool do_posConfs_first_notAddPtSolve;
	// 2. Actual solution for adding new points for crossing or extremum are done
	bool do_posConfs_AddPtSolve;
	// 3. After all solutions are done, yet, another time no new point is added and extremum / crossings are obtained
	bool do_posConfs_second_notAddPtSolve;

	// print primary points of different options into one file
	bool b_print_PrimaryPoints;

	// for mode 1 solution where extremum and crossings are calculated in a restart of the general solve, 
	// often resolution of the underlying grid (delx) is not changing but a new solution is sought.
	// If resolution is not going to change (often the case) --> use			b_recalculate_solutions = false
	// else b_recalculate_solutions = true
	bool b_recalculate_solutions;

	vector<double> paras1D; // like loading rate for the fragmentation problem
	vector<unsigned int> indices4paras1D;
	map<string, string> str_map;
	string baseName;
	// in extremum search is extremum stays the same for "maxNumIterExtremumNotChanging" number of iterations, the point is reported as extremmum
	unsigned int maxNumIter, maxNumIterExtremumNotChanging;

	// read results means it reads results if already have been computed
	//	read_results: if true if reopens the log file and starts from scratch
	// write results means that it writes results of x, yVec calculations in log files
	// equality_by_x: means that too check if a solution exists, the x coordinate of previous solutions is checked rather than iMain, iSecondary indices
	bool read_Results, write_Results, equality_by_x;

	void InitializeFunction(gFx2y* functionIn);
	// functions called after InitializeFunction
	// solve crossing for the function
	// these two functions search for the crossing and extremum and potentially add to the point set
	ConvergenceLog Solve_x4_Crossing(genIndexVal& sln, double crossing_y, unsigned int y_pos, double x0, double tol_y);
	ConvergenceLog Solve_x4_Extremum(genIndexVal& sln, bool isMin, unsigned int y_pos, double tol_y);

	// similar to above, but no point is added
	ConvergenceLog Solve_x4_Crossing_NoPointAdded(genIndexVal& sln, double crossing_y, unsigned int y_pos, double x0, double tol_y);
	ConvergenceLog Solve_x4_Extremum_NoPointAdded(genIndexVal& sln, bool isMin, unsigned int y_pos);

private:
	void Store_pts(bool print_dbAsIs = true);
	void Store_pts_unsorted_SlnEnd(int y_pos = -1);

	void Restore_pts(bool read_dbAsIs = true, bool compute_ys_not_computed = true);
	genIndexVal_DBH Compute_Add_pt_value(genIndexVal& giv, bool& pointWasComputed, bool checkIfPtExists);
	// x0: is the location we want to find the solution closest to. Sometimes there are multiple crossings and want to choose the solution closes to x0
	// if x0 > xM (below) x0 = xM
	// if x0 < xm (bel0w) x0 = xm
	bool Solve_x4_Crossing_Aux(unsigned int cntrAddedOnSides, double crossing_y, unsigned int y_pos, double x0, double tol_y_zero, int& jCrossingm1, int& jCrossing, sgnT& sgnjCrossing, ConvergenceLog& cl);
	bool Solve_x4_Extremum_Aux(unsigned int cntrAddedOnSides, bool isMin, unsigned int y_pos, int& jExtremum, int& jExtremumm1, int& jExtremump1, ConvergenceLog& cl);

	unsigned int getValidValues(unsigned int y_pos, vector<int>& poss);

	string name_log, name_WOExt, name_log_x_unsorted;
	gFx2y* function;
	double xMin, xMax, del_secondary_x, tol_x, tol_x_eq_check;
	vector<double> tol_ys;
	vector<double> primary_xs, secondary_xs;
	int num_y;
	fstream out_error;

	genIndexVal_DB pts;
	int startInds[4];
	void SolveOneSolveMode_AfterInitialization(vector<genIndexVal>& slns, vector<ConvergenceLog>& cls, gFx2y* functionIn, bool addAdditionalPoints, int mode_counter);
	void SolveAllSolveMode_AfterReadingConfig(gFx2y* functionIn);

	void Initialize_AfterReading();
	unsigned int sz_indices4paras1D;
	unsigned int sz_posConfs2Solve;
	unsigned int sz_solve_mode;
	vector<bool> solve_mode_addAdditionalPoints;
	vector<int> solve_mode_counter;

	void PrintSln_Header(ostream& out);
	void PrintSln_Values(ostream& out, const genIndexVal& sln, ConvergenceLog& cl, bool isExtremum);
	unsigned int Cross_Helper(double crossing_y, unsigned int y_pos, double tol_y_zero,
		vector<int>& poss, vector<sgnT>& signs, vector<cross_optimal_helper>& helpers);
 
	void Print_PrimaryPoints(genIndexVal& giv, unsigned int secondaryPos);
};


void MAIN_Solver1D(string confName = "default", gFx2y* functionIn = NULL);
void Test_SolveYCrossing();
void Test_SolveYExtremum();

#endif