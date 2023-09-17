#ifndef TSR_1D__H
#define TSR_1D__H

#include "Dims.h"

// delu_n+1 = del_n + delt * (alpha delv_n+1 + (1-alpha) delv_n)
// alpha = 1: backward euler			-> consistent but lower order	
// alpha = 0.5: central					-> results in some inconsistencies
// alpha = 0: forward euler				-> invalid here	

#define TSR_BEuler	1	// chooses alpha = 1 (backward euler) - otherwise  central is used
#define g_tol_TSR		1e-8

#define THROW_TSR(msg){ Write_MainInputs_NonConvergentCase();	THROW(msg); }


typedef enum {
	TSR_loadingStages_none,
	tsrl_dam0_sigI_neg, tsrl_dam0_sigI_less_sigC, tsrl_dam0_sigI_over_sigC, // undamaged stage (delu) = 0
	tsrl_damMid_loadingBranch, tsrl_damMid_unloadingBranchPosDelv, tsrl_damMid_unloadingBranchNegDelv, tsrl_damMid_delun_neg_sigI_neg, tsrl_damMid_delun_neg_sigI_pos, // mid-damaged
	tsr1_dam1_delun_pos, tsr1_dam1_delun_neg, TSR_loadingStages_SIZE
} TSR_loadingStages; // fully damaged delu > delC
// intentionally, want to reduce storage size for 1D problems (slip is moved to the back)
string getName(TSR_loadingStages dat);
void name2Type(string& name, TSR_loadingStages& typeVal);
ostream& operator<<(ostream& out, TSR_loadingStages dat);
istream& operator>>(istream& in, TSR_loadingStages& dat);

// Ortiz branch
// A:			delU = 0, sigma < sigmaC
// B:			softening: sigma = sigmaC (1 - delta/deltaC)
// C:			uloading: sigma = delta/deltaM * sigmaM
// D:			fullDamage, tensile: sigma = 0, delta > 0
// E:			full contact, sigma < 0.0: deltaU = 0
typedef enum { TSR_OrtizStages_slnFound = -2, TSR_OrtizStages_invalid, tsro_A, tsro_B, tsro_C, tsro_D, tsro_E, TSR_OrtizStages_SIZE} TSR_OrtizStages;


class TSR1d_fields_1side
{
public:
	TSR1d_fields_1side();
	// s: stress, u: displacement, v: velocity
	double u, v, s;
};

class TSR1d_2Sides
{
public:
	TSR1d_2Sides();
	TSR1d_fields_1side sides[NUM_SIDES];
	double delv, delu;
};

// inputs: 
//		step_n_tsrStage ->	
//								u, v (not sigma) for both sides
//								step_n_max_delu
//		Zs					impednces
//		ws					incoming characteristics from the two sides
//		delT				time step
//		sigmaC, deltaC		cohesive strength and displacement scale

class TSR1d_stepn_to_np1
{
public:
	TSR1d_stepn_to_np1();
	// true if successful read
	bool Read_MainInputs(istream& in);
	void Write_MainInputs(ostream& out);
	void Write_MainInputs_NonConvergentCase();

	// return sI
	double Compute_step_np1_Solution(double& sStar, double& uLStar, double& vLStar, double& uRStar, double& vRStar, double& maxDelu, TSR_loadingStages& solutionStage, TSR_OrtizStages& stageShort);

	double delT;

	/// material properties
	double sigmaC, deltaC;
	double ws[NUM_SIDES];
	double Zs[NUM_SIDES];

	/// step n
	TSR_loadingStages step_n_tsrStage;
	TSR1d_2Sides step_n;
	double step_n_max_delu;

	//// computed
	/// step n + 1
	double tolsc, toldc;
	double sI; // Riemann 1
	TSR1d_2Sides step_np1[TSR_OrtizStages_SIZE]; // various solutions are possible 
	int calcType[TSR_OrtizStages_SIZE]; // -1: not computed, 0: computed and value is valid, computed and value is invalid
	TSR_OrtizStages validShortStage;
	TSR_loadingStages step_np1_tsrStage;
	double step_np1_max_delu;

private:
	// this boolean is used to gradually increase the tolerances
	int iterN, max_iterN;
	bool maxDelUEqZero;
	string message;
	double tol_TSR;
	double Zeff;
	double delv_np1, delu_np1;
	double step_n_max_sigma;
	// calculate delu, v and Zeff
	void Initialize_TSR1d_stepn_to_np1(double tol_TSRFactor);
	// this is called after Initialize_TSR1d_stepn_to_np1
	bool Compute_Stage_Aux();


	// this is testing if zeroDelU is a valid solution
	// return:		TSR_OrtizStages_slnFound (the input choice stageShort) is correct
	//				else	the other choice to test
	TSR_OrtizStages Test_Stages_DelU_eqZero_A_or_E(TSR_OrtizStages stageShort);
	TSR_OrtizStages Test_Stages_DelU_eqZero_A_or_E_Aux(double& sStar);

	// loading branch: returns TSR_OrtizStages_slnFound if this is the correct case, if not, it sends the case that should be solved
	TSR_OrtizStages Test_Stage_B();
	// unloading branch
	TSR_OrtizStages Test_Stage_C();
	// full damage branch
	TSR_OrtizStages Test_Stage_D();

	// delu is not used and is only used for double-checking calculation.
	// validShortStage shoule be set
	void Set_StarValues_From_sStar(double sStar, double delu);
	void ComputeStarValues_From_delu(double delu, double& sStar, double& delv);
	TSR_OrtizStages Test_newcase(TSR_OrtizStages stageShort2Test);
	void Zero_TSR1d_stepn_to_np1();
};

class TSR_XuNeedleman
{
public:
	TSR_XuNeedleman();
	double a; 
	double Z;
	double E;
	double deltaC;
	double sigmaC;
	double delTFactor; // as a factor of tauC
	double sigmaCFactor4Zero; // sigma never becomes zero, but a fraction like 0.01 or 0.001 is taken as value considered zero past peak
	// loading rate: problem sigma(delta) + Z/2 * deltaDot = at
	// sigma(delta) = sigmaC * ndel exp(-ndel), ndel = delta/deltaC
	void Compute(ostream* outPtr = NULL);
	void OutMain(ostream& out);

	// computed
	double tauC;
	double delT;
	double tSigmaMax;
	double delnSigmaZero;
	double tSigmaZero;
};

// return number of cases tested
int Test_TSR_Ortiz_1D(string nameWOExtIn = "TestFiles/TestOrtiz1D_debug");

#endif