#ifndef SL_INTERFACE_CALCULATOR__H
#define SL_INTERFACE_CALCULATOR__H

#include "globalMacros.h"
#include "globalTypesClasses.h"
#include "SLBulk_Properties.h"
#include "SL_Interface_PtData.h"
#include "SLInterfaceFracturePF.h"
#include "TSR1D.h"

// this file combines things from other files:
// 1. mechanical fields (permanent and temporary) for calculation of star values
//		1.a SL_interfacePPtData				-> permanent
//		1.b SL_interface_Temp_PPtData		-> temporary
// 2. SL_Interface_Fracture_PF:			fracture / contact flags & properties for the interface 
// 3. SL_Elastic_InterfaceProperties:	impedances, etc. needed to calculate target values
// 4. Characteristics: left & right goings
// 5. SLInterfaceCalculator: takes all the above (plus a few below) to calculate target values
// 6. Earlier values that form the following:			ptTimeSequenceSlns
//				OR
//	  directly providing them as
//			interface damage value				tPtSlns->interface_damageLatestValue
//			side vel, traction, displacement	sl_side_temp_ptData[SDL], [SDR]->v_downstream_latestValue	, sigma_downstream_latestValue , u_downstream_latestValue	
//			v_r component for ring problem		v_r_latestValue, ...

// 7. factors that modify interface properties: like fracture strength and displacement factors

// data is encapsulated in SLInterfaceCalculator
// IMPORTANT: INPUTS:	
//	a.						wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ
//	b.						tPtSlns:
//									most-recent traces of side (v, s, u)		interface Damage (see item 6 above)	v_r (ring problem)
//  c.						sigmaCFactor, deltaCFactor (for inhomogeneous media)

// Data assisting in computation
//	a.						interfacePFs: interface properties and flags 
//  b.						ts_bulkProps: elastic properties (impedance, etc) from the two sides that assist in computing Riemann solutions 

// Final result
// pPtSlns (traces of (v, s, u) on the two sides		interface damage )

// this class is only used to update downstream characteristics for problems with source term
class SL_OneInterfaceAllTimes;
class Periodic1IntrFrag;

class SLInterfaceCalculator
{
public:
	SLInterfaceCalculator();
	~SLInterfaceCalculator();
	// return value is the type of adaptivity flag for this point
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//				INITIALIZATION functions. Afterwards the main function below can be called
	//
	//------------------
	// A)	Interface ELastic (ts_bulkProps) and Contact/Fracture (interfacePFs)	
	//					For 1 point calculations they can be read from input file. But in general, they are read outside for fixed special positions and are passed along to this function for calculations based on their flags/parameters.
	void Initialize_setFromOutside_ElsticFractureProperties(SL_Elastic_InterfaceProperties*	ts_bulkPropsIn, SL_Interface_Fracture_PF* interfacePFsIn = NULL);

	//------------------
	// B) permanent Pt solutions
	//		B1) - time step(n +1) for a given fixed position, there will be many of these (e.g. vector indexed by time position).
	//								The next point (i.e. step(n+1) is send to the function
	void Initialize_setFromOutside_Current_Step_EmptyPermanentPt_Storage(SL_interfacePPtData* pPtSlnsIn); // , double xIn = 0);

	//		B2) - Sequence of solutions is stored in SL_interfacePPtData_Time_Seq
	//			For time step 1 - Initial condition (sz 1) is sent to the function. For next steps a a few previous steps are sent to the function

	void Initialize_setFromOutside_Earlier_Steps_NonEmptyPermanentPts_Storage(SL_interfacePPtData_Time_Seq& ptTimeSequenceSlnsIn);

	//------------------
	// C) - traces of characeristics OR stress velocity traces from the two sides. There are two ways of doing this
	//		- in both approaches - traces of displacement are also provided, in case fracture is concerned and no earlier solution is given (for example in step B2)
	//		- inputs are sent as pointers to check if something is NULL, it means it's not given	
	//		- prerequisite: step A an B (specially B1) are called earlier
	//		- ALL values (characteristics OR stresses) are given WITHOUT in-situ part

	//		For many point solutions ONLY incoming characeristics (OR stresses) need to be provided and traces of current displacements and velocity will be obtained from earlier solutions
	//		For one-point computations traces of current displacement and velocity from the two sides should be given

	//		C1)	providing left and right characteristics (and possibly displacement and/or velocities)
	//
	//		Suitable:					METHOD of CHARACTERISTICS
	//	
	//					(for BCs, only the characteristics from inside the domain needs to be meaningful)

	void Initialize_setFromOutside_Incoming_Characteristics_etc(const VEC* wlSide_rGoing_WO_in_situPtrIn, const VEC* wrSide_lGoing_WO_in_situPtrIn, const VEC* traceCurrentDisplacementLeftIn = NULL, const VEC* traceCurrentDisplacementRightIn = NULL, const VEC* traceCurrentVelLeftIn = NULL, const VEC* traceCurrentVelRightIn = NULL);

	//		C2)	providing left and right stresses and velocities (and possibly displacements)
	//
	//		Suitable:					Time Marching DG and spacetime DG methods; FV, etc. where s, v, u and traces from the two side of the interface
	//		
	void Initialize_setFromOutside_Traces_vel_stress_etc(const VEC* traceCurrentStress_WO_in_situ_LeftIn, const VEC* traceCurrentStress_WO_in_situ_RightIn, const VEC* traceCurrentVelLeftIn, const VEC* traceCurrentVelRightIn, VEC* traceCurrentDisplacementLeftIn = NULL, VEC* traceCurrentDisplacementRightIn = NULL);

	//------------------
	//	D) If applicable -> set different from 1 sigmaCFactor and deltaCFactor (inhomogeneous fracture properties)
	void Initialize_setFromOutside_Fracture_InhomogeneousFactors(double sigmaCFactorIn = 1.0, double deltaCFactorIn = 1.0);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//				MAIN computation function to be called from outside
	//
	AdaptivityS Main_Compute_OnePoint(bool& accept_point, IOF_type iof_finalVals, IOF_type iof_scalarVals, double space_or_time, ostream* outScalars = NULL, ostream* outFinalSln = NULL, ostream* outAdaptivity = NULL, ostream* outIterationConv = NULL);

	void Read(istream& in);
	void Output_SLInterfaceCalculator(bool accept_point, IOF_type iof_finalVals, IOF_type iof_scalarVals, double space_or_time, ostream* outScalars, ostream* outFinalSln, ostream* outAdaptivity, ostream* outIterationConv);

	VEC wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ;

	// incident wave
	InterfaceLocation1DT incidentSide; // ilt_noSided: does not have incident, ilt_left, ilt_right wave enters from the left/right sides
	VEC incientIncomingStress;

	// characteristics from left (rGoing) and right (lGoing), respectively
	// mechanical fields, permanent (v, s, u)
	// SOURCE: Supposed to come from outside -> from the next points permanent data -> fixed space position -> there are many of these
	SL_interfacePPtData* pPtSlns;
	// temporary data that will be either written to file or is just for temporary computation. It will not be stored in time sequence of solutions for different space values
	// SOURCE: UNLIKE other pointers in this class, it's supposed to be created by new in this object and deleted upon completion of calculations and output of results.
	SL_interface_Temp_PPtData* tPtSlns;

	// ELastic (ts_bulkProps) and Contact/Fracture (interfacePFs) 
	// SOURCE: Supposed to come from outside ->											fixed space position -> only one of these for all points
	// two sides bulk properties
	SL_Elastic_InterfaceProperties*	ts_bulkProps;
	// properties and flags at the interface
	SL_Interface_Fracture_PF* interfacePFs;

	// only needed for problems with source term
	SL_OneInterfaceAllTimes * oneIntAlltimesPtr;

	// factors for strength and length scales
	double sigmaCFactor, deltaCFactor;
	// this must be initialized in iteration zero from previous time step
	double dcont;
//	double x;
	// these are earlier solutions used to start the solution for this point for the first iteration or later ones. It uses alpha's and other members in SLFractureGlobal_Configuration
	// they are indexed in revered order (0 -> step n, 1 -> step n - 1, ...)
	SL_interfacePPtData_Time_Seq* ptTimeSequenceSlns;

	//////////////////////////////////////////////////////////////////////////////////////////
	////////////////  Non-input data that are computed as a part of the solution of this point
	/// this is computed from w's WO_in_situ
	VEC wlSide_rGoing, wrSide_lGoing;
	// it's computed from difference of current and previous times
	double current_delT;
	int current_timeIndex;
	// error checks 
	// iteration -> iteration
	Pt_Error_Data within_step_convergence_check;

	// between steps -> if the current point should be accepted
	Pt_Error_Data between_steps_timestep_check;
	AdaptivityS	  between_steps_adaptivity_s;
	TSR1d_stepn_to_np1* ortiz1DTSRptr;

	void Set_tPtSlns(SL_interface_Temp_PPtData* tPtSlnsIn, bool tPtSlnsDeletableIn);
	void Set_pPtSlns(SL_interfacePPtData* pPtSlnsIn, bool pPtSlnsDeletableIn);
private:
	////////////////////////////////////   level 1 functions called in Main_compute_OnePoint ...
	// except wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ other data needed for star calculation (e.g. traces of v, stress, ...) is either directly given or obtained 
	// from ptTimeSequenceSlns (and earlier solutions for this point itern > 0)
	void InitializeOrUpdate_Most_Current_Values();
	// return value:
	//			true if no more iterations is needed. This happens if
	//					a. convergence is achieved.
	//					b. convergence is not achieved BUT at maximum number of iterations refine_after_max_num_iter.
	//							In this case the acceptance of this point is stored in between_steps_adaptivity_s.a_flag (above) based on the following:
	//								g_slf_conf->refine_after_max_num_iter == false
	//								b1. a_none:		
	//								g_slf_conf->refine_after_max_num_iter == true
	//								b2. a_refine	and proposed time (current_delT * g_slf_conf->refine_after_max_num_iter_del_t_factor < g_slf_conf->min_del_t)
	//								b3. a_terminate_run_prematurely					opposite of condition above
	bool Check_withinstep_iteration_convergence();
	// this function is called after Check_withinstep_iteration_convergence returns true (iterations for time step n + 1, have ended)
	// returns value is adaptivity flag request
	void Copy_temp_to_final_vals_Check_adaptivityFlag(bool& accept_point);

	//////////////////////////////// functions needed to compute star values

	// updates characteristics, compute star values, and put things in permanent storage pPtSlns
	void Compute_1Cycle_vsStar_Calculation();
	void Compute_1ConvergedImplcit_vsStar_TSR_OrtizSimple_Calculation();

	void Compute_Star_sv_from_RiemannItoIII_toStar();

	// I, II, III solutions
	void Compute_Riemann_I_Bonded();
	// return value: relCont 
	SLFF_ContactType Compute_Riemann_III_Separation(SLFF_InterfacialDamageModeType damageOffOnMix, double& relCont);
	onOffT Compute_Riemann_II_Slip(double& relStick, SLFF_SlipType& slipOffOnMix);

	void SetStressComputeVelocityStarFromStressStar(const VEC& sigmaStar, RiemannMode_StorageT rMode);

	void Create_TemporaryStorage();

	// iteration number - the point can go through multiple iterations until convergence is achieved
	int itern;

	inline bool DirectlyCalculateBondedInFinalSln() const;
	bool b_directlyCalculateBondedInFinalSln;

	// handling pointers
	bool tPtSlns_Deletable;
	bool pPtSlns_Deletable;
	bool ts_bulkProps_Deletable;
	bool interfacePFs_Deletable;
	bool ptTimeSequenceSlns_Deletable;

	bool terminateConditionMet;
	SLInterfaceCalculator(const SLInterfaceCalculator& other);
	SLInterfaceCalculator& operator=(const SLInterfaceCalculator& other);

	// Incident functions:
	void UpdateCharacteristics_FromIncident();
	void UpdateStarValues_FromIncident();
	// periodic problem, ring type
	void Update_Velocity_Left_Ring_Opened1D();
};

// intentionally the function is moved out of SLInterfaceCalculator so it can more easily be called from outside
// Inputs: lines 1 & 2
// Line 1: mechanical fields sigmaI, traces of u and v from either side and current damage value; maxEffDelU is the maximum effective deltaU so far
// Line 2: sigmaC and deltaCfactors are possible inhomogeneity multiplier of deltaC and sigmaC in interfacePFs, interfacePFs contains all important interface properties; YAve_n, YAve_t are inverse impedance parameters used in RBH damage rate model
// Line 3 are the outputs:
//		sigmaI, delu, delv complexes (normal, shear, etc. needed along the way in computations)
//		DDot, DTarget, DSrc_Allix, DSrc_RBH: parameters related to damage model (if active)
//		sigmaIII, actual target traction for separation

/////	separation functions
void Compute_Riemann_III_Separation_Damage_Rate_Aux(SLFF_InterfacialDamageModeType damageOffOnMix, SLFF_ContactType contactOffOnMix, Vec_nt& sigma_I_nt_parts, Vec_nt& del_u_nt_parts, const VEC& vL, const VEC& vR, double maxEffDelU, TSR_loadingStages tsrStagePrevStepValue, TSR_loadingStages& tsrStageLatestValue, double damage,
	double sigmaC, double deltaC, double tau, const SL_Interface_Fracture_PF* interfacePFs, double dn, double dt,
	Vec_nt& del_v_nt_parts, double& effS, double& DDot, double& DTarget, double& DSrc_Allix, double & DSrc_RBH, VEC& sigmaIII, bool& terminateConditionMet);
void ComputeTSR_Tractions(SLFF_TSRType tsrModel, Vec_nt& del_u_nt_parts, Vec_nt& sigma_I_nt_parts, double beta_delU, double deltaC, double sigmaC, double kCompressive, VEC& S, double max_delU, TSR_loadingStages tsrStagePrevStepValue, TSR_loadingStages& tsrStageLatestValue, bool& terminateConditionMet);
double Compute_Interface_EffectiveStress(const SL_Interface_Fracture_PF* interfacePFs, double sigmaC, Vec_nt& sigma_I_nt_parts, double &strengthNormalizer);



///// slip functions
onOffT Compute_Riemann_II_Slip_Aux(SLFF_SlipType slipTypeIn, Vec_nt& sigma_I_nt_parts, Vec_nt& del_u_nt_parts,
	const SL_Interface_Fracture_PF* interfacePFs, const SL_Elastic_InterfaceProperties*	ts_bulkProps,
	VEC& sigmaII, double& k, double& thetaTauFriction, double& delvIIAbs, double& thetaDelv, double& powerIDissipationSlip);

class Slip3D_1Dir
{
	friend ostream& operator<<(ostream& out, const Slip3D_1Dir& dat);
public:
	double n[2];
	double theta;
	double betaTheta;
	double alphaThetaTheta;
	double damaged_sn;
	double damaged_c;
	double k;
	double delvII;
	double sIInVal;
	double tauIhetaVal;
	double tauIIhetaVal;
	double power;
	bool validSlip;
};

class Slip3D
{
	friend ostream& operator<<(ostream& out, const Slip3D& dat);
public:
	Slip3D();
	double tauItVal[2];
	double sInval;
	int numPhiSegments;
	double delPhi;

	vector< Slip3D_1Dir> slp3ds;
	int indexMaxPower;
	double theta4MaxPower;
	double maximumPower;
};

class Periodic1IntrFrag_TimeStageStorage
{
public:
	Periodic1IntrFrag_TimeStageStorage();
	void CopyData(const Periodic1IntrFrag_TimeStageStorage& other);
	void Initialize(PITSS pitssIn, const string& nameBase, Periodic1IntrFrag* perConfIn);
	void Print();

	int b_set;
	PITSS pitss;
	Periodic1IntrFrag* perConf;
	unsigned int timeIndex;
	vector<double> vals; // U, K, delU, ...
	vector<double> spatial_stress_vals, spatial_vel_vals;
	fstream out;
};

class Periodic1IntrFrag
{
public:
	Periodic1IntrFrag();
	~Periodic1IntrFrag();
	void Output_Periodic1IntrFrag_Header(ostream& out);
	void Output_Periodic1IntrFrag(ostream& out);
	double relTol; // to make small values
	// inputs:
	// all quantities are noramlized to result in 
	int aIndex;
	int lIndex;
	double a; // normalized loading rate
	double l; // normaliezd fragment size used
	bool isActive;
	bool isExtrinsic;
	double energyScale; // factor of deltaC * sigmaC that gives sigmaC
	string nameBase;
	string nameOne_a_SharedAll_l;
	string nameOne_a_One_l;
	// bulk interface PP
	int numSpatialSubsegments_BulkInterfacePoints_Print_4PP;
	bool useRepeatedSimpsonRuleForHigherOrders_4PP;
	int step4_segment_vsigma_output;

	void Initialize_Periodic1IntrFrag();
	// returns PITSS_none if the run should continue, otherwise the run terminates
	PITSS UpdateStats(double timeNew, double delu, double delv, double sigma);
	void Finalize_Periodic1IntrFrag();

	///// related to Periodic1IntrFrag_TimeStageStorage
	void Updata_Periodic1IntrFrag_TimeStageStorage(PITSS pitss);

	// computed
	//Zhou_2006_Effects of material properties on the fragmentation of brittle materials.pdf
	// my time scale is twice larger so equations are different

	double t0; // = 1/a
	double t_SigmaMax_dilute;
	double t_SigmaZero_minus_SigmaMax_dilute;
	double t_SigmaMax_real;
	double t_SigmaZero_minus_SigmaMax_real;
	double t_dilute_approx, t_dilute;
	double t_dilute_approx_p, t_dilute_p; // _p means values are normalized by Zhu06 and Drugan's paper 
	double l_dilute_approx_p, l_dilute_p;
	TSR_XuNeedleman tsr_xn;

	double delu4Sigma0;
	double l_dilute_approx; // equations (17 & 19)
	double l_dilute; // equation (16) = time_dilute (length and time needed for fragments to not interact with each other)
		
	// Drugan reference:
	// Xu-Needleman	: REF = MAX : correponds to max stress in TSR, l_Drugan is the length for which solid speed is zero
	// Ortiz		: REF = FIN : correponds to final stress in TSR (delta = 1, sigma = 0) as for initial stage vSolid is always going to be > 0
	double l_SigmaRef_Drugan_Dilute;
	bool diluteFractureModel; // l >= l_dilute
	double a_p, l_p; // zhu6_a, zhu6_l; // epsilonDotNonDimensional; // (11)
	double zhu6_epsilonDotScale; // (13)
	double zhu6_sBarScale; // (21) fragment size scale
	double l_Zhu6a; // (28)
	double l_Zhu6b; // (29a)
	double l_Grady;
	double l_Glenn;

	double la;
	double sigma_t0;
	double aInv;
	double log10a;
	double log10l;
	double suggeste_delt;
	double suggested_final_time4ZeroSigmaC;
// see the formula how predicted sigma max is computed

	double tFailure;

	double energy_diss_per_length;
	double energy_diss_per_length_At_t_dilute;

	bool energy_diss_per_length_At_t_dilute_set;
	// this is another way to normalize energy dissipation of different models. 
	// It's normalized by input energy throught stress source term, from time 0 to t
	// int 0 to t of I = a * spatialMean(stress) dtau -- t is the final time
	double inputEnergyFromMeanStress_PL;
	// max value that I could have taken above, if sigma = at -> IMax = a^2t^2 E / 2 -- t is the final time
	double inputEnergyFromMeanStress_MaxPossVal_PL;

	double energy_diss_real_inp_E_t_final; // energy_diss_per_length / inputEnergyFromMeanStress_PL
	double energy_diss_max_inp_E_t_final;  // energy_diss_per_length / inputEnergyFromMeanStress_MaxPossVal_PL


	// bulk interface PP
	vector<double> xs_wrCenter;
	vector<double> weights;
	unsigned int sz_xs;
	double half_l;

	///
	unsigned int timeIndexNew;
	vector<double> currentStepVals;
	double max_bsigma, max_bsigma_pw;
	vector<statHolder> stat4Vals;
	fstream outAllVals, outStress, outVel;
	vector<string> velHeader, stressHeader;
	vector<double> velVals, stressVals;
	bool delvsZeroObserved;

	// the time to check if we have reacheed dilute limit
	double t_dilute_zeroStressCheck;

	///// related to Periodic1IntrFrag_TimeStageStorage
	Periodic1IntrFrag_TimeStageStorage stageSlns[PITSS_SIZE];
	PITSS terminateState;

private:
	void setEmpty_Periodic1IntrFrag();
	void set_energy_diss_per_length_4_t_dilute(double time, double deltau);
};

extern Periodic1IntrFrag per_if;

bool SetSlip3DParameters(const SL_Interface_Fracture_PF* interfacePFs, const SL_Elastic_InterfaceProperties*	ts_bulkProps, Slip3D& s3d, Vec_nt& del_u_nt_parts, Slip3D_1Dir& s3d1Dir);

void Compute_Directional_FrictionParameters(Vec_nt& del_u_nt_parts, double delV_Angle, const SL_Interface_Fracture_PF* interfacePFs, double& k, double& damaged_sn, double& damaged_c);
void GetDirectionalFrictionCoefficientFromRaw_k(const SL_Interface_Fracture_PF* interfacePFs, double k_upper, double k_lower, double directional_tangential_delU, double& k);
void Compute_Directional_k(double delV_Angle, const SL_Interface_Fracture_PF* interfacePFs, double& k_upper, double& k_lower);
void C1_interpolateC(const double p0, const double m0, const double x0, const double p1, const double m1, const double x1, const double x, double & value, double & slope);


#endif

