#ifndef SL_INTERFACE_PT_DATA__H
#define SL_INTERFACE_PT_DATA__H

// the classes in this file are for (permanent) and temporary storage of information for a point at a (potential contact/fracture) solid interface
// the permanant ones must be stored at discrete times for a given interface: They are used for calculating characteristics for later times at other interfaces
// the temporary storage is for calculations for a given point and will not be stored

#include "globalMacros.h"
#include "globalTypesClasses.h"
#include "LAfuncsFinalStep.h"
#include "SLFracTypes.h"
#include "TSR1D.h"

// last 4 entries are dissipated powers
// dp: dissipated power
// n -> normal, s -> shear
// s1i_dpn_mean, s1i_dps_mean -> area fraction averaged normal and shear dissipation powers
// s1i_dpn_star, s1i_dps_star -> normal and shear dissipation powers from final computed star values. These are not equal to mean ones from aI, aII, aII
typedef enum { s1i_absStick = 0, s1i_absSlip = 1, s1i_absSep = 2, s1i_D, s1i_contactRel, s1i_stickRel, s1i_effectiveStrs, s1i_Dtarget, s1i_DsrcAllix, s1i_DsrcRBHV, s1i_Dsrc, s1i_fricCoef, s1i_slip_dvelTheta, s1i_slip_tauTheta, s1i_tau_minus_delv_Theta, s1i_magnitude_delv, s1i_dpn_separation, s1i_dps_separation, s1i_dps_slip, s1i_dpn_mean, s1i_dps_mean, s1i_dpn_star, s1i_dps_star, s1i_dp_star, s1i_SIZE } SolidInterfaceShortType;
// for storing Riemann solutions
typedef enum { rmode_stick, rmode_sep, rmode_slip} RiemannMode_StorageT;


// prints TSR stage of the point in the final solution files
#define TSR_STAGE_IO 0

#if DiM2a3_F
constexpr auto NUMSLRMN = 3;
#else
constexpr auto NUMSLRMN = 2;
#endif

// the first 3 are also used for labeling target values

// stores main properties of one side of the interface.
//		
//							downstream_final -> these are final star values obtained by combination of Riemann values
//
class SL_Interface_PtData_OneSide
{
public:
	void Read(istream& in);
	void SetZero();

	VEC v_downstream_final;
	VEC sigma_downstream_final;
	VEC u_downstream_final;
};

class SL_interfacePPtData
{
public:
	SL_interfacePPtData();
	void Read(istream& in);
	void Output_FinalSolution(ostream& out, IOF_type iot, double space_or_time, double x, double t, bool has_ring_opened1D_al);
	void SetZero();
	// function needed to output 1 value for the interface for visualization
	void get1DValues4Visualization(vector<double>& vals, double delC, InterfaceLocation1DT side, double t, double ax = 0.0, int unsigned dir = 0, bool periodicBoundaryPt = false);

	SL_Interface_PtData_OneSide sl_side_ptData[NUM_SIDES];
	double interface_damage_final;
	double interface_damage_source_final; // source term for damage equation
	// can add "inside" interface stuff later
	//  maximum effective delU (used for keep tracking of TSR based damage)
	double maxEffDelU;
	TSR_loadingStages tsrStage;

#if RING_PROBLEM
	// for ring problem this stores the final "converged" v_r value
	double v_r_final;
	double v_r_source_final; // p - sigma_theta/rho/R
	double sigma_theta_source_final; // E v_r / R (same on both sides)
#endif
	// this is the time of this specific point
	double interface_time;
};

// this is a class that stores time sequence of point dat (at one fixed spatial point)
// I intentionally want to make the storage private, so it can be hidden later
class SL_interfacePPtData_Time_Seq
{
public:
	SL_interfacePPtData_Time_Seq();
	~SL_interfacePPtData_Time_Seq();
	void Read(istream& in);
	void RemoveLastPoint();
	inline double GetMaxTime() const { return timeSeqPtrs[curPos]->interface_time; }
	inline bool IsEmpty() const { return (sz == 0); }
//	inline unsigned int size() const { return sz; }
	inline unsigned int Earlier_size() const { return sz - 1; }
	SL_interfacePPtData* AddNewPoint(double timeIn, unsigned int& pos);
	// relativeBackwardPos == 0, first point back relative to currect head (curPos), ...
	SL_interfacePPtData* GetBackwardPosition(unsigned int relativeBackwardPos);
	SL_interfacePPtData* GetCurrentPosition();
	bool Interpolate_Pt_Solution_At_time(double timeIn, SL_interfacePPtData*& ptSlnPtr, bool& ptSlnPtr_Deletable);

	SL_interfacePPtData_Time_Seq(const SL_interfacePPtData_Time_Seq& other);
	SL_interfacePPtData_Time_Seq& operator=(const SL_interfacePPtData_Time_Seq& other);

private:
	// if ptB == NULL, ptA has the exaxt time of timeIn and factor_ptA == 1.0
	// else ptA.interface_time < timeIn < ptB.interface_time & factor_ptA = (ptB.interface_time - timeIn)/(ptB.interface_time - ptA.interface_time), factor_ptB = 1.0 - factor_ptA
	bool Find_Pts_Around_Time(double timeIn, SL_interfacePPtData*& ptA, double& factor_ptA, SL_interfacePPtData*& ptB, double& factor_ptB);

	SL_interfacePPtData* timeSeqPtrs[SIZE_PT_TIME_SEQUENCE];
	unsigned int curPos;
	unsigned int sz;
};

class SL_Interface_Temp_PtData_OneSide
{
public:
	SL_Interface_Temp_PtData_OneSide();
	void Read(istream& in);
	void SetZero();
	// eventually downstream "final" (in SL_Interface_PtData_OneSide; computed from combination of SL_Interface_PtData_OneSide) 
	// should match "latestValue" below.
	///////////////////////////////////////////////////////////////////////////////////////////
	//
	//
	// in cSDG code, and alike "latestValue" can be interior trace and "final" is star solution. There are other variants, say for Time Marching DG methods ... 
	//
	//
	///////////////////////////////////////////////////////////////////////////////////////////
	VEC v_downstream_latestValue;
	VEC sigma_downstream_latestValue;
	VEC u_downstream_latestValue;


	// mode here refers to 
	// bonded/contact-stick , contact-slip , and separation solutions
	VEC v_mode_Star[NUMSLRMN];
	VEC sigma_mode_Star[NUMSLRMN];
	
	// these are the star values computed based on v_mode_Star, sigma_mode_Star and area fractions
	VEC v_Star;
	VEC sigma_Star;
};

// _Int -> refers to interface
// this stores important things for the interior of the interface (i.e. the shared region in between)
class SL_Interface_Temp_PtData_IntContFrac
{
public:
	SL_Interface_Temp_PtData_IntContFrac();
	void SetZero();
	void Read(istream& in);
	double interface_scalar_Vals[s1i_SIZE];

	Vec_nt sigma_I_nt_parts;
	Vec_nt del_u_nt_parts;
	Vec_nt del_v_nt_parts;
};

class SL_interface_Temp_PPtData
{
public:
	SL_interface_Temp_PPtData();
	~SL_interface_Temp_PPtData();
	void MakeReady_For_ContactFractureRuns(double beta_delU, double beta_traction, const VEC& sigmaI, double delu0Change = 0.0);
	void SetZero();
	void Read(istream& in);
	void Output_ScalarValues(ostream& out, IOF_type iot, double space_or_time);

	SLFF_InterfacialDamageModeType damageOffOnMix;
	SLFF_ContactType			   contactOffOnMix;
	SLFF_SlipType				   slipOffOnMix;

	SL_Interface_Temp_PtData_OneSide sl_side_temp_ptData[NUM_SIDES];
	SL_Interface_Temp_PtData_IntContFrac* interfacePropPtr;
	// the reason this is not stored in "interfacePropPtr"is that it can be used to make the pointer null or not
	double interface_damageLatestValue;
	TSR_loadingStages tsrStageLatestValue;

	// these are displacement normal scale, and dimensional separation, contact scales and are used to check separation-contact transition and limit interpenetration
	// dcont is moved to SLInterfaceCalculator so it can be used in controlling value of next un
	double tau, dn, dt, dsep;
	double sigmaC, deltaC;
	// can add "inside" interface stuff later
#if RING_PROBLEM
	// for ring problem this stores the most recent v_r value
	double v_r_latestValue;
	double v_r_source_latestValue; // p - sigma_theta/rho/R
	double sigma_theta_source_latestValue; // E v_r / R (same on both sides)
#endif
};

// this class stores errors (see SLFractureGlobal_Configuration items 1 and 2 for further discussion and relevant tolerances)
// 1. Iteration to iteration within one time step: step n+1, iterations k and k + 1
//				Used to continue iterations until solution for time step n + 1 is finalized
// 2. Time step to the next: step n+1, 
//				Used to adjust time step, it acts as refinement check to reduce time step for an already analyzed point


class Pt_Error_Data
{
public:
	Pt_Error_Data();
	void SetZero();
	void Read(istream& in);
	void PrintShort(ostream& out, double time, int itern);
	double norm_del_v;
	double norm_del_s;
	double del_damage;
	double damage_source_tau;
	double del_sep2cont_c;
};

void Output_ScalarValuesAux(ostream& out, IOF_type iot, double space_or_time, double* interface_scalar_Vals, SLFF_InterfacialDamageModeType damageOffOnMix, SLFF_ContactType	contactOffOnMix, SLFF_SlipType slipOffOnMix);

#endif