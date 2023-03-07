#ifndef SLDESCRIPTORDATA__H
#define SLDESCRIPTORDATA__H

#include "globalMacros.h"
#include "globalTypesClasses.h"
#include "FourierSimple1D.h"

// this is a class for computing:
//	1. IC
//	2. BC
//	3. Source term

//	if Ring problem				RING_PROBLEM 1 
//	ring_R
//	ring_p
//  are used
//	if AXT problem				LOAD_NUMBER_AXT	1 (u_i = a_i xt is used)
//	else						load number is used

#define AXT_LN	1
// cyclic loading, tensile, compressive
#define LC_CYCLIC_TEN_COMP	2
#define SQUARE_PULSE	3

// piece-wise loading paras are
// 0 -> t0, 1 -> l0
// 2 -> t1, 2 -> l1
// ...
// load values are l0, ....
// after final t solution is zero
#define LC_PIECE_WISE_LIN 4 

class SL_Bulk_Properties;

class SL_Elastic_InterfaceProperties;

class SLDescriptorData
{	
public:
	SLDescriptorData();
	~SLDescriptorData();
	void Finalize_SLDescriptorData();
	void Read(string configNameIn);

	// return: true if source term is nonzero
	// souce_v and _sigma are source terms corresponding to v and sigma equations.
	//	They do NOT include Ring source terms and zeros order term (e.g. D_vv, ... -based source terms in SL_Bulk_Properties
	bool GetNonRingNonZeroTerm_SourceTerm(double x, double t, GID bulkFlag, VEC& source_v, VEC& source_sigma) const;
	double GetRing_p(double x, double t);
	/////////////////// general problem with multiple interfaces
	// eps = (u0,0, u1,0, u2,0) = (gamma_0, gamma_1, gamma2) = (eps_00, 2 eps_01, 2 eps_02 (3D)) only contains strain components that are nonzero
	// damage0 can be nonzero
	// damage_source_0 is often 0

	// the function may be called for 
	//			bulk (one region) in which v_0 and eps_0 are needed
	//			interface (between two bulk regions) in which case all quantities are needed.
	//			flags are sent to the function to allow tunning the values
	void Get_IC(SL_Bulk_Properties* bulk_Ptr, GID bulkFlag, double x, VEC& u_0, VEC& v_0, VEC& eps_0
#if RING_PROBLEM
	, double& v_r
#endif
	) const;
	// this is the actual boundary condition value: ts_bulkProps is the interface
	void GetBoundaryConditionValue_LeftRight(bool isLeft, SL_Elastic_InterfaceProperties* ts_bulkProps, double x, double t, VEC& BC_val) const;

	/////////////////// One interface problem
	// this function is useful for one interface problems wherein characteristics imping from left and right on it
	void Get1InterfaceImpingingCharacteristics(SL_Elastic_InterfaceProperties* ts_bulkProps, double t, VEC& wlSide_rGoing, VEC& wrSide_lGoing);
	double GetLoadingTimeScale() const;
	void Finalize_tdLoadParameters(double finalTime, double delT, double amientProjLength, double Ein, double rhoin, double Eout = -1, double rhoout = -1);
	void Print_tdLoad();
	void Read_tdLoad();

	// axt parameters
	VEC a_xt_prob;
	// ring problem parameters
	unsigned int load_number;
	unsigned int sz_load_parameters;
	vector<double> load_parameters;

	//////////////////////////////////////// incident / impact / etc.
	// this is used to enforce incident, impact, and time-dependent Dirichlet, Neumann BC
	// read
	TD_FD* tdLoadComputer;
	LoadModeType tdLoadType;
	int tdLoadSide; // SDL -> wave goes from left to right; SDR -> wave goes from right to left
	double tdLoad_stressScale, tdLoad_velScale;

	// computed:
	// Z inside and outside correspond to Z's of main domain (part being investigated) and ambient/projectile (Zo) - For Dirichlet/Neumann Zo = Zi
	// ZEffective = (2 Zi Zo) / (Zi + Zo) -> impact, Zi -> Dirichlet & Neumann (Zi + Zo)^2 / 4Zi -> incident
	double tdLoad_Zi, tdLoad_Zo, tdLoad_ZEffective;
	double bndryLoad_inputEnergy, tdLoad_loadScale, tdAmbientProjectileLength;
	bool b_tdLoad_stressScale, b_tdLoad_velScale;

private:
	void Read(istream& in);
	double ring_p0;
};


// data members moved from SLDescriptorData to other classes:
// -> Domain_All_Interfaces_All_Times: directionalBCTypeLeftSide, directionalBCTypeRightSide, rign_R (L/2/PI): use g_domain

// ring_p -> is changed to a function that takes space-time coordiante

class CyclicLoading
{
public:
	void Initialize_FromParamaters(const vector<double>& paras);
	double getValue(double time);
	double sigma0;	//magnitude of applied load								- paras[0]
	double T;	// period of sin wave										- paras[1]
	double beta;	//to compute the magnitude of uniform compression		- paras[2] (if available)
	double tr;	//the relative ramp time (w.r.t T)							- paras[3] (if available)
};

extern SLDescriptorData g_SL_desc_data;


#endif