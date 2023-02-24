#ifndef DOMAIN_POST_PROCESSING__H
#define DOMAIN_POST_PROCESSING__H

#include "SLBulk_Properties.h"
#include "SL_Interface_PtData.h"

//////////////////////////////////////////
// This file is Level 1
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


// Definitions:
//			Subdomain (sdi): it's a collection of SL_Bulk_Properties that forms one subdomain of the entire domain
//						e.g. in projectile target example, projectile and target are each one a subdomain

//					sd_segment (segi): a segment within a subdomain, it is this class			SL_Bulk_Properties

//							seg_pt (pti): these are the points within a segment
//	Indices for these are shown in paranthesis

typedef enum {sdit_left, sdit_interior, sdit_right, subdomainInterfaceType_SIZE} subdomainInterfaceType;

#define NUM_AVE_SUM_TIMES 3

// STORAGE of SPACIAL or SPACETIME Integrals or average for 
//				ONE subdomain
// ave: average					a_ave = 1/|V| int_space a dV
// sum: integral over space		a_sum =		  int_space a dV
class Subdomain_spatial_ave_sum
{
public:
	Subdomain_spatial_ave_sum();
	void Subdomain_spatial_ave_sum_set_zero();

//	data print and header
	void Subdomain_spatial_ave_sum_Read_Data(istream& in);
	void Subdomain_spatial_ave_sum_Write_Data(ostream& out);
	static void Subdomain_spatial_ave_sum_Data_Write_Header(ostream& out);

	int timeIndex;
	double timeVal;
	////////////////////////////////////////////////////////////
	// Below are quantities that are integrated (and averaged) over temporal domain from t = 0 -> time
	// BC-based definition of strain and stress averages
	VEC eps_ave_bc, sigma_ave_bc;
	VEC eps_ave_bc_minus_bulk, sigma_ave_bc_minus_bulk;
	double energy_eps_total_bc, energy_eps_diss_bc, energy_eps_recoverable_bc;

	// BC linear momentum
	// stresses from left and right
	VEC u_L, u_R;
	VEC sigman_L, sigman_R;
	// the integrals of stresses in time -> impulse
	VEC impulse_L, impulse_R, impulse_BC;

	// BC energy
	// powers on the left and right side of the domain | power = sigma . v
	double power_L, power_R;
	// energy input from left, right, and both sides
	double energy_L, energy_R, energy_BC; 

	// interior contact/fracture surfaces
	// I stands for Interface
	VEC powerIDiss_vec;
	VEC eneIDiss_vec; // integral of above
	double powerIDiss, eneIDiss; // these are the time integrals of the above
	double energyIDiss_Recoverable; // 0.5 * delta U * sigma is the recoverabe part
	double phys_diss_interface_lost;// = eneIDiss - energyIDiss_Recoverable
	double max_interface_damage, min_interface_damage, mean_interface_damage, sdiv_interface_damage;
	double max_interface_damage_source;
	double max2EffDelU;

	////////////////////////////////////////////////////////////
	// Below are quantities that are integrated (and averaged) over spatial domain
	// eps_ave_bulk = eps_ave_bulk_intfc + eps_ave_bulk_blk, interface is the strain from interface displacement jumps and eps_ave_bulk_blk is raw bulk strains
	VEC	u_ave, v_ave, eps_ave_bulk, eps_ave_bulk_intfc, eps_ave_bulk_blk, sigma_ave_bulk;
	double energy_eps_total_bulk, energy_eps_diss_bulk, energy_eps_recoverable_bulk;

	// energy and linear momentum density stuff
	// rho vDot - sigma,x = rho source_v
	// sigmaDot - C v,x = source_sigma

	// Linear momentum balance law: TEMPORAL: (rho v)			SPATIAL: sigma		SOURCE: rho source_v
	// Energy:						TEMPORAL: phi = K + U		SPATIAL: v sigma    SOUCE: source_e_v + source_e_sigma
	//				K = 0.5 * rho * v * v,			U = 0.5 * sigma * eps
	//				source_e_v = v (rho v)			souce_e_sigma = sigma * C^-1 * source_sigma = eps * source_sigma (since C is symmetric)

	VEC linMomentum_dsum;
	double K_dsum, U_dsum, phi_dsum;

	////////////////////////////////////////////////////////////
	// Below are quantities that are integrated over spacetime up to time
	// dsum are domain sum (not physically important but needed for calculating spacetime ones)
	// dtsum are space time ones
#if HAVE_SOURCE
#if RING_PROBLEM
	double v_r_ave;
#endif
	// source term (other than damping) contribution to linear momentum and energy
	// space:
	VEC source_linMomentum_dsum;
	double source_e_v_dsum, source_e_sigma_dsum, source_e_dsum;
	// spacetime (real important one)
	VEC source_linMomentum_dtsum;
	double source_e_v_dtsum, source_e_sigma_dtsum, source_e_dtsum;

	// damping contribution to linear momentum and energy
#if HAVE_SOURCE_ORDER0_q
	// space:
	VEC damping_linMomentum_dsum;
	double damping_e_v_dsum, damping_e_sigma_dsum, damping_e_dsum;

	// spacetime (real important one)
	VEC damping_linMomentum_dtsum;
	double damping_e_v_dtsum, damping_e_sigma_dtsum, damping_e_dtsum;
#endif
#endif

	// numerical errors
	// phi_t = phi_0 + Energy_BC + Energy_source  - Energy_damping - Energy_Inerface_diss - Energy_numerical_diss
	// -> Energy_numerical_diss = (phi_0 + Energy_BC + Energy_source) - (Energy_damping + Energy_Inerface_diss)	- phi_t 	
	double input_energy; // = phi_0 + Energy_BC + Energy_source
	double phys_diss_tot; // = Energy_damagping + Energy_Inerface_diss
	double phys_diss_lost;// = phys_diss_tot - energyIDiss_Recoverable
	double numerial_energy_diss; // = input_energy - phi_t - phys_diss_tot
	// similar for linear momentum with no contribution from fracture interfaces (stress is continuous)
	VEC numerical_linmom_error;


	// this is a temporary vector to compute energy dissipation of domain by storing sigma and delu vectors
	vector<pair<VEC, VEC> > interface_sigma_delu_vec;
	void Size_interface_sigma_delu_vec(unsigned int numInterfaces);

private:
	static void Subdomain_spatial_ave_sum_Data_Get_Header_Labels(vector<HeaderLabels>& hLabels);
};

// this class can be used to refer to and easily access bulk and interfaces on the sides of one segment in post-processing. This Is used in Domain_AllInterfacesAllTimes
class OneBulktwoSideInterfaceInfo
{
public:
	OneBulktwoSideInterfaceInfo();
	SL_Bulk_Properties *bulkPtr;
	SL_OneInterfaceAllTimes *interfaceLeftOfBulkPtr, *interfaceRightOfBulkPtr;

	// cntrs refer to Domain_All_Interfaces_All_Times::bulks (1st line) and ::interfaces (2nd line)
	unsigned int bulk_cntr;
	unsigned int interfaceLeftOfBulk_cntr, interfaceRightOfBulk_cntr;
};

class OneSubdomain_All_bulksConnectivityInfo
{
public:
	void PrintIndicesLengthsKeyRunParameters(ostream& out) const;
	vector<OneBulktwoSideInterfaceInfo> subdomain_bulk_segments;
	vector<SL_OneInterfaceAllTimes*> subdomain_interfaces;

	unsigned int subdomain_number;

	// the min and max of bulk_cntrs inside all subdomain_bulk_segments
	unsigned int bulk_cn_st, bulk_cn_en;
	int numSegments;
	vector<double> subdomain_interface_xs;


	/////// computed
	double xm;
	vector<double> segment_lengths;
	double xM, length, inv_length;
	void OneSubdomain_All_bulksConnectivityInfo_Initialize();
};

// ONE subdomain
//			ONE time
// stores spatial point values 
class Subdomain_oneTime_spatial_points
{
public:
	void Subdomain_oneTime_spatial_points_Print_Header_Data(ostream& out, bool printHeader = true);
	// this is indexed as spatialPoints[segi][pti], where segi is the segment (SL_Bulk_Properties) number in the subdomain and pti is the point number within the segment
	vector<vector<OnePoint_inBulk_Fields> > spatialPoints;
};

// This stores all post-processed data for one subdomain
class Subdomain_spacetime_pp_data
{
public:
	Subdomain_spacetime_pp_data();
	~Subdomain_spacetime_pp_data();
	void Initialize_Subdomain_spacetime_pp_data(unsigned int numTimesIn, OneSubdomain_All_bulksConnectivityInfo* sdciPtrIn, bool uniform_deltIn, int numTimeStep_BulkInterfacePoints_PrintIn = -20, int numTimesStep4Print_SegmentIn = -20,
		int numSpatialSubsegments_BulkInterfacePoints_PrintIn = 2, bool useRepeatedSimpsonRuleForHigherOrders = false,
		unsigned int maxIterIn = 10, double relTol4ConvIn = 1e-4);
	void AddComputeTimeStep(int timeIndexIn, double timeValue);

	//////////////////////////////////////////////////////////////////////////
	// provided from outside
	OneSubdomain_All_bulksConnectivityInfo* sdciPtr;
	unsigned int numTimes;
	unsigned int maxIter;
	double relTol4Conv;
	bool uniform_delt;
//	int numTimes;
	// every numTimesStep4Print steps spatial data is printed
	int numTimeStep_BulkInterfacePoints_Print;
	int numTimeStep_Interface_DSU_Fragment_Print;
	// how many segments are there in each segment (SL_Bulk_Properties)
	int numSpatialSubsegments_BulkInterfacePoints_Print;
	///////////////////////////////////////////////////////////////////////////
private:

	// this is sized two in reverse order (0 is current time step, 1 prev one, ...)
	vector< Subdomain_oneTime_spatial_points*> allSpatialPointsRevOrder;
	int allSpatialPointsRevOrder_size;

	// index by time index in reverse order (current time -> 0, previous time 1, etc)
	Subdomain_spatial_ave_sum* spatial_ave_sums[NUM_AVE_SUM_TIMES];
	double linMomentumZero;
	double phi0;

	//////////////////////////////// times and number of times
	vector<double> times;

	//////////////////////////////// segments (SL_Bulk_Properties) and their number
	// indexed as xs[segi][[pti] are the x coordinate of the domain
	vector<vector<double> > xs, x_weights;


	//////////////////////////////// set of points and subsegments within each segment

	// = numSpatialSubsegments_BulkInterfacePoints_Print + 1, the one above is given to the class
	int numSpatialPointsPerSegment;
	// below is of size numSpatialPointsPerSegment and stores the Newton-Cotes weights for a segment (the actual integration weights are these times segment_lengths[segi])
	vector<double> spatialIntegrationPoints, spatialIntegrationWeights;

	int timeIndex;
	// time integration weights
	// these are the weights and size of values that multiplied by velocity from current time step backward
	// rti means index zero is current time step, 1 one back, ...
	vector<double> timeIntWeights_rti;
	unsigned int timeIntWeights_rti_size;
	double delt;

	bool print_space_points, print_fragmentation;

	///////////////////////////////////////////
	// related to initialization stage
	// prerequisite: 	xm, numSegments, segment_lengths -> xM, xs, length, inv_length
	void Set_xs();

	///////////////////////////////////////////
	// Point computation
	// timeIndex is already known
	// 0 < pti < numSpatialSubsegments_BulkInterfacePoints_Print (so it's inside the segment)
	void Compute_Inside_Segment_pt(int segi, int pti, SL_Bulk_Properties& segment);
	// interfacei == 0 -> left side
	// interfacei == numSegments -> right side
	// else interior interface with fracture
	subdomainInterfaceType Compute_End_Segment_pt(SL_OneInterfaceAllTimes* subdomain_interface, int interfacei, SL_interfacePPtData& ptSln, SL_Bulk_Properties* segmentPtrLeft, SL_Bulk_Properties* segmentPtrRight);

	// NoInterfaceParts: means dissipation at interior interfaces and left and right side of the domain (contributing to BC linear momentum and energy) are not handled in this function
	// timeIndex is already known
	void Update_Domain_space_spacetime_Integrals_fromPoint_NoInterfaceParts(int segi, int pti, SL_Bulk_Properties& segment);
	void Integrate_v4u_nonIC(int segi, int pti);

	/////////////////
	// Finalizing computations and printing
	// this is called after all bulks and interfaces of the domain are calculated (2 functions above)
	void Finalize_Subdomain_spatial_ave_sum_One_TimeStep();
	void Print_Interface_DSU_Fragment_OneTimeStep();

	std::fstream* out_sd_summary;
};

#endif