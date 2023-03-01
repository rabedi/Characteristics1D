#ifndef SL_BULK_PROPERTIES__H
#define SL_BULK_PROPERTIES__H

#include "globalMacros.h"
#include "globalTypesClasses.h"

double	getCdIsotropic(double rho, double nu, double E);
double	getCsIsotropic(double rho, double nu, double E);

class SL_OneInterfaceAllTimes;

////////// Characteristic Calculation
// The following 3 classes are used to calculate
//		CHARCTERISTIC vector downstream at ONE point ONE side of an interface
//					from DiM Upstream Point q's (v, sigma), source_q, and delT (each goes down in time with a certain time lag)

// this data upstream a downstream point along ONE of the characteristics lines in negative or positive x direction
class OnePoint_Data_Upstream_w_lOr_r
{
public:
	OnePoint_Data_Upstream_w_lOr_r();
	// Input entire q (v, sigma) and corresponding force (if applicable)
	// delT is how earlier this upstream point is from the downstream point for which characteristics are calculated
	double delT;
	VEC	v, sigma;
#if HAVE_SOURCE
	VEC source_v, source_sigma;
#endif
};

// collection of all upstream points (each at a different direction)
class AllPoints_Data_Upstream_w_lOr_r
{
public:
	OnePoint_Data_Upstream_w_lOr_r  up_Pts[DiM];
	VEC wUpstream;
#if HAVE_SOURCE
	VEC source_w_Upstream;
#endif
};

//
class OnePoint_Data_Downstream_w_lOr_r
{
public:
	OnePoint_Data_Downstream_w_lOr_r();

	// OUTPUT: characteristic is computed
	VEC wDownstream;

	// Note: v and sigma on this side of the interface is not known and can keep being updated.
	// It is needed if allDiagonal_Ds_are_zero == false or if other source terms depend on [v, sigma]
	VEC last_estimate_v, last_estimate_sigma;
#if HAVE_SOURCE
	// Input (and output, as the forces can be updated)
	VEC source_v, source_sigma;
#endif
};

/// this class stores u, velocity, eps, sigma, etc. for 1 point in the bulk
class OnePoint_inBulk_Fields
{
public:
	OnePoint_inBulk_Fields();
	void OnePoint_inBulk_Fields_Print_Data(ostream& out, bool print_x = true, bool print_time = false);
	static void OnePoint_inBulk_Fields_Print_Header(ostream& out, bool print_segi_pi, bool print_x, bool print_time);

	double x;
	double time;
	VEC	u, v, eps, sigma;
	double damage, delu0;
#if HAVE_SOURCE
	VEC source_v, source_sigma;
#if RING_PROBLEM
	double v_r;
#endif
#endif
private:
	// print_segi_pi: print segi and pi in class Subdomain_oneTime_spatial_points
	static void OnePoint_inBulk_Fields_Data_Get_Header_Labels(vector<HeaderLabels>& hLabels, bool print_segi_pi = true, bool print_x = true, bool print_time = false);
};


class SL_Bulk_Properties
{
public:
	SL_Bulk_Properties();
	void Initialize_FromInputParas();
	// this function is called only when even 1 factor is != 1, so instead of using precomputed bulk properties, the components are scaled
	void Initialize_FromOther_withFactors(const SL_Bulk_Properties& other, double factorC, double factorRho, double factor_damping = 1.0);
	void Form_Iso_E1Rho1Default();

	void Read_SL_Bulk_Properties(istream& in, int bulkFlag);// does not read computed parts
	void Print(ostream& out, bool printComputedVals) const;

	// q = [v, sigma] primary fields
	// w = [wL, wR]; wL is left-going characteristics (c < 0), wR right going
	// q2w means from primary to characeristic transformation
	// w2q is the opposite
	// w_index is the index of w (0 max wave speed ... to min wave speed)
	// is_wr == true (right-going wave)
	// return is wr_(w_index) (r or l)
	// it works for q OR source of q
	double q_to_w(const VEC& v, const VEC& sigma, int w_index, bool is_wr);
	// this computes the entire right or left characteristic
	void q_to_w_vec(const VEC& v, const VEC& sigma, VEC& w_rOrlGoing, bool is_wr);
	void q_to_w_2vecs(const VEC& v, const VEC& sigma, VEC& w_rGoing, VEC& w_lGoing);
	// reverse direction
	void w_to_q_2vecs(const VEC& w_rGoing, const VEC& w_lGoing, VEC& v, VEC& sigma);

	// Input: w and v are given. is_wr = true -> right going wave (i.e. left side of an interface)
	// Output: sigma
	void Compute_stress_from_v_and_characteristics(const VEC& w, const VEC& v, bool is_wr, VEC& sigma);
	// This function is used to compute downstream characteristic (left or right going) from upstream w and up/downstream forces (if any force)
	// Prerequisite (computed)
	//		1. delTs: delta times from upstream (earlier time) to current time for different points (inside upstream_pts)
	//		2. q's of all upstream points (v, sigma's for all upstream points)
	//		3. Up&Downstream q forces (if problem has source term)
	//				3a. upstream.source_v, upsteam.source_sigma (from ring & other sources BUT not q0Terms are computed)
	//				3b. downstream.source_v, source_sigma computed
	// 		useDSCharsAsUSChars: For problems with source term sometimes current v and sigma are not available we can acceptable some level of error and use upstream characterisitcs as downstream ones
	//		Output: downstreamPt.w computed
	// for periodic domain, e.g. Ring with constant vr (see section 2.3.1. in Zhou_2006_Molinari_Ramesh_Analysis of the brittle fragmentation of an expanding ring.pdf)	a constant velocity of magnitude of La (a loading rate, L domain size) is added to the velocity of the left side before Riemann solution, then the velocity is subtracted from all star values
	// this jump if present is sent as ring_opened1D_alPtr
	// x is only needed for opened up ring problem
	void Compute_Downstream_characteristic_From_Upstream_qs_and_Forces(AllPoints_Data_Upstream_w_lOr_r& upstream_pts, OnePoint_Data_Downstream_w_lOr_r& downstream_pt, bool is_wr, double x, bool useDSCharsAsUSChars = false, double* ring_opened1D_alPtr = NULL);
	// in-situ stresses:
	//	Step 1: they are added to downstream chacateristics from the two sides
	// after the function above, we need to update characteristics with in-situ stresses (is present)
	void Compute_Downstream_characteristic_From_in_situ_Stress(const VEC& in_situ_stress, const VEC& wDownstream_WO_in_situ, VEC& wDownstream_W_in_situ);
	//	Step 2: After the calculation sigma Star values, they are subtracted from star value, but since this is a simple operation independent of members of this class, it's not added here

	// strain and stress are _0i components, not the full Voigt strain and stress
	void Compute_Stress_from_Strain(const VEC& eps, VEC& sigma);
	void Compute_Strain_from_Stress(const VEC& sigma, VEC& eps);

	/////////////////////////////////////////////////////////////////////////////////////////
	// functions related to computing values at a given spaceime location
	// see Compute_Bulk_vseps_NonIC ...
	void Compute_Bulk_Values(double x, SL_OneInterfaceAllTimes *interfaceLeftOfBulk, SL_OneInterfaceAllTimes* interfaceRightOfBulk, OnePoint_inBulk_Fields& point_fieldsNT, bool isIC, OnePoint_inBulk_Fields* point_fieldsPT_or_NTPIPtr = NULL,
		unsigned int maxIter = 10, double relTol4Conv = 1e-4);

	// input
	double rho;
#if HAVE_SOURCE_ORDER0_q
	// equations are 
	//				vDot		- 1/rho sigma,x + D_vv v + D_vsigma		sigma = otherSource_v
	//				sigmaDot	-C		v,x		+ D_qv v + D_sigmasigma sigma = otherSource_sigma
	// for example if the problem has nonzero damping D_vv = damping / rho
	// values below will be provided (defaults are 0). The most likely nonzero value is D_vv
	// D_wlwl, vectors are computed from these, see below
	double D_vv, D_vsigma, D_sigmav, D_sigmasigma;
	bool has_DTerms;
#endif
	//////////////////// iso  (1D and 2D, iso and aniso)
	// input
	double E_iso;
	double nu_iso;
	GID	flag;
	// computed
	double cd_iso;
	double cs_iso;

	//////////////////// aniso  (1D and 2D, iso and aniso)
	// input
#if COMPUTE_ANISO_BULK
	MAT C;
	// computed
	////////////////////////////////////////////////////////////////////////////////////
	// E is marrix of left eigenvectors of C_over_rho: E[i][:] left eigenvalue i
	// F = (c rho) * E, where c is diagonal matrix of wave speeds (ws) and rho is density. F is related to Z
	// H = Inv(F)
	// Z = E^-1 F = E^t (c rho) E
	// Y = Z^-1
	MAT C_over_rho, E, F, H, Z, Y;
	MAT inv_C;
#endif
	// wave speeds for both cases
	VEC ws;
	// This is not called impedance as for anisotropic case Z = E^t (c rho) E
	VEC c_rhos;
	VEC effective_CsIso;
	VEC inv_c_rhos;
	VEC inv_ws;
	bool is_iso;

// computed
#if HAVE_SOURCE_ORDER0_q
	//	wl:  left-going characteristics (same as w minus)
	//	wr: right-going characteristics (same as w  plus)
	// 0th order source term in characteristics coordinate is broken down into two parts
	// Part 1: can be directly handled with exponentional change (decay) or solution along characteristics
	// e.g. wDot + c w,x + d w = ...
	VEC D_wlwl, D_wrwr;
	// Part 2: Off-diagonal terms that need to be approximated. Good solutions are obtained only in the limit of Delta-t (time advance) along characeristics -> 0
	VEC D_wlwr, D_wrwl;
	// D_omega = [diag(D_wlwl), diag(D_wlwr); diag(D_wrwl), diag(D_wrwr)]; 
	// formulas for individual components are in Initialize_FromInputParas
	// if D_wlwr, D_wrwl are zero, computations become exact and simple, if not, a current estimate of q is needed in Compute_Downstream_characteristic_From_Upstream_qs_and_Forces
	bool allDiagonal_Ds_are_zero;
#endif

private:
	// x and time of point_fields and 
	// NT: new time (x, time)
	// PT: previous time (x, time - del t)
	// NTPI: new time but the previous iteration of it (when we have source term multiple iterations may be needed)

	void Compute_Bulk_vseps_IC(OnePoint_inBulk_Fields& point_fieldsNT);

	// SOURCE term != 0
	// Why is prior/current					point_fieldsPT_or_NTPIPtr			is	needed?
	//		ONLY for source term != 0, this info is needed as the source terms are to be integrated along the characteristics
	// If point_fieldsPT_or_NTPIPtr == NULL, there is no previous time or prior info for this point and characeristics of downstream is used for upstream too (some error will be introduced)

	// ELSE
	// IF source term == 0, this is point_fieldsPT_or_NTPIPtr used at all
	void Compute_Bulk_vseps_NonIC(double x, SL_OneInterfaceAllTimes *interfaceLeftOfBulk, SL_OneInterfaceAllTimes* interfaceRightOfBulk, OnePoint_inBulk_Fields& point_fieldsNT, OnePoint_inBulk_Fields* point_fieldsPT_or_NTPIPtr = NULL);
};

// the following class is the interface between two bulks and is capable of computing Riemann I (bonded) and 
// III (or others where sigmaStar) is given
// it does NOT contain fracture properties
class SL_Elastic_InterfaceProperties
{
public:
	SL_Elastic_InterfaceProperties();
	~SL_Elastic_InterfaceProperties();

	void FormMaps();
	void Read(istream& in);// does not read computed parts
	void Print(ostream& out, bool printComputedVals) const;

	// w's are incoming characteristics 
	//										OR
	// BC for the appropriate size based on directionalBCType and side of the domain. See description of Compute_sigmaI_vI_from_wsAndBC__BoundaryCase
	void Compute_sigmaI_vI_from_ws(const VEC& wlSide_rGoing, const VEC& wrSide_lGoing, VEC& vI, VEC& sigmaI);

	// type of BC in each direction (Dirichlet, Neumann, Characteristics) in each direction is given in directionalBCType
	//		w[i]	(wlSide_rGoing for left boundary;	wrSide_lGoing for right boundary)	stores
	//		Dirichlet		(directionalBCType[i] == bct_Dirichlet)								velBar_i
	//		Neumman			(directionalBCType[i] == bct_Neumann)								tracBar_i
	//		Characteristics	(directionalBCType[i] == bct_Characteristics)						characteristicsBar_i
	//		ISOTROPIC case: each direction (i) can have a different BC
	//		ANISOTROPIC case: only the flag of the normal direction (0) is used. That is only, pure Dirichlet, Neumann, and Characteristics BCs are supported.
	void Compute_sigmaI_vI_from_wsAndBC__BoundaryCase(const VEC& w_OR_BC_lSide, const VEC& w_OR_BC_rSide, VEC& vI, VEC& sigmaI);

	void Compute_vStars_from_sigmaStar_ws(const VEC& wlSide_rGoing, const VEC& wrSide_lGoing, const VEC& sigmaStar, VEC& vStarLeft, VEC& vStarRight);

	void Compute_Downstream_characteristic_From_in_situ_Stress(bool has_in_situ, const VEC& in_situ_stress, const VEC& wlSide_rGoing_WO_in_situ, const VEC& wrSide_lGoing_WO_in_situ,
		VEC& wlSide_rGoing, VEC& wrSide_lGoing);

	// if incoming characteristics and velocity traces are given, stress traces can be computed
	// velocities are given in pointer, so if NULL, they are not provided
	// outputs are in line 2
	void Compute_stresses_from_vs_and_characteristicss(const VEC& wlSide_rGoing, const VEC& wrSide_lGoing, const VEC* vlPtr, const VEC* vrPtr,
		VEC& sigmal, VEC& sigmar);
	// this time the traces of v,s are given for either side and characteristics are calculated
	void Computecharacteristics_from_Traces_vel_stress_etc(const VEC* traceCurrentStress_WO_in_situ_LeftIn, const VEC* traceCurrentStress_WO_in_situ_RightIn, const VEC* traceCurrentVelLeftIn, const VEC* traceCurrentVelRightIn, VEC& wlSide_rGoing, VEC& wrSide_lGoing);

	void Get_Bulks2Side(SL_Bulk_Properties*& bulk_left_outPtr, SL_Bulk_Properties*& bulk_right_outPtr);
	InterfaceLocation1DT interfaceLoc;
	// directional BC type is only relevant for left and right side of the domain. For interior, this is not used
	BoundaryConditionT directionalBCType[DiM];

	// left side: characteristic move in the positive direction (wlSide_rGoing)
	// right side: haracteristic move in the negative direction (wrSide_lGoing)
	SL_Bulk_Properties *bulk_leftPtr, *bulk_rightPtr, *bulk_inside_domain4BoundaryPtr;
	bool bulk_leftPtr_deletable, bulk_rightPtr_deletable;
	bool isIsoInterface;

	/// maps from characeristics to Riemann solutions for v and sigma 
	/////////////////
	// iso
	/// I bonded
	// sigmaI = iso_wlSide_rGoing_2_sigmaI * wlSide_rGoing + iso_wrSide_lGoing_2_sigmaI * wrSide_lGoing;
	VEC iso_wlSide_rGoing_2_sigmaI;
	VEC iso_wrSide_lGoing_2_sigmaI;
	// vI = iso_wlSide_rGoing_2_vI * wlSide_rGoing + iso_wrSide_lGoing_2_vI * wrSide_lGoing;
	VEC iso_wlSide_rGoing_2_vI;
	VEC iso_wrSide_lGoing_2_vI;
	/// III (and others with specified sigma*)
	// vStar_lSide = iso_sigmaStar_2_vStar_lSide * sigmaStar + iso_wlSide_rGoing_2_vStar_lSide * wlSide_rGoing
	// vStar_rSide = iso_sigmaStar_2_vStar_rSide * sigmaStar + iso_wrSide_lGoing_2_vStar_rSide * wrSide_lGoing
	VEC iso_sigmaStar_2_vStar_lSide;
	VEC iso_sigmaStar_2_vStar_rSide;
	VEC iso_wlSide_rGoing_2_vStar_lSide;
	VEC iso_wrSide_lGoing_2_vStar_rSide;
	// these are the averages of normal and shear Y's from the two sides of the interface. They are used to compute displacement scales of RBH damage evolution model
	double YAve_n, YAve_t;

#if !USE_ISO_ASSUMPTION
	// aniso ones are similar but the maps are matrices
	MAT wlSide_rGoing_2_sigmaI;
	MAT wrSide_lGoing_2_sigmaI;
	MAT wlSide_rGoing_2_vI;
	MAT wrSide_lGoing_2_vI;

	MAT sigmaStar_2_vStar_lSide;
	MAT sigmaStar_2_vStar_rSide;
	MAT wlSide_rGoing_2_vStar_lSide;
	MAT wrSide_lGoing_2_vStar_rSide;

	MAT Yl_plus_Yr, Zl_plus_Zr, inv_Yl_plus_Yr, inv_Zl_plus_Zr;

#if DiM2a3_F
	double sumYnn;
	VECm1 sumYnt, sumYtn;
	MATm1 sumYtt, sumYtt_eff;
	MATm1 alphaMat; // = Inv(sumYtt_eff)
	// matrices for friction
	VECm1 gammaVec; // = sumYnn^-1 Ynt
	VECm1 betaVec; // = gammaVec * alphaMat
#endif

#endif
	SL_Elastic_InterfaceProperties(const SL_Elastic_InterfaceProperties& other);
	SL_Elastic_InterfaceProperties& operator=(const SL_Elastic_InterfaceProperties& other);
};


#endif
