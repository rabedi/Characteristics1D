#ifndef SL_ONE_INTERFACE_ALL_TIMES__H
#define SL_ONE_INTERFACE_ALL_TIMES__H

#include "SLInterfaceCalculator.h"
#include "globalTypesClasses.h"


// This class contains the time sequence of solutions at ONE spatial position

class SL_OneInterfaceAllTimes
{
public:

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// 1-interface problem
	// returns 
	//			1 : If reaches the terminal condition (e.g. final time, max damage, ...); 
	//			0 : If the run fails to reach the final time
	//			-1: If neither condition is met for a very large number of time steps.
	int Main_One_InterfaceProblem();

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// main functions

	// functions for computing
	// 1. Initial step
	// this function sets the very first storage in timeSeqData, using ICs
	void InitialStep_Use_IC(SLInterfaceCalculator& slic);
	// 2. Next steps
	AdaptivityS NonInitialStep(double deltaT, bool &accept_point, double& maxTime, int currentTimeIndex, SLInterfaceCalculator& slic);

	void Set_EF_Properties(SL_Interface_Fracture_PF* interfacePFsIn, SL_Elastic_InterfaceProperties* ts_bulkPropsIn = NULL);
	double getDeltaC() const;
	double getSigmaC() const;
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// auxiliary functions

	SL_OneInterfaceAllTimes();
	~SL_OneInterfaceAllTimes();
	// Inputs: 

	//				Fnterface location and BC flags (only for left/right interfaces): 	interfaceLoc, directionalBCType

	//				Bulk:	left and right side bulks
	//						ts_bulkPropsDeletableIn:		whether elastic interface formed is deletable
	//						if either bulkLeft or bulkRight is inhomogeneous or NULL (i.e. boundary of the domain), ts_bulkPropsDeletableIn is set to true
	//							interfacePFs is formed inside the function and is not set equal to ts_bulkPropsInOut, interfacePFs_Deletable becomes true
	//						else
	//							interfacePFs = ts_bulkPropsInOut (if ts_bulkPropsInOut == NULL, it will be formed). In this case interfacePFs_Deletable becomes false

	//				Side time sequences: left_allTimes, right_allTimes Not null values provided for the side that has a neighbor time sequence

	//				periodic_totalLength : for periodic domain (e.g. ring problem) length needs to be sent in to calculate the correct distance on the sides of the domain

	//	Outputs: based on the discussion above
	//			ts_bulkPropsInOut can be an output ONLY if 1. ts_bulkPropsDeletableIn == false, and 2. ts_bulkPropsDeletableIn = NULL

	//			also	sides_timeSeqData, sides_x, sides_delx, sides_delts are set

	void Set_left_right_TimeSequenceData(InterfaceLocation1DT interfaceLoc, BoundaryConditionT directionalBCType[],
											SL_Bulk_Properties* bulkLeft, SL_Bulk_Properties* bulkRight, 
											bool ts_bulkPropsDeletableIn, SL_Elastic_InterfaceProperties*& ts_bulkPropsInOut,
											SL_OneInterfaceAllTimes* left_allTimes = NULL, SL_OneInterfaceAllTimes* right_allTimes = NULL, 
											double* periodic_totalLength = NULL);

	void Open_fixed_x_files_SL_OneInterfaceAllTimes(IOF_type iofFinalSolutionIn, IOF_type iofScalarValsIn, int subdomainLeft = -1, int subdomainRight = -1);

	// current values of downstream source and q are obtained from
	//  1) current_ptData != NULL	itern == 0 for SLInterfaceCalculator				current_ptData			
	//  2) current_ptData == NULL	itern >  0 for SLInterfaceCalculator				last previously solved permanent point from timeSeqData	if 
	void Compute_DownStream_Characteristics_wo_in_situ(VEC& wlSide_rGoing_WO_in_situ, VEC& wrSide_lGoing_WO_in_situ,
		double currentTime, int currentTimeIndex, const SL_interface_Temp_PPtData* current_ptData = NULL);

	double getEffectiveDamage_4_InterfacialDamage_TSRs(const SL_interfacePPtData& pt);

	void Set1DOrtizType();
	void Print(ostream& out) const;

	// inputs
	// for periodic domain, e.g. Ring with constant vr (see section 2.3.1. in Zhou_2006_Molinari_Ramesh_Analysis of the brittle fragmentation of an expanding ring.pdf)
	// a constant velocity of magnitude of La (a loading rate, L domain size) is added to the velocity of the left side before Riemann solution, then the velocity is subtracted from all star values
	void Set_ring_opened1D_left_side_jump_handling_true();
	bool has_ring_opened1D_al;

	int interface_pos;
	double interface_x;
	GID interface_flag;
	SL_Elastic_InterfaceProperties*	ts_bulkProps;
	SL_Interface_Fracture_PF* interfacePFs;
	// an implicit (backward Euler) solution procedue for Ortiz modle is given in TSR1D.* files. If this boolean is true, such method is going to be use
	bool only1DOrtizModel;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if there are no sides_timeSeqData (below), the interface is not surrounded by any other interfaces and it's much easier to obtain downstream characteristics
	bool b_single_interfaceProblem;
	// Set by Set_left_right_TimeSequenceData
	// this interface can have left and right neighbors or none (for example an interface in an infinite domain)
	// for finite size domains boundaries of the domain only one is provided (for example the right one for the left boundary)
	// Since the pointers are coming from outside, they are not deleted by destructor
	SL_interfacePPtData_Time_Seq* sides_timeSeqData[NUM_SIDES];
	double sides_x[NUM_SIDES];
	double sides_delx[NUM_SIDES];
	VEC sides_delts[NUM_SIDES];
	// stores the position of the bulk index
	int sides_bulk_index[NUM_SIDES];
	// stores distint subdomains of this point
	vector<unsigned int> subDomainNos, relPos_wrt_subDomainStartPoints;
	unsigned int sz_subDomainNos;
	// for computing incident wave. If not left or right, nothing is done
	InterfaceLocation1DT incidentSide;

	// will be computed
	bool ts_bulkProps_Deletable;
	bool interfacePFs_Deletable;
	SL_interfacePPtData_Time_Seq timeSeqData;

	double min_delT, max_delT; // minimum delta T computer from either side for domain problems

	// g_SL_desc_data is used for this
	double sigmaCFactor;
	double deltaCFactor;
	double iniDamage;
	IOF_type iofFinalSolution, iofScalarVals;
	ostream *outScalars, *outFinalSln, *outAdaptivity, *outIterationConv;
	SL_OneInterfaceAllTimes(const SL_OneInterfaceAllTimes& other);
	SL_OneInterfaceAllTimes& operator=(const SL_OneInterfaceAllTimes& other);

	int Main_One_InterfaceProblem_Aux();
};

// this is useful in generating contour plots of data
class OneVisContour_xInfo
{
public:
	SL_OneInterfaceAllTimes* interfacePtr;
	unsigned int subdomainNo;
	unsigned int interface_index;
	double interface_x;
	InterfaceLocation1DT side4output;
};

void MAIN_SL_OneInterfaceAllTimes_ONE_Interface(string configNameIn = "config/OneInterface/SampleConfig.txt");


#endif