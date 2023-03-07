#ifndef DOMAIN_ALL_INTERFACES_ALL_TIMES__H
#define DOMAIN_ALL_INTERFACES_ALL_TIMES__H

// this is the domain for characteristics method

// Main output: vector of SL_OneInterfaceAllTimes, one for each spatial position
// Main input: Subdomain_ElasticFractureModifier for each subdomain, the sequence of these instructions form the whole domain

#include "Subdomain_ElasticFractureModifier.h"
#include "SL_OneInterfaceAllTimes.h"
#include "DomainPostProcessing.h"

class Domain_All_Interfaces_All_Times
{
public:
	Domain_All_Interfaces_All_Times();
	~Domain_All_Interfaces_All_Times();
	// Reading the domain
	void Read_Initialize(string configNameIn, int serialNumberIn = -1);
	int Main_Domain_Solution();
	void Print(ostream& out) const;

	////////////////////////////////////////////////////////////////
	//						Inputs
	//////////////////////	serial number
	int serialNumber;
	// -1, fracture on, off decided based on run parameter, 0 turns it off
	// read from file
	int fractureMode;
	// set after reading all interfaces
	bool hasFracture;
	//////////////////////	Geometry
	vector<string> subdomain_config_names;
	// main subdomain number
	unsigned int mainSubdomainNo;

	// possible other read ins
	// the user may skip providing x_min and individual subdomain / 1 unit sizes. If so, they'll be obtain from individual subdomains, if not, the subdomain's start and end are computed directly
	double x_min;
	// each subdomain individually can be repeated Subdomain_ElasticFractureModifier::numRepeatSequence. The size of one of those may be provided below
	vector<double> subdomain_one_part_size;

	//////////////////////	Bulk Elastic properties
	map<GID, SL_Bulk_Properties*> bulk_elastic_map;
	//////////////////////	Interface Fracture properties
	// flag 0 will always be filled with no fracture interface (fully bonded one)
	map<GID, SL_Interface_Fracture_PF*> interface_fracture_map;

	//////////////////////	Left and Right BC types
	// for the ring problem, it will automatically turned to true
	bool isPeriodic;
	// type of left and right boundary conditions for (semi)-finite problems
	BoundaryConditionT directionalBCTypeLeftSide[DiM];
	BoundaryConditionT directionalBCTypeRightSide[DiM];
	// by checking the flags above, it determines if either side uses symmetric / anti-symmetric BCs
	// if so, for that side, we don't have to use 100% bonded interface condition
	bool hasSymOrAntiSymBC[NUM_SIDES];

	////////////////////////////////////////////////////////////////
	// post-process space / spacetime data and print it (summary of energies) space solutions at certain time steps
	bool do_space_spacetime_PP;
	// when source term is nonzero to get the solution at a point inside a segment (between interfaces) iterations are needed
	// pointSolution_maxIter_4PP -> what is the maximum of iterations needed
	// pointSolution_relTol4Conv_4PP -> relative tolerance between previous solution and current solution to claim convergence
	unsigned int pointSolution_maxIter_4PP;
	double pointSolution_relTol4Conv_4PP;
	// how often space solution is printed (timeIndex % numTimesStep_4PP_Print == 0) 
	// cases:
	//	numTimesStep_4PP_Print > 0
	//			e.g.	numTimesStep_4PP_Print = 1 (all time steps printed)
	//					numTimesStep_4PP_Print = 10 (every time step is printed (step 0, 10, 20, ...)
	//	numTimesStep_4PP_Print = 0
	//					only the last step is printed
	//	numTimesStep_4PP_Print < 0
	//					the time sequence is divided into (-numTimesStep_4PP_Print) segments
	//			e.g.	numTimesStep_4PP_Print = -10, numTimes = 400 -> 
	//						numTimesStep_4PP_Print = 400/10 = 40, so everytime step is printed (with a total of 11 prints)
	// first one is for detailed output of space (intermediate points within the bar), the second one is only brief data for [Damage, delU] between segments
	int numTimeStep_BulkInterfacePoints_Print_4PP, numTimeStep_Interface_DSU_Fragment_Print_4PP;

	//	for spatial segments (interval between interfaces) how many segment(s) is/are used.
	//			1 -> only the two end points at the side interfaces that already exist are used.
	//			2 -> 2 segments -> 3 points, 2 already on the sides + 1 in the middle.
	//			n -> n + 1 points are used with 2 already existing on the two sides
	int numSpatialSubsegments_BulkInterfacePoints_Print_4PP;

	// For higher orders (numSpatialSubsegments_BulkInterfacePoints_Print >= 3), it uses repeated Simpson (+Trapezoidal rule)
	bool useRepeatedSimpsonRuleForHigherOrders_4PP;

	// output option
	////// all space for one time
	IOF_type io_type_InterfaceRawFinalSln_AllSpace_Print_4PP;
	IOF_type io_type_InterfaceRawScalars_AllSpace_Print_4PP;
	// The number of time steps at which we write raw data (scalar values, final solutions). 
	// This is different from numTimeStep_BulkInterfacePoints_Print_4PP in that the former prints raw data, and latter prints processed data even for points inside the element
	// = 0 -> inactive
	// > 0 -> every this number of steps data is written
	// < 0 -> -numTimeStep_InterfaceRawFinalSlnScalars_AllSpace_Print_4PP time steps overall chosen (e.g. -20 total of 20 time steps are used)
	int numTimeStep_InterfaceRawFinalSlnScalars_AllSpace_Print_4PP;
	// for one time multiple space output (corresponding to numTimeStep_InterfaceRawFinalSlnScalars_AllSpace_Print_4PP) whether to print scalar values or not

	////// all time for one space
	IOF_type io_type_InterfaceRawFinalSln_AllTime_Print_4PP;
	IOF_type io_type_InterfaceRawScalars_AllTime_Print_4PP;
	// this is similar to above, but this time data is printed for a fixed time value
	int numSpaceStep_InterfaceRawFinalSlnScalars_AllTime_Print_4PP;

	////////////////////////////////////////////////////////////////
	/// 1D visualization files generated
	int visualization_dir;
	int visualization_numSpaceStep;
	int visualization_numTimeStep;

	////////////////////////////////////////////////////////////////
	//						Computed
	// maximum x and length
	double x_max, L;
	// for ring problem
	double ring_R;
	// subdomains
	unsigned int num_subdomains;
	vector<Subdomain_ElasticFractureModifier*> subdomains;

	// interfaces and bulks that are used in its creation
	// bulk size = num_interfaces - 1
	unsigned int num_bulks;
	vector<SL_Bulk_Properties*> bulks;
	// some of the bulks directly take the pointer from bulk_elastic_map so they are not deletable, some are created by new so are deletable
	vector<bool> bulks_deletable;
	// bulk-bulk interfaces: All will be deletable
	map<pair<GID, GID>, SL_Elastic_InterfaceProperties*> lr_IDS_2_ts_bulkProps;

	unsigned int num_interfaces;
	vector<double> interface_xs;
	vector<SL_OneInterfaceAllTimes*> interfaces;

	// min delT for interface calculations over all interfaces
	double min_domain_del_t;
	double max_domain_del_t;

	///////////////////////////////////////////////////////////////////////////////////////////
	// data needed for postprocessing
	// this is a vector of size num_bulks + 1 and pos i, i + 1 contains the start and end bulks numbers for subdomain i
	vector<int> subdomain_bulk_start_nos;
	// includes bulks and interfaces related to individual subdomains
	vector<OneSubdomain_All_bulksConnectivityInfo> bulk_interfaces_subdomains;
	// postprocessing data for each subdomain
	vector<Subdomain_spacetime_pp_data> postProcessing_subdomains;
	vector<int> subdomainNo4AllBulks;

	///////////////////////////////////////////////////////////////////////////////////////////
	// opened up ring problem
	// In Zhou_2006_Molinari_Ramesh_Analysis of the brittle fragmentation of an expanding ring.pdf, section 2.3.1. 
	// the ring problem is opened up and treated as a 1D problem with specific BC (eqs 11 & 12)
	// v = vTheta + ax, a = vr/R 
	// b_ring_opened1D = true for such treatment
	bool b_ring_opened1D;
	//	this one when on, turns on fracture on the end point of the periodic domain to avoid complications pertained to modeling fracture on this boundary
	bool b_ring_open_turn_fracture_on_periodic_end;
	// below: true	-> damping term / rho = Dvv * vTheta 
	//			(physically correct, but for elastic response, vTheta = 0, and for fracture this is mostly zero, so damping for ring problem is a 2nd order effect with this set-up)
	//		  false: not recommended (non-physical) -> damping / rho = Dvv * v
	bool b_ring_opened1D_damping_on_full_vTheta;
	//	below true:			Kinetic energy is computed on vTheta = v - ax, and vr
	//						kinetic energy of vr =  ring_opened1D_kinetic_energy_vr = rhoAverage * a^2 * L^2/8/PI^2 is added to any times kinetic energy
	//						and kinetic energy of vTheta = v - ax is considered
	bool b_ring_opened1D_kinetic_energy_on_full_vTheta;
	double ring_opened1D_kinetic_energy_vr;
	// = aL
	double ring_opened1D_al;

private:
	//////////////////////////////////////////
	// Initialize
	// input serial number < 0, means that the run does not have a serial number
	void Read_BaseData(istream& in, int serialNumberIn = -1);
	void PrepareForImpactIncidentEtcLoading();
	void FinalizeImpactIncidentEtcLoading();
	void Form_subdomains();
	bool Is_ActiveSubdomain(int sii, Subdomain_ElasticFractureModifier*& subdomain, unsigned int& si, unsigned int num_subdomains_tmp);
	void Form_Bulks_Interfaces_WithoutFormingConnections();
	void Connect_Interfaces_Set_BC_Types__Form_PP();
	// subdomainNo4AllBulks[bulki] is the subdomain number for bulks[bulksi]
	void Generate_subdomain_nos_for_all_bulks(vector<int>& subdomainNo4AllBulks);

	//////////////////////////////////////////
	// Initial condition
	void Set_InitialCondition_step();
	// time step until reaching the terminal condition
	int TimeStepsNonAdaptive();
	// open raw data files [scalarValues and finalSolutions] all interfaces (spatial points) for one time 
	// returns true if print should happen
	bool OpenFiles_RawData_OneTimeAllSpatialPoints(unsigned int timeIndex, double time, unsigned int numTimes);
	void Print__RawData_OneTimeAllSpatialPoints(SL_OneInterfaceAllTimes* interfacePtr, SLInterfaceCalculator& slic);
	void Close_Files_RawData_OntTimeAllSpatialPoints();

	vector<ostream*> outScalars_fixed_time, outFinalSln_fixed_time;

	bool b_x_min,b_subdomain_one_part_size, b_serialNumber, b_ring_R;
	unsigned int interface_offset;
	string configName;

	unsigned int GetInterfaceBulkSide_Subdomains_RelIndices(SL_OneInterfaceAllTimes* interfacePtr, vector<unsigned int>& subDomainNos, vector<unsigned int>& relPos_wrt_subDomainStartPoints);

	// to compute 1D averages of E, rho, ...
	void Compute1D_Averages();

	////////////////////////////////////////////////////////////////
	/// 1D visualization files generated
	bool b_visualization1D;
	unsigned int visualization1D_numFlds;
	// indexed by subdomain, then inteface number
	vector< vector<OneVisContour_xInfo> > visualization1D_xInfo;
	double visualization_TimeStep;
	void Delete_v1DFiles();
	void Size_v1DFiles();
	void Initialize_v1D();
	void Print_v1D(double time);
	
	vector<vector<ostream*> > v1DOutPtr;
	fstream v1Dtout; // output of time values
};


extern Domain_All_Interfaces_All_Times* g_domain;

// if serialNumberIn >= 0, it will be used in reading input material files

// There are TWO ways to read config files
//	A) configBC == "none"
//		config1: includes 3 separate config file names in it: 1. domain, 2. BC, 3. Time step and adaptivity
//  B) configBC != "none"
//		config1 -> Domain;	configBC -> BC;		config_TS_Adapt		-> Time step and adaptivity
int MAIN_Domain(string config1 = "config/Domain/MM8_Left_N_Right_N_Square_pulse_main.txt", int serialNumberIn = -1, string configBC = "none", string config_TS_Adapt = "none");

void Configure_sfcm_sfcm_gen();

#endif