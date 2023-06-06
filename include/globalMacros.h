#ifndef FSI_GLOBAL__H
#define FSI_GLOBAL__H

#include "Dims.h"
#include "commonMacros.h"
#include "SimpleConfigMaker.h"

// CFL
// See the private class in SL_interfacePPtData_Time_Seq, the size - timeSeqPtrs
//#define SIZE_PT_TIME_SEQUENCE	7

//////////////////////#define SIZE_PT_TIME_SEQUENCE	10000
/// mainConfig_resolution_x_New.txt, configMaker_axt_resolution_x_Fracture_New.inst, resolutionFactors
// resolutionFactor * 4 + 2
// 1
//#define SIZE_PT_TIME_SEQUENCE	6
// 2
//#define SIZE_PT_TIME_SEQUENCE	10
// 4
//#define SIZE_PT_TIME_SEQUENCE	18
// 8
//#define SIZE_PT_TIME_SEQUENCE	34
// 16
#define SIZE_PT_TIME_SEQUENCE	66
// 32
//#define SIZE_PT_TIME_SEQUENCE	130
// 64
//#define SIZE_PT_TIME_SEQUENCE	258
// 128
//#define SIZE_PT_TIME_SEQUENCE	514
// 256
//#define SIZE_PT_TIME_SEQUENCE	1026
// 512
//#define SIZE_PT_TIME_SEQUENCE	2050
// 1024
//#define SIZE_PT_TIME_SEQUENCE	4098

//#define SIZE_PT_TIME_SEQUENCE	162
//#define SIZE_PT_TIME_SEQUENCE	2002
//#define SIZE_PT_TIME_SEQUENCE	50000

// prints velocity fields in DSU files (useful for elstic problems)
#define DSU_PRINT_VS	1
// Ring problem is 1D ring fragmentation problem. It has a source terms from v_r
#define RING_PROBLEM	0

#define USE_ISO_ASSUMPTION_PRE	0	// if iso assumption is used, normal and shear modes (2D, 3D) are decoupled and C, Z, Y, ..., in matrix form are not used
#define USE_PLANE_STRAIN	1	// for 2D isotropic model uses plane strain
// turn the flag below one, if the problem has source term that is 0th order in q = [v, sigma]
// For example, when damping is nonzero (D_vv = damping / rho) the flag below is 1
// see D matrix components in SL_Bulk_Properties
#define HAVE_SOURCE_ORDER0_q	0 //1
#define HAVE_SOURCE_TERM_OTHER	0	// source term other than 0th order q and ring problem source (e.g. constant body force ...)


#if HAVE_SOURCE_ORDER0_q
#define HAVE_SOURCE	1
#else
	#if RING_PROBLEM
	#define HAVE_SOURCE	1
	#else
		#if HAVE_SOURCE_TERM_OTHER
		#define HAVE_SOURCE	1
		#else
		#define HAVE_SOURCE	0
		#endif
	#endif
#endif

#if DiM1
#define USE_ISO_ASSUMPTION	1
#else
#define USE_ISO_ASSUMPTION	USE_ISO_ASSUMPTION_PRE
#endif

// C, Z, Y, etc. are only necessary for anisotropic C, but still can be computed for isotropic C
// for debugging, we can still compute them (second line), but first line should be the default
//#define COMPUTE_ANISO_BULK	!USE_ISO_ASSUMPTION
#define COMPUTE_ANISO_BULK	!USE_ISO_ASSUMPTION

#if USE_DEALII
#include <deal.II/base/tensor.h>
using namespace dealii ;
// uses deal.ii I/O format for vectors and matrices
// include the correct deal.ii files
typedef Tensor<1, DIM> VEC;
typedef Tensor<2, DIM> MAT;

#if DiM2a3_F
typedef Tensor<1, DIMm1> VECm1;
typedef Tensor<2, DIMm1> MATm1;
#endif

#else
//#include "LAfuncs.h"
//typedef VECTOR VEC;
//typedef MATRIX MAT;

#include "LAfuncsFixed.h"
typedef Vc_dd VEC;
typedef Mtrx_dd MAT;

#if 1 //DiM2a3_F
typedef Vc_dm1 VECm1;
typedef Mtrx_dm1 MATm1;
#endif

typedef long GID;
#endif
#endif

// don't use this
extern double s_invalidNum;

extern string g_prefileName;
extern string g_prefileNameWOSlash;
extern string g_prefileNamePPS2;

// these values are written to files when entry is invalid so that in read in it's clear the values are invalid
extern double invalidNonnegativeNum;
extern double invalidNum;
#define IS_INVALID(x) (fabs(x) > 1e39)

// serial number correspond to the input random domain to be read (dealing with input domain)
extern int g_serialNumber;
extern double g_time;

// if 1, it does not generate large files, 0 it does, 2 -> only prints fixed space data
extern int g_low_disk_space;

#define DB_ON 0
extern fstream dbout;
// can turn on or off printing for certain things
extern bool b_db_p;

// to print important things for the log of a particular run
extern fstream g_logout;
#if DB_ON
#define DBF(x) x
#define DBCHK(x) { if (b_db_p) { x } }
#else
#define DBF(x)
#define DBCHK(x)
#endif

/// this is a vert strict error check and even for things that can work but a bit iffy, the code exits
#define DB_STRICT_EXIT	0


// for incident, impact, Dirichlet/Neumann BC, ... cases there are potentially two ambient regions around a zone to be characterized
// they are in order outside/inside interfaces on the left, inside and outside interfaces on the right
#define LO_INT	0
#define LI_INT	1
#define RI_INT	2
#define RO_INT	3
extern int g_interfaceFlags_IncImp[4];
extern int g_bulkFlags_IncImp[NUM_SIDES];
// indices for printing for directions
extern vector<string> pINDs;
// this one has a comma too (if needed)
extern vector<string> pINDcommas;
bool Is_Valid_Num(double val);
void setGlobalMembers();
void CopyFile2OutputDirectory(const string& baseFileName);

