#ifndef SUBDOMAIN_ELASTIC_FRACTURE_FACTORS__H
#define SUBDOMAIN_ELASTIC_FRACTURE_FACTORS__H

// this file is used to read modifications (in forms of factors) to elastic (bulk) and fracture (interface) properties
#include "InhomogeneousElasticFracture.h"

// currently there are two ways to read Subdomain_ElasticFractureModifier:
// 1. bisr_direct:	Direction:	the two vectors in "Subdomain_ElasticFractureModifier" are directly read from the text file
// 2. Inhomogeneous (or homogeneous) with one material flag. It reads it based on ElasticFractureInhomogField
typedef enum {bisr_undecided, bisr_direct, bisr_homogeneous, bisr_inhomogeneous, BulkInterfaceSeqReadT_SIZE} BulkInterfaceSeqReadT;

string getName(BulkInterfaceSeqReadT dat);
void name2Type(string& name, BulkInterfaceSeqReadT& typeVal);
ostream& operator<<(ostream& out, BulkInterfaceSeqReadT dat);
istream& operator>>(istream& in, BulkInterfaceSeqReadT& dat);


class Bulk_Elastic_Modifier
{
public:
	Bulk_Elastic_Modifier();
	// returns true if it reads one instance, false if not (i.e. reading })
	bool Read_Bulk_Elastic_Modifier(istream& in);
	void Print(ostream& out) const;

	// read data
	double length; // length of the segment
	GID	bulk_flag;
	// property factors
	double CFactor, rhoFactor, dampingFactor;

	// computed data
	// if any of the factors != 1, this modifier modifies reference elastic property
	bool b_modifies;

	//auxiliary functions
	void Initialize_Bulk_Elastic_Modifier();
};

class Interface_Fracture_Modifier
{
public:
	Interface_Fracture_Modifier();
	// returns true if it reads one instance, false if not (i.e. reading })
	bool Read_Interface_Fracture_Modifier(istream& in);
	void Print(ostream& out) const;
	// read data
	GID	interface_flag; // flag 0 means that the interface is bonded
	// property factors
	double sigmaFactor, deltaFactor, iniDamage;
};

// the numbering is as follows:
// |  bulk_0 |(interface_0)   bulk_1  |(interface_1) .....   bulk_(n-1) |(interface_(n-1)
// that is the very left interface is not given in intefaceMs below 
// the last interface is used to either end the domain OR connect to the next sequence
class Subdomain_ElasticFractureModifier
{
public:
	Subdomain_ElasticFractureModifier();
	bool Read_Subdomain_ElasticFractureModifier(unsigned int subdomainNumber, string configFileNameIn, int* serialNumberPtrIn = NULL, bool* isPeriodicPtrIn = NULL, double* xMPtrIn = NULL, double* xmPtrIn = NULL);
	bool Read_Subdomain_ElasticFractureModifier(istream& in, int* serialNumberPtrIn = NULL, bool* isPeriodicPtrIn = NULL, double* xMPtrIn = NULL, double* xmPtrIn = NULL, unsigned int subdomainNumber = 0);
	void Print(ostream& out) const;

	// the sequence below can be repeated numRepeatSequence (typically just 1)
	// read
	// 0, 1 -> boolean (-1: Domain decides on this)
	int isActive;
	bool isPeriodic;
	//// to form the vector of positions, interface_xs, the following 4 are used. Some or few of those are read as needed.
	// this can be read from the file and used for inhomogeneous data in which the number of segment from inhomogeneous data is multiplied by this number (if < 0). If > 0 and there is no vertices in homogeneous data this number is used to divide the domain.
	unsigned int num_bulk_one_seq;
	vector<double> interface_xs;
	// bound of the domain
	double xm, xM;
	///// end of data for forming interface_xs

	unsigned int numRepeatSequence;

	// direct read option bisrt == bisr_direct
	// if data is directly read, the file name w/o serial number is given below. From that file, bulk MS, and interfaceMs are read
	string directDataFileNameWOSerialExt;
	// directSpaceSizeModifier == 0 do nothing
	// > 0	-> bulk sizes should be smaller than this
	// < 0  -> refinement ratio (e.g. -2 -> each segment is divided into two)
	double directSpaceSizeModifier;
	// read from data file name
	vector<Bulk_Elastic_Modifier> bulkMs;
	vector<Interface_Fracture_Modifier> intefaceMs;

	// (in)homogeneous file read option bisrt == bisr_inhomogeneous
	string inhomogeneousConfigName;
	GID default_bulk_flag, default_interface_flag;
	
	inline double getTotalLength() const { return length_one_seq * numRepeatSequence; }
	inline unsigned int getTotalNumBulk() const { return num_bulk_one_seq * numRepeatSequence; }

	// computed
	int serialNumber;
	string serialNumber_str;
	double length_one_seq;
	BulkInterfaceSeqReadT bisrt;

private: // whether such data is read
	bool b_num_bulk_one_seq;
	bool b_xm, b_xM;
	bool b_serialNumber;
	bool b_isPeriodic;
	ElasticFractureInhomogField efif;
	string configFileName;
	unsigned long subdomainNumber;

	bool Read_BulkInterfaces(istream& in, double directSpaceSizeModifier);
	void Initialize_Subdomain_ElasticFractureModifier();
	void Initialize_Subdomain_ElasticFractureModifier_Direct();
	void Initialize_Subdomain_ElasticFractureModifier__In_homogeneous();
};

#endif