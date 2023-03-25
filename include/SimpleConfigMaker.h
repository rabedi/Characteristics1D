#ifndef SIMPLE_CONFIG_MAKER__H
#define SIMPLE_CONFIG_MAKER__H

#include "Dims.h"

class SimpleFormatConfigMaker
{
public:
	SimpleFormatConfigMaker();
	//	returns true if successfull read
	bool Read(string configMakerin, int cntrIn);
	int getIndexValueString(string groupNameIn, string& valOut);
	int getIndexValueDouble(string groupNameIn, double& valOut);
	int getIndexValueInt(string groupNameIn, int& valOut);
	bool success;
	int cntr;
	string configMaker;
	vector<string> names;
	vector<int> indices;
	vector<string> sVals;
	string addedName;
	    
	// added for characteristic code
	double cfl_factor;
	double tFinal; 
	int number_of_elements; // how many elements the input mesh should have - this is number of sampling for KL mesh
	int direct_resolutionFactor;
private:
	void Zero_SimpleFormatConfigMaker();
};

void generateArbitraryIndexListAux(const vector<int>& minVals, const vector<int>& maxVals, const vector<int>& incVals, vector< vector<int> >& resultsOut, vector< vector<int> >& resultsIn, int pos, bool inclusiveIntervals);
int generateArbitraryIndexList(const vector<int>& minVals, const vector<int>& maxVals, const vector<int>& incVals, vector< vector<int> >& resultsOut, bool inclusiveIntervals);

class SimpleFormatConfigMaker_Generator
{
public:
	SimpleFormatConfigMaker_Generator();
	void CreateSimpleFormatConfigMakerFromInstructions(string configMakerInstructionNameIn, bool forceRewrite = false);
	vector<string> names;	
	vector<map<string, string> > s2sMap;
	vector<vector<string> > svalVec;
	unsigned int numGroups;
	vector<int> minVals, maxVals, incVals;
	string fileNameOut, nameSpecific;
	vector< vector<int> > indicesOut;
	// for many runs we deal with specific subdomains that need to be changed, bulk and inteface properties to be changed
	// the positions of these are given below
	int subdomainNo;
	int bulkFlag;
	int interfaceFlag;
	// this string is used to set a few specific things for different problems, e.g. axt problem
	string specificProblemName;
	// We want to output fields in summary files (e.g. PPS3 of characteristics code) in groupOutputOrder order. 
	// This may not match the order in names. The two maps below go 
	// from output position to those read in order in names
	// and opposite direction
	vector<int> outputPos_2readPos, readPos_2outputPos;
	vector<string> groupOutputOrder;
private:
	void Clear_SimpleFormatConfigMaker_Generator();
	// simple functions for creating list of options
	void CreateSimpleFormatConfigMakerFromInstructions_Aux(string configMakerInstructionNameIn);
	// forceRewrite: if true, always will write the sfcm iterations
	void Write_SimpleFormatConfigMakerIterations(bool forceRewrite = false);
};
// bool is for success
bool ReadSimpleFormatConfigMaker(string configMakerIn, int cntr, vector<string>& names, vector<int>& indices, vector<string>& sVals, string& addedName);

extern SimpleFormatConfigMaker sfcm;
extern SimpleFormatConfigMaker_Generator sfcm_gen;

// version number corresponds to run parameters other than input random field number (which is handled by g_serialNumber)
// examples are: looping over deltaC, correlation length and range of values for random field, loading rate etc.
// this is taken care of by sfcm and sfcm_gen
extern int g_versionNumber;
extern int g_versionNumber_wOffset;
extern string g_versionNumber_str;

// returns if the value is stored in sfcm and it has a double value
// value is the coresponding value stored in sfcm
// mpPtr is extra map for this key
bool Find_Version_Value(const string& key, double& value, map<string, string>*& mpPtr);
bool Find_Version_String(const string& key, string& value, map<string, string>*& mpPtr);

#endif