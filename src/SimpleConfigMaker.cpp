#include "SimpleConfigMaker.h"
#include "commonMacros.h"
#include "globalTypesClasses.h"
#include "globalFunctions.h"

SimpleFormatConfigMaker sfcm;
SimpleFormatConfigMaker_Generator sfcm_gen;

int g_versionNumber = -1;
int g_versionNumber_wOffset = -1;
string g_versionNumber_str = "";

#define DBBIN	0//1
#if DBBIN
fstream dbbin("dbbin.txt", ios::out);
#endif

SimpleFormatConfigMaker::SimpleFormatConfigMaker()
{
	success = false;
}

bool SimpleFormatConfigMaker::Read(string configMakerin, int cntrIn)
{
	Zero_SimpleFormatConfigMaker();
	configMaker = configMakerin;
	cntr = cntrIn;
	success = ReadSimpleFormatConfigMaker(configMaker, cntr, names, indices, sVals, addedName);
	return success;
}
int SimpleFormatConfigMaker::getIndexValueString(string groupNameIn, string& valOut)
{
	if (!success)
		return -1;
	int index = Find(names, groupNameIn);
	if (index >= 0)
		valOut = sVals[index];
	return index;
}
int SimpleFormatConfigMaker::getIndexValueDouble(string groupNameIn, double& valOut)
{
	if (!success)
		return -1;
	string valOutS;
	int index = getIndexValueString(groupNameIn, valOutS);
	if (index >= 0)
	{
		if (fromString(valOutS, valOut) == false)
			index = -1;
	}
	return index;
}

int SimpleFormatConfigMaker::getIndexValueInt(string groupNameIn, int& valOut)
{
	if (!success)
		return -1;
	string valOutS;
	int index = getIndexValueString(groupNameIn, valOutS);
	if (index >= 0)
	{
		if (fromString(valOutS, valOut) == false)
			index = -1;
	}
	return index;
}

void SimpleFormatConfigMaker::Zero_SimpleFormatConfigMaker()
{
	success = 0;
	cntr = 0;
	configMaker = "none";
	names.clear();
	indices.clear();
	sVals.clear();
	addedName = "";

	cfl_factor = -1.0;
	tFinal = -10.0;
	number_of_elements = -1024;
}

void generateArbitraryIndexListAux(const vector<int>& minVals, const vector<int>& maxVals, const vector<int>& incVals, vector< vector<int> >& resultsOut, vector< vector<int> >& resultsIn, int pos, bool inclusiveIntervals)
{
	int sz = resultsIn.size();
	if (sz == 0)
		THROW("sz == 0");
	int indsz = resultsIn[0].size();
	int maxV = maxVals[pos] + (int)inclusiveIntervals * incVals[pos];
	int numTotal = sz * ((maxV - minVals[pos]) / incVals[pos]);
	resultsOut.resize(numTotal);
	int cntr = 0;
	for (int i = minVals[pos]; i < maxV; i = i + incVals[pos])
	{
		for (int j = 0; j < sz; ++j)
		{
			resultsOut[cntr].resize(indsz + 1);
			resultsOut[cntr][0] = i;
			for (int k = 0; k < indsz; ++k)
				resultsOut[cntr][k + 1] = resultsIn[j][k];
			cntr = cntr + 1;
		}
	}
}

int generateArbitraryIndexList(const vector<int>& minVals, const vector<int>& maxVals, const vector<int>& incVals, vector< vector<int> >& resultsOut, bool inclusiveIntervals)
{
	resultsOut.resize(1);
	int sz = minVals.size();
	for (int pos = sz - 1; pos >= 0; --pos)
	{
		vector< vector<int> > resultsOutTemp;
		resultsOutTemp = resultsOut;
		generateArbitraryIndexListAux(minVals, maxVals, incVals, resultsOut, resultsOutTemp, pos, inclusiveIntervals);
	}
	return resultsOut.size();
}

void SimpleFormatConfigMaker_Generator::CreateSimpleFormatConfigMakerFromInstructions(string configMakerInstructionNameIn, bool forceRewrite)
{
	CreateSimpleFormatConfigMakerFromInstructions_Aux(configMakerInstructionNameIn);
	Write_SimpleFormatConfigMakerIterations(forceRewrite);
}

void SimpleFormatConfigMaker_Generator::Clear_SimpleFormatConfigMaker_Generator()
{
	subdomainNo = 0;
	bulkFlag = 1;
	interfaceFlag = 1;
	names.clear();
	s2sMap.clear();;
	svalVec.clear();;
	numGroups = 0;
	minVals.clear();
	maxVals.clear();
	incVals.clear();;
	fileNameOut = "";
	nameSpecific = "";
	specificProblemName = "";
	indicesOut.clear();
}

void SimpleFormatConfigMaker_Generator::CreateSimpleFormatConfigMakerFromInstructions_Aux(string configMakerInstructionNameIn)
{
	fstream in(configMakerInstructionNameIn.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "configMakerNameIn\t" << configMakerInstructionNameIn << "\ncannot be opened\n";
		THROW("Invalid file name\n");
	}
	Clear_SimpleFormatConfigMaker_Generator();
	string buf;
	READ_NSTRING(in, buf, buf);
	READ_NSTRING(in, buf, fileNameOut);
	if (fileNameOut == "default")
	{
		fileNameOut = removeExtension(configMakerInstructionNameIn);
		fileNameOut += ".txt";
	}
	READ_NSTRING(in, buf, buf);
	if (buf == "specificProblemName")
	{
		READ_NINTEGER(in, buf, specificProblemName);
		READ_NSTRING(in, buf, buf);
	}
	if (buf == "subdomainNo")
	{
		READ_NINTEGER(in, buf, subdomainNo);
		READ_NSTRING(in, buf, buf);
	}
	if (buf == "bulkFlag")
	{
		READ_NINTEGER(in, buf, bulkFlag);
		READ_NSTRING(in, buf, buf);
	}
	if (buf == "interfaceFlag")
	{
		READ_NINTEGER(in, buf, interfaceFlag);
		READ_NSTRING(in, buf, buf);
	}
	READ_NSTRING(in, buf, nameSpecific);

	READ_NSTRING(in, buf, buf);
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf != "{")
			THROW("istream should start with {");
		string name;
		READ_NSTRING(in, buf, name);
		names.push_back(name);
		READ_NSTRING(in, buf, buf);
		map<string, string> mp;
		if (buf == "map")
		{
			string skey, sval;
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
				THROW("istream should start with {");
			READ_NSTRING(in, buf, buf);
			while (buf != "}")
			{
				if (buf != "(")
				{
					cout << buf << '\n';
					THROW("map key should be precided by (\n");
				}
				READ_NSTRING(in, buf, skey);
				READ_NSTRING(in, buf, buf);
				if (buf != ",")
				{
					cout << buf << '\n';
					THROW("map key and value should be separated by ,\n");
				}
				READ_NSTRING(in, buf, sval);
				mp[skey] = sval;
				READ_NSTRING(in, buf, sval);
				if (buf != ")")
				{
					cout << buf << '\n';
					THROW("map value should be followed by )\n");
				}
				READ_NSTRING(in, buf, buf);
			}
			READ_NSTRING(in, buf, buf);
		}
		s2sMap.push_back(mp);
		vector<string> svals;
		while (buf != "}")
		{
			svals.push_back(buf);
			READ_NSTRING(in, buf, buf);
		}
		svalVec.push_back(svals);
		READ_NSTRING(in, buf, buf);
	}
	numGroups = names.size();
	READ_NSTRING(in, buf, buf);
	readPos_2outputPos.resize(numGroups);
	outputPos_2readPos.resize(numGroups);
	if ((!in.eof()) && (buf == "groupOutputOrder"))
	{
		fill(readPos_2outputPos.begin(), readPos_2outputPos.end(), -1);
		string name;
		READ_NSTRING(in, buf, buf);
		if (buf != "{")
			THROW("istream should start with {");
		int cntr = 0;
		while (buf != "}")
		{
			READ_NSTRING(in, buf, name);
			int pos = Find(names, name);
			if (pos >= 0)
			{
				readPos_2outputPos[pos] = cntr++;
			}
		}
		for (unsigned int j = 0; j < numGroups; ++j)
			if (readPos_2outputPos[j] == -1)
				readPos_2outputPos[j] = cntr++;
		for (unsigned int j = 0; j < numGroups; ++j)
			outputPos_2readPos[readPos_2outputPos[j]] = j;
	}
	else
	{
		for (unsigned int j = 0; j < numGroups; ++j)
		{
			readPos_2outputPos[j] = j;
			outputPos_2readPos[j] = j;
		}
	}
}

void SimpleFormatConfigMaker_Generator::Write_SimpleFormatConfigMakerIterations(bool forceRewrite)
{
	if (!forceRewrite)
	{
		fstream in(fileNameOut.c_str(), ios::in);
		if (in.is_open())
			return;
	}
	fstream out(fileNameOut.c_str(), ios::out);
	minVals.resize(numGroups), maxVals.resize(numGroups), incVals.resize(numGroups);
	for (unsigned int i = 0; i < numGroups; ++i)
	{
		minVals[i] = 0;
		maxVals[i] = svalVec[i].size();
		incVals[i] = 1;
	}
	bool inclusiveIntervals = false;
	int totSize = generateArbitraryIndexList(minVals, maxVals, incVals, indicesOut, inclusiveIntervals);
	for (int cntr = 0; cntr < totSize; ++cntr)
	{
		vector<int> indices;
		indices = indicesOut[cntr];
		vector<string> strs(numGroups);
		out << "cntrV\t" << cntr;
		out << "\tnameSpec\t" << nameSpecific;
		out << "\t{";
		for (unsigned int i = 0; i < numGroups; ++i)
		{
			strs[i] = svalVec[i][indices[i]];
			out << "\tgrp\t" << names[i] << "\tind\t" << indices[i] << "\tstr\t" << strs[i];
		}
		out << "\t}\n";
	}
}

bool ReadSimpleFormatConfigMaker(string configMakerIn, int cntr, vector<string>& names, vector<int>& indices, vector<string>& sVals, string& addedName)
{
	fstream in(configMakerIn.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "configMakerIn\t" << configMakerIn << '\n';
		THROW("Cannot open the file\n");
	}
	string buf = "none";
	READ_NSTRING(in, buf, buf);
	bool success = false;
	string nameSpec;
	while ((buf == "cntrV") && (!in.eof()))
	{
		int cntrTmp;
		READ_NINTEGER(in, buf, cntrTmp);

		READ_NSTRING(in, buf, buf);
		READ_NSTRING(in, buf, nameSpec);

		if (cntrTmp == cntr)
		{
			success = true;
			break;
		}
		getline(in, buf);
		READ_NSTRING(in, buf, buf);
	}
	if (nameSpec == "none")
		addedName = "";
	else
		addedName = nameSpec;
	if (success == false)
		return false;

	READ_NSTRING(in, buf, buf);
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf != "grp")
		{
			cout << "buf\t" << buf << '\n';
			THROW("Line must start with group\n");
		}
		READ_NSTRING(in, buf, buf);
		names.push_back(buf);
		addedName += "_";
		addedName += buf;
		READ_NSTRING(in, buf, buf);
		int ind;
		READ_NINTEGER(in, buf, ind);
		indices.push_back(ind);
		READ_NSTRING(in, buf, buf);
		READ_NSTRING(in, buf, buf);
		sVals.push_back(buf);
		string add_nm;
		//		if (fromString(buf, tmpi) || fromString(buf, tmpd))
		//		if (fromString(buf, tmpi))
		//		add_nm = "v" + buf;
		//		else
		{
			toString(ind, buf);
			add_nm = "i" + buf;
		}
		addedName += add_nm;
		READ_NSTRING(in, buf, buf);
	}
	return true;
}

bool Find_Version_Value(const string & key, double & value, map<string, string>*& mpPtr)
{
	if (!sfcm.success)
		return false;
	unsigned int pos;
	unsigned int sz = sfcm.names.size();
	for (pos = 0; pos < sz; ++pos)
		if (sfcm_gen.names[pos] == key)
			break;
	if (pos != sz)
	{
		if (fromString(sfcm.sVals[pos], value) == false)
			return false;
		mpPtr = &sfcm_gen.s2sMap[pos];
		return true;
	}
	// can add new features for sfcm later
	return false;
}

bool Find_Version_String(const string& key, string& value, map<string, string>*& mpPtr)
{
	if (!sfcm.success)
		return false;
	unsigned int pos;
	unsigned int sz = sfcm.names.size();
	for (pos = 0; pos < sz; ++pos)
		if (sfcm_gen.names[pos] == key)
			break;
	if (pos != sz)
	{
		value = sfcm.sVals[pos];
		mpPtr = &sfcm_gen.s2sMap[pos];
		return true;
	}
	// can add new features for sfcm later
	return false;
}

