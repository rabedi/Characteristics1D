#include "DomainPostProcessingS3.h"
#include "SimpleConfigMaker.h"
#include "SLDescriptorData.h"

/////////////////
string getName(pps3_om_T dat)
{
	if (dat == pps3_om_scalars_vectors)
		return "scalars_vectors";
	if (dat == pss3_om_scalars)
		return "scalars";
	if (dat == pps3_om_vectors)
		return "vectors";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, pps3_om_T& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= pps3_om_T_SIZE)
			THROW("too large of a number\n");
		typeVal = (pps3_om_T)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < pps3_om_T_SIZE; ++i)
	{
		typeVal = (pps3_om_T)i; // casting integer to pps3_om_T, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading pps3_om_T\n");
}

//operator for output
ostream& operator<<(ostream& out, pps3_om_T dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, pps3_om_T& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}


string getAddedName4OutputFile(pps3_om_T dat)
{
	if (dat == pps3_om_scalars_vectors)
		return "";
	if (dat == pss3_om_scalars)
		return "_scalars";
	if (dat == pps3_om_vectors)
		return "_vectors";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

/////////////////
string getName(PPS3_IOT dat)
{
	if (dat == pps3_i)
		return "input";
	if (dat == pps3_o)
		return "output";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, PPS3_IOT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= PPS3_IOT_SIZE)
			THROW("too large of a number\n");
		typeVal = (PPS3_IOT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < PPS3_IOT_SIZE; ++i)
	{
		typeVal = (PPS3_IOT)i; // casting integer to PPS3_IOT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading PPS3_IOT\n");
}

//operator for output
ostream& operator<<(ostream& out, PPS3_IOT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, PPS3_IOT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getName(PPS3_TimeOT dat)
{
	if (dat == pps3_stageStat)
		return "pps3_stageStat";
	if (dat == pps3_timeStep)
		return "pps3_timeStep";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, PPS3_TimeOT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= PPS3_TimeOT_SIZE)
			THROW("too large of a number\n");
		typeVal = (PPS3_TimeOT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < PPS3_TimeOT_SIZE; ++i)
	{
		typeVal = (PPS3_TimeOT)i; // casting integer to PPS3_TimeOT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading PPS3_TimeOT\n");
}

//operator for output
ostream& operator<<(ostream& out, PPS3_TimeOT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, PPS3_TimeOT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}


DomainPostProcessS3::DomainPostProcessS3()
{
	out_baseName = "runSummaries";
	out_ext = "csv";
	spatialFieldResolutionCorrector = 0;
	temporalFieldTimeStepValOrNumSteps = -400;
	addInvalidData = false;
	timeStep4_DSU_outputs = 1;
//	sep = "\t,\t";
//	sep = "\t";
	sep = ",";

	version_seperatePP3Folders = true;
	print_uncomputed_vals_as_nan = true;
	version_print_version_columns = 1;
	version_print_before_accept_serial = false;
	// prints version Number from the list of versions in generator file
	version_print_version_No = true;
	version_print_version_No_wOffset = true;
	version_print_indices = true;
	version_print_values = true;

	output_modes.clear();
	output_modes.push_back(pss3_om_scalars);
	output_modes.push_back(pps3_om_scalars_vectors);
	output_modes.push_back(pps3_om_vectors);

	for (unsigned int i = 0; i < ScalarFieldOutputT_SIZE; ++i)
	{
		scalarFieldAddOn[i] = getName((ScalarFieldOutputT)i);
		scalarFieldAddOn_Latex[i] = getLatexName((ScalarFieldOutputT)i);
	}
}

void DomainPostProcessS3::MAIN_DomainPostProcessS3(const string configFileName)
{
	fstream inconfig(configFileName.c_str(), ios::in);
	if (!inconfig.is_open())
	{
		cout << "configFileName\t" << configFileName << '\n';
		THROW("cannot open file\n");
	}
	CopyFile2OutputDirectory(configFileName);

	if (DomainPostProcessS3_Read_WO_Initialization(inconfig) == false)
		return;
	// ensuring configPPS2 reads data
	configPPS2.Main_DomainPostProcessS2();

	string out_baseNameBK = out_baseName;

	bool accebleOverAlVals;
	int time_outputIndex = -1;
	//! A. Printing stage/stat data
	PPS3_TimeOT tot = pps3_stageStat;
	root = "_PPS3";
	fileOperation(makeD, root);
	string rt2 = root + "/" + "timeVals";
	MakeDir(rt2);

	if (version_seperatePP3Folders)
		out_baseName += g_versionNumber_str;

	if (outputTypeActive[tot] == 1)
	{
		string fileNameLog = root + "/" + out_baseName + "_invalid_runs.txt";
		fstream outerr(fileNameLog.c_str(), ios::out | ios::app);
		bool runPrinted;
		for (unsigned int j = 0; j < output_modes.size(); ++j)
		{
			runPrinted = ComputePrint_Data(tot, accebleOverAlVals, time_outputIndex, output_modes[j]);
			if ((j == 0) && (!accebleOverAlVals))
				outerr << g_versionNumber << sep << g_versionNumber_wOffset << sep << g_serialNumber << sep << accebleOverAlVals << sep << runPrinted << '\n';
		}
	}
	tot = pps3_timeStep;
	out_baseName = out_baseNameBK;
	if (version_seperatePP3Folders)
	{
		root += g_versionNumber_str;
//		out_baseName += g_versionNumber_str;
	}
	fileOperation(makeD, root);

	if (outputTypeActive[tot] == 1)
	{
		int time_outputIndex_sz = configPPS2.onesubdomainPPS2[configPPS2.mainSubdomainNo].segmentInfo.maxIndex_Interface_DSU_Fragment_Print / timeStep4_DSU_outputs;
		bool notPrintField = false; // ((g_low_disk_space) && (time_outputIndex_sz > 8));

		string fileNameLog = root + "/" + out_baseName + "_timeIndex_invalid_runs.txt";
		for (time_outputIndex = 1; time_outputIndex <= time_outputIndex_sz; ++time_outputIndex)
		{
			fstream outerr(fileNameLog.c_str(), ios::out | ios::app);
			bool runPrinted;
			for (unsigned int j = 0; j < output_modes.size(); ++j)
			{
				if (notPrintField && (output_modes[j] != pss3_om_scalars))
					continue;
				runPrinted = ComputePrint_Data(tot, accebleOverAlVals, time_outputIndex, output_modes[j]);
				if ((j == 0) && (!accebleOverAlVals))
					outerr << g_versionNumber << sep << g_versionNumber_wOffset << sep << g_serialNumber << sep << accebleOverAlVals << sep << runPrinted << sep << time_outputIndex << '\n';
			}
		}
	}
}


bool DomainPostProcessS3::DomainPostProcessS3_Read_WO_Initialization(istream& in)
{
	string key;
	map<string, string>* mpPtr;
	double value = -1;

	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "output_modes")
		{
			READ_NSTRING(in, buf, buf); // should be {
			if (buf != "{")
				THROW("Should start with {\n");
			READ_NSTRING(in, buf, buf);
			output_modes.clear();
			while (buf != "}")
			{
				pps3_om_T omode;
				fromString(buf, omode);
				output_modes.push_back(omode);
				READ_NSTRING(in, buf, buf);
			}
		}
		else if (buf == "version_seperatePP3Folders")
		{
			READ_NBOOL(in, buf, version_seperatePP3Folders);
		}
		else if (buf == "print_uncomputed_vals_as_nan")
		{
			READ_NBOOL(in, buf, print_uncomputed_vals_as_nan);
		}
		else if (buf == "version_seperatePP3Folders")
		{
			READ_NBOOL(in, buf, version_seperatePP3Folders);
		}
		else if (buf == "version_print_version_columns")
		{
			READ_NINTEGER(in, buf, version_print_version_columns);
		}
		else if (buf == "version_print_before_accept_serial")
		{
			READ_NBOOL(in, buf, version_print_before_accept_serial);
		}
		else if (buf == "version_print_version_No")
		{
			READ_NBOOL(in, buf, version_print_version_No);
		}
		else if (buf == "version_print_version_No_wOffset")
		{
			READ_NBOOL(in, buf, version_print_version_No_wOffset);
		}
		else if (buf == "version_print_indices")
		{
			READ_NBOOL(in, buf, version_print_indices);
		}
		else if (buf == "version_print_values")
		{
			READ_NBOOL(in, buf, version_print_values);
		}
		else if (buf == "configPPS2")
		{
			string fileName5PPS2;
			READ_NSTRING(in, buf, fileName5PPS2);
			if (fileName5PPS2 == "inline")
				configPPS2.DomainPostProcessS2_Read_WO_Initialization(in);
			else
			{
				if (fileName5PPS2 == "default")
					fileName5PPS2 = "config/Domain/PPS2ConfigTest.txt";
				fstream inpps2(fileName5PPS2.c_str(), ios::in);
				if (!inpps2.is_open())
				{
					cout << "fileName5PPS2\t" << fileName5PPS2 << '\n';
					THROW("Cannot open file\n");
				}
				configPPS2.DomainPostProcessS2_Read_WO_Initialization(inpps2);
			}
		}
		else if (buf == "out_baseName")
		{
			READ_NSTRING(in, buf, out_baseName);
		}
		else if (buf == "out_ext")
		{
			READ_NSTRING(in, buf, out_ext);
		}
		else if (buf == "addInvalidData")
		{
			READ_NBOOL(in, buf, addInvalidData);
		}
		else if (buf == "spatialFieldResolutionCorrector")
		{
			READ_NINTEGER(in, buf, spatialFieldResolutionCorrector);
			key = "output_resolution";
			bool found = Find_Version_Value(key, value, mpPtr);;
			if (found)
				spatialFieldResolutionCorrector = (int)round(value);
		}
		else if (buf == "temporalFieldTimeStepValOrNumSteps")
		{
			READ_NDOUBLE(in, buf, temporalFieldTimeStepValOrNumSteps);
		}
		else if (buf == "timeStep4_DSU_outputs")
		{
			READ_NINTEGER(in, buf, timeStep4_DSU_outputs);
		}
		else if (buf == "outputTypeActive")
		{
			for (unsigned int tot = 0; tot < PPS3_TimeOT_SIZE; ++tot)
			{
				READ_NINTEGER(in, buf, outputTypeActive[tot]);
				if (outputTypeActive[tot] == -1)
				{
					if (g_low_disk_space != 0)
						outputTypeActive[tot] = 0;
					else
						outputTypeActive[tot] = 1;
				}
			}
		}
		else if (buf == "timeStampOverwriters")
		{
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("istream should start with {");
			}
			while (true)
			{
				string nm;
				READ_NSTRING(in, buf, nm);
				if (nm == "}")
					break;
				PPS2_dataPointer_StageOverwritter datPtrOverwriter;
				if (datPtrOverwriter.PPS2_dataPointer_StageOverwritter_Read(in) == false)
					break;
				timeStampOverwriters[nm] = datPtrOverwriter;
			}
		}
		else if (buf == "dataPointers")
		{
			READ_NSTRING(in, buf, buf);
			if (buf == "TimeOT")
				READ_NSTRING(in, buf, buf);
			PPS3_TimeOT tot;
			name2Type(buf, tot);

			READ_NSTRING(in, buf, buf);
			if (buf == "IOT")
				READ_NSTRING(in, buf, buf);
			PPS3_IOT iot;
			name2Type(buf, iot);
			vector<PPS2_dataPointer>* dataPointersPtr = &dataPointers[tot][iot];

			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("istream should start with {");
			}
			while (true)
			{
				PPS2_dataPointer datPtr;
				if (datPtr.PPS2_dataPointer_Read(in) == false)
					break;
				if (datPtr.overwriterName != "")
				{
					map<string, PPS2_dataPointer_StageOverwritter>::iterator it = timeStampOverwriters.find(datPtr.overwriterName), ite = timeStampOverwriters.end();
					if (it != ite)
						it->second.Overwrite_PPS2_dataPointer(datPtr);
				}
				if (datPtr.isActive)
					dataPointersPtr->push_back(datPtr);
			}
		}
		else
		{
			cout << "Invalid buf\t" << buf << '\n';
			THROW("exit\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	if (out_ext != "csv")
		sep = "\t";

	for (unsigned int tot = 0; tot < PPS3_TimeOT_SIZE; ++tot)
	{
		unsigned int sz = dataPointers[tot][pps3_i].size() + dataPointers[tot][pps3_o].size();
		if (sz == 0)
			outputTypeActive[tot] = 0;
	}
//	if (g_low_disk_space == 1)
//		outputTypeActive[pps3_timeStep] = 0;
	return true;
}

void DomainPostProcessS3::Get_Version_Header_Names(vector<string>& names, vector<string>& namesLatex)
{
	if (version_print_version_columns == 0)
		return;
	string preName = "inp_s_", tmp;
	if (version_print_version_No)
	{
		names.push_back("verNo");
		namesLatex.push_back("verNo");
	}
	if (version_print_version_No_wOffset)
	{
		names.push_back("verWOffsetNo");
		namesLatex.push_back("verWOffsetNo");
	}
	for (unsigned int gi = 0; gi < sfcm_gen.numGroups; ++gi)
	{
		int numVals = sfcm_gen.svalVec[gi].size();
		if ((numVals == 1) && (version_print_version_columns == 2))
			continue;
		unsigned int gii = sfcm_gen.outputPos_2readPos[gi];
		if (version_print_indices)
		{
			tmp = "ind" + sfcm_gen.names[gii];
			names.push_back(tmp);
			namesLatex.push_back(tmp);
		}
		if (version_print_values)
		{
			tmp = preName + sfcm_gen.names[gii];
			names.push_back(tmp);
			namesLatex.push_back(tmp);
		}			
	}
}

void DomainPostProcessS3::Print_Version_Values(ostream& out)
{
	if (version_print_version_columns == 0)
		return;

	extern SimpleFormatConfigMaker sfcm;
	out << setprecision(22);
	if (version_print_version_No)
		out << g_versionNumber << sep;
	if (version_print_version_No_wOffset)
		out << g_versionNumber_wOffset << sep;
	for (unsigned int gi = 0; gi < sfcm_gen.numGroups; ++gi)
	{
		int numVals = sfcm_gen.svalVec[gi].size();
		if ((numVals == 1) && (version_print_version_columns == 2))
			continue;
		unsigned int gii = sfcm_gen.outputPos_2readPos[gi];
		if (version_print_indices)
			out << sfcm.indices[gii] << sep;
		if (version_print_values)
		{
			double tmpv = 0.0;
			string str = sfcm.sVals[gii];
			if (!fromString(sfcm.sVals[gii], tmpv))
			{
				setStatOp_type typeVal;
				if (name2Type(str, typeVal) == true)
					tmpv = (double)typeVal;
				else
					tmpv = sfcm.indices[gii];
			}
			out << tmpv << sep;
		}
	}
}

bool DomainPostProcessS3::ComputePrint_Data(PPS3_TimeOT tot, bool& accebleOverAlVals, int time_outputIndex, pps3_om_T output_mode)
{
	string added_name = getAddedName4OutputFile(output_mode);
	string root2 = root;
	if (added_name != "")
	{
//		root2 = root + "/" + added_name;
		root2 = root + "/" + added_name;
		MakeDir(root2);
	}
	g_SL_desc_data.Read_tdLoad();
	string rt = root2 + "/";

	accebleOverAlVals = true;
	string fileName = rt + out_baseName;
	OneTimeValuePPS2Data modifiableOneTimeSlns;
	PPS2_TimeStamp timeStamp4FixedTime;
	vector<double> vals;
	vector<string> names, names_latex;
	vector<int> names_indices;
	string fileNameWOExt = fileName;
	if (tot == pps3_timeStep)
	{
		string ser;
		toString(time_outputIndex, ser);
		fileName += "_timeIndex_";
		fileName += ser;
		fileNameWOExt += "_timeIndex";

		timeStamp4FixedTime.timeStampType = lst_actualTime;
		timeStamp4FixedTime.actualTime_index = time_outputIndex * timeStep4_DSU_outputs;
		timeStamp4FixedTime.actualTime_indexType = ti_DSU_Fragment;
		timeStamp4FixedTime.actualTime_val = -1.0;
	}
	fileName += ".";
	fileName += out_ext;
	// check if the file outputFile_exists
	fstream in(fileName.c_str(), ios::in);
	bool outputFile_exists = in.is_open();
	if (outputFile_exists)
		in.close();

	// first taking care of input parameters
	for (int iiot = 0; iiot < PPS3_IOT_SIZE; ++iiot)
	{
		PPS3_IOT iot = (PPS3_IOT)iiot;
		vector<PPS2_dataPointer>* datPtrs = &dataPointers[tot][iot];
		unsigned int sz = datPtrs->size();
		bool acceptableVals;
		ScalarFieldOutputT sfot;

		double scalarVal;
		string preName = "inp_", nameBase, name;
		if (iot == pps3_o)
			preName = "out_";
		unsigned int fldSz;
		for (unsigned int i = 0; i < sz; ++i)
		{
			PPS2_dataPointer* datPtr = &(*datPtrs)[i];
			if (tot == pps3_timeStep)
				datPtr->timeStamp = timeStamp4FixedTime;
			vector<double> vecVal;
			acceptableVals = configPPS2.Get_Scalar_Or_Vector_Output(*datPtr, scalarVal, vecVal, sfot, modifiableOneTimeSlns, temporalFieldTimeStepValOrNumSteps, spatialFieldResolutionCorrector);
			if (!acceptableVals)
			{
				accebleOverAlVals = false;
				if (!addInvalidData)
					return false;
			}
			if (sfot == sfo_scalar)
			{
				if (output_mode != pps3_om_vectors)
				{
					vals.push_back(scalarVal);
					if (!outputFile_exists)
					{
						name = preName + scalarFieldAddOn[sfot] + datPtr->name_In_CSV_file;
						names.push_back(name);
//						string nl = "{" + datPtr->name_Latex + "}" + scalarFieldAddOn_Latex[sfot];
						string nl = datPtr->name_Latex;
						names_latex.push_back(nl);
						names_indices.push_back(-1);
					}
				}
			}
			else
			{
				if (output_mode != pss3_om_scalars)
				{
					fldSz = vecVal.size();
					if (fldSz == 0)
					{
						accebleOverAlVals = false;
						if (!addInvalidData)
							return false;
					}
					for (unsigned int j = 0; j < fldSz; ++j)
						vals.push_back(vecVal[j]);
					if (!outputFile_exists)
					{
						nameBase = preName + scalarFieldAddOn[sfot] + datPtr->name_In_CSV_file;
						for (unsigned int j = 0; j < fldSz; ++j)
						{
							string ser;
							toString(j, ser);
							name = nameBase + ser;
							names.push_back(name);
//							string nl = "{" + datPtr->name_Latex + "}" + scalarFieldAddOn_Latex[sfot];
							string nl = datPtr->name_Latex;
							names_latex.push_back(nl);
							names_indices.push_back(j);
						}
					}
				}
			}
		}
	}
	unsigned int szVals = vals.size();
	fstream out;
	if (!outputFile_exists)
	{
		out.open(fileName.c_str(), ios::out);
		vector<string> namesVersionAcceptable, namesLatexVersionAccepable;
		if (version_print_before_accept_serial)
		{
			Get_Version_Header_Names(namesVersionAcceptable, namesLatexVersionAccepable);
			namesVersionAcceptable.push_back("runNo");				namesVersionAcceptable.push_back("accepable");
			namesLatexVersionAccepable.push_back("runNo");			namesLatexVersionAccepable.push_back("accepableo");
		}
		else
		{
			namesVersionAcceptable.push_back("runNo");				namesVersionAcceptable.push_back("accepable");
			namesLatexVersionAccepable.push_back("runNo");			namesLatexVersionAccepable.push_back("accepableo");
			Get_Version_Header_Names(namesVersionAcceptable, namesLatexVersionAccepable);
		}
		unsigned int szVH = namesVersionAcceptable.size(), szVHm1 = szVH - 1;
		for (unsigned int i = 0; i < szVH; ++i)
		{
			out << namesVersionAcceptable[i];
			if (i < szVHm1)
				out << sep;
		}
		for (unsigned int i = 0; i < szVals; ++i)
			out << sep << names[i];
		out << "\n";
		string fileHeaders = fileNameWOExt + ".header";
		fstream outh(fileHeaders.c_str(), ios::out);

		for (unsigned int i = 0; i < szVH; ++i)			{			outh << namesVersionAcceptable[i];			if (i < szVHm1)				outh << '\t';		}
		for (unsigned int i = 0; i < szVals; ++i) { outh << names[i] << '\t'; }		outh << "\n";
		for (unsigned int i = 0; i < szVH; ++i) { outh << namesLatexVersionAccepable[i];			if (i < szVHm1)				outh << '\t'; }
		for (unsigned int i = 0; i < szVals; ++i) { outh << names_latex[i] << '\t'; }		outh << "\n";
		for (unsigned int i = 0; i < szVH; ++i) { outh << -1;			if (i < szVHm1)				outh << '\t'; }
		for (unsigned int i = 0; i < szVals; ++i) { outh << names_indices[i] << '\t'; }		outh << "\n";
	}
	else
	{
		out.open(fileName.c_str(), ios::out | ios::app);
	}
	out << setprecision(22);
	if (version_print_before_accept_serial)
	{
		Print_Version_Values(out);
		out << g_serialNumber << sep << accebleOverAlVals << sep;
	}
	else
	{
		out << g_serialNumber << sep << accebleOverAlVals << sep;
		Print_Version_Values(out);
	}
	int indLast = szVals - 1;
	out << setprecision(22);
	if (!print_uncomputed_vals_as_nan)
	{
		for (unsigned int i = 0; i < szVals; ++i)
		{
			out << vals[i];
			if (i != indLast)
				out << sep;
		}
	}
	else
	{
		double vl;
		for (unsigned int i = 0; i < szVals; ++i)
		{
			vl = vals[i];
			if (!IS_INVALID(vl))
				out << vl;
			else
				out << "nan"; // "NaN";
			if (i != indLast)
				out << sep;
		}
	}
	out << "\n";
	return true;
}

void MAIN_DomainPostProcessS3(const string configFileName)
{
	DomainPostProcessS3 dpps3;
	dpps3.MAIN_DomainPostProcessS3(configFileName);
}

