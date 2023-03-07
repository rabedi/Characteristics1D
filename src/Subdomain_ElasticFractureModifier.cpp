#include "Subdomain_ElasticFractureModifier.h"
#include "LAfuncsFinalStep.h"
#include "commonFunctions.h"
#include "SLDescriptorData.h"

string getName(BulkInterfaceSeqReadT dat)
{
	if (dat == bisr_undecided)
		return "undecided";
	if (dat == bisr_direct)
		return "direct";
	if (dat == bisr_homogeneous)
		return "homogeneous";
	if (dat == bisr_inhomogeneous)
		return "inhomogeneous";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, BulkInterfaceSeqReadT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= BulkInterfaceSeqReadT_SIZE)
			THROW("too large of a number\n");
		typeVal = (BulkInterfaceSeqReadT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < BulkInterfaceSeqReadT_SIZE; ++i)
	{
		typeVal = (BulkInterfaceSeqReadT)i; // casting integer to BulkInterfaceSeqReadT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading BulkInterfaceSeqReadT\n");
}

//operator for output
ostream& operator<<(ostream& out, BulkInterfaceSeqReadT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, BulkInterfaceSeqReadT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

Bulk_Elastic_Modifier::Bulk_Elastic_Modifier()
{
	length = 1.0;
	bulk_flag = 1;
	CFactor = 1.0, rhoFactor = 1.0, dampingFactor = 1.0;
	b_modifies = false;
}

bool Bulk_Elastic_Modifier::Read_Bulk_Elastic_Modifier(istream& in)
{
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
		if (buf == "l")
		{
			READ_NDOUBLE(in, buf, length);
		}
		else if (buf == "f")
		{
			READ_NINTEGER(in, buf, bulk_flag);
		}
		else if (buf == "C")
		{
			READ_NDOUBLE(in, buf, CFactor);
		}
		else if (buf == "r")
		{
			READ_NDOUBLE(in, buf, rhoFactor);
		}
		else if (buf == "d")
		{
			READ_NDOUBLE(in, buf, dampingFactor);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Initialize_Bulk_Elastic_Modifier();
	return true;
}

void Bulk_Elastic_Modifier::Print(ostream & out) const
{
	out << "bulk_flag\t" << bulk_flag << '\t';
	out << "length\t" << length << '\t';
	out << "CFactor\t" << CFactor << '\t';
	out << "rhoFactor\t" << rhoFactor << '\t';
	out << "dampingFactor\t" << dampingFactor << '\t';
	out << "b_modifies\t" << b_modifies << '\n';
}

void Bulk_Elastic_Modifier::Initialize_Bulk_Elastic_Modifier()
{
	b_modifies = false;
	static double tol = 1e-7;
	if (fabs(CFactor - 1.0) > tol)
	{
		b_modifies = true;
		return;
	}
	if (fabs(rhoFactor - 1.0) > tol)
	{
		b_modifies = true;
		return;
	}
	if (fabs(dampingFactor - 1.0) > tol)
	{
		dampingFactor = true;
		b_modifies = true;
		return;
	}
}

Interface_Fracture_Modifier::Interface_Fracture_Modifier()
{
	interface_flag = 1;
	sigmaFactor = 1.0;
	deltaFactor = 1.0;
	iniDamage = 0.0;
}

bool Interface_Fracture_Modifier::Read_Interface_Fracture_Modifier(istream& in)
{
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
		if (buf == "f")
		{
			READ_NINTEGER(in, buf, interface_flag);
		}
		else if (buf == "s")
		{
			READ_NDOUBLE(in, buf, sigmaFactor);
		}
		else if (buf == "d")
		{
			READ_NINTEGER(in, buf, deltaFactor);
		}
		else if (buf == "i")
		{
			READ_NDOUBLE(in, buf, iniDamage);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	return true;
}

void Interface_Fracture_Modifier::Print(ostream & out) const
{
	out << "interface_flag\t" << interface_flag << '\t';
	out << "sigmaFactor\t" << sigmaFactor << '\t';
	out << "deltaFactor\t" << deltaFactor << '\t';
	out << "iniDamage\t" << iniDamage << '\n';
}

Subdomain_ElasticFractureModifier::Subdomain_ElasticFractureModifier()
{
	numRepeatSequence = 1;
	length_one_seq = 0.0;
	num_bulk_one_seq = 0;

	b_num_bulk_one_seq = false;
	b_xm = false;
	b_xM = false;
	b_serialNumber = false;
	b_isPeriodic = false;
	xm = 0.0;
	xM = 0.0;
	num_bulk_one_seq = 0;
	directDataFileNameWOSerialExt = "none";
	inhomogeneousConfigName = "none";
	default_bulk_flag = 1;
	default_interface_flag = 0;

	serialNumber = -1;
	serialNumber_str = "";
	bisrt = bisr_undecided;
	isPeriodic = false;
	isActive = -1;
	directSpaceSizeModifier = 0.0;
}

bool Subdomain_ElasticFractureModifier::Read_Subdomain_ElasticFractureModifier(unsigned int subdomainNumberIn, string configFileNameIn, int* serialNumberPtrIn, bool* isPeriodicPtrIn, double* xMPtrIn, double* xmPtrIn)
{
	configFileName = configFileNameIn;
	subdomainNumber = subdomainNumberIn;
	fstream in(configFileName.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "configFileName\t" << configFileName << '\n';
		THROW("Cannot open file\n");
	}
	return Read_Subdomain_ElasticFractureModifier(in, serialNumberPtrIn, isPeriodicPtrIn, xMPtrIn, xmPtrIn, subdomainNumber);
}

bool Subdomain_ElasticFractureModifier::Read_Subdomain_ElasticFractureModifier(istream& in, int* serialNumberPtrIn, bool* isPeriodicPtrIn, double* xMPtrIn, double* xmPtrIn, unsigned int subdomainNumber)
{
	bool versionChange = ((sfcm.success) && (sfcm_gen.subdomainNo == subdomainNumber));
	string str_subdomainNumber;
	toString(subdomainNumber, str_subdomainNumber);
	double value;
	string key;	map<string, string>* mpPtr;

	b_serialNumber = false;
	if (serialNumberPtrIn != NULL)
	{
		serialNumber = *serialNumberPtrIn;
		b_serialNumber = true;
	}
	b_isPeriodic = false;
	if (isPeriodicPtrIn != NULL)
	{
		isPeriodic = *isPeriodicPtrIn;
		b_isPeriodic = true;
	}
	b_xM = false;
	if (xMPtrIn != NULL)
	{
		xM = *xMPtrIn;
		b_xM = true;
	}
	b_xm = false;
	if (xmPtrIn != NULL)
	{
		xm = *xmPtrIn;
		b_xm = true;
	}
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			string toCheck = "infile_subdomain";
			string subdomainNumber_str;
			toString(subdomainNumber, subdomainNumber_str);
			toCheck += subdomainNumber_str;
			while ((buf != toCheck) && (!in.eof()))
				READ_NSTRING(in, buf, buf);
			if (in.eof())
			{
				cout << "toCheck\t" << toCheck << '\n';
				THROW("Reached end of file looking for toCheck string\n");
			}
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
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "serialNumber")
		{
			int tmpi;
			READ_NINTEGER(in, buf, tmpi);
			if (!b_serialNumber)
			{
				serialNumber = tmpi;
				b_serialNumber = true;
			}
		}
		else if (buf == "isPeriodic")
		{
			bool tmpb;
			READ_NBOOL(in, buf, tmpb);
			if (!b_isPeriodic)
			{
				isPeriodic = tmpb;
				b_isPeriodic = true;
			}
		}
		else if (buf == "xM")
		{
			double tmpd;
			READ_NDOUBLE(in, buf, tmpd);
			if (!b_xM)
			{
				xM = tmpd;
				if (versionChange)
				{
					key = "L";
					if (Find_Version_Value(key, value, mpPtr) == true)
						xM = 0.5 * value;
					else
					{
						key = "lL";
						if (Find_Version_Value(key, value, mpPtr) == true)
							xM = 0.5 * pow(10.0, value);
					}
				}
				b_xM = true;
			}
		}
		else if (buf == "xm")
		{
			double tmpd;
			READ_NDOUBLE(in, buf, tmpd);
			if (!b_xm)
			{
				xm = tmpd;
				if (versionChange)
				{
					key = "L";
					if (Find_Version_Value(key, value, mpPtr) == true)
						xm = -0.5 * value;
					else
					{
						key = "lL";
						if (Find_Version_Value(key, value, mpPtr) == true)
							xm = -0.5 * pow(10.0, value);
					}
				}
				b_xm = true;
			}
		}
		else if (buf == "num_bulk_one_seq")
		{
			READ_NINTEGER(in, buf, num_bulk_one_seq);
			key = "num_bulk_one_seq";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
			{
				num_bulk_one_seq = (int)value;
			}
			else
			{
				key = "num_bulk_one_seq" + str_subdomainNumber;
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					num_bulk_one_seq = (int)value;
				}
				else if (g_SL_desc_data.tdLoadComputer != NULL)
				{
					if ((subdomainNumber == 0) || (subdomainNumber == 2))
					{
						key = "num_bulk_one_seqAmbient";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							num_bulk_one_seq = (int)value;
						}
					}
				}
			}
		}
		else if (buf == "interface_xs")
		{
			ReadVector(in, interface_xs);
		}
		else if (buf == "bisrt")
		{
			in >> bisrt;
		}
		else if (buf == "numRepeatSequence")
		{
			READ_NINTEGER(in, buf, numRepeatSequence);
			key = "numRepeatSequence";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
			{
				numRepeatSequence = (int)value;
			}
			else
			{
				key = "numRepeatSequence" + str_subdomainNumber;
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					numRepeatSequence = (int)value;
				}
				else if (g_SL_desc_data.tdLoadComputer != NULL)
				{
					if ((subdomainNumber == 0) || (subdomainNumber == 2))
					{
						key = "numRepeatSequenceAmbient";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							numRepeatSequence = (int)value;
						}
					}
					else if (subdomainNumber == 1)
					{
						key = "numRepeatSequenceInside";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							numRepeatSequence = (int)value;
						}
					}
				}
			}
		}
		else if (buf == "directDataFileNameWOSerialExt")
		{
			READ_NSTRING(in, buf, directDataFileNameWOSerialExt);
		}
		else if (buf == "directSpaceSizeModifier")
		{
			READ_NDOUBLE(in, buf, directSpaceSizeModifier);
			key = "directSpaceSizeModifier";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
			{
				directSpaceSizeModifier = value;
			}
			else
			{
				key = "directSpaceSizeModifier" + str_subdomainNumber;
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					directSpaceSizeModifier = value;
				}
				else if (g_SL_desc_data.tdLoadComputer != NULL)
				{
					if ((subdomainNumber == 0) || (subdomainNumber == 2))
					{
						key = "directSpaceSizeModifierAmbient";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							directSpaceSizeModifier = value;
						}
					}
					else if (subdomainNumber == 1)
					{
						key = "directSpaceSizeModifierInside";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							directSpaceSizeModifier = value;
						}
					}
				}
			}
		}
		else if (buf == "inhomogeneousConfigName")
		{
			READ_NSTRING(in, buf, inhomogeneousConfigName);
		}
		else if (buf == "default_bulk_flag")
		{
			READ_NINTEGER(in, buf, default_bulk_flag);
		}
		else if (buf == "default_interface_flag")
		{
			READ_NINTEGER(in, buf, default_interface_flag);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Initialize_Subdomain_ElasticFractureModifier();
	return true;
}

void Subdomain_ElasticFractureModifier::Print(ostream & out) const
{
	out << "isActive\t" << isActive << '\t';
	out << "isPeriodic\t" << isPeriodic << '\t';
	out << "configFileName\t" << configFileName << '\t';
	out << "subdomainNumber\t" << subdomainNumber << '\t';
	out << "serialNumber\t" << serialNumber << '\t';
	out << "serialNumber_str\t" << serialNumber_str << '\t';
	out << "directDataFileNameWOSerialExt\t" << directDataFileNameWOSerialExt << '\t';
	out << "inhomogeneousConfigName\t" << inhomogeneousConfigName << '\n';

	out << "b_num_bulk_one_seq\t" << b_num_bulk_one_seq << '\t';
	out << "b_xm\t" << b_xm << '\t';
	out << "b_xM\t" << b_xM << '\t';
	out << "b_serialNumber\t" << b_serialNumber << '\t';
	out << "b_isPeriodic\t" << b_isPeriodic << '\t';
	out << "efifNOTPRINTED\n";


	out << "bisrt\t" << bisrt << '\t';
	out << "xm\t" << xm << '\t';
	out << "xM\t" << xM << '\t';
	out << "directSpaceSizeModifier\t" << directSpaceSizeModifier << '\n';

	out << "numRepeatSequence\t" << numRepeatSequence << '\t';
	out << "length_one_seq\t" << length_one_seq << '\t';
	out << "num_bulk_one_seq\t" << num_bulk_one_seq << '\n';
	out << "interface_xs\n";
	WriteVectorDouble(out, interface_xs);
	out << "\n";

	out << "default_bulk_flag\t" << default_bulk_flag << '\n';
	out << "bulkMs\t" << bulkMs.size() << '\n';
	for (unsigned int i = 0; i < bulkMs.size(); ++i)
	{
		out << "bulkMs[" << i << "]\n";
		bulkMs[i].Print(out);
	}

	out << "default_interface_flag\t" << default_interface_flag << '\n';
	out << "intefaceMs\t" << intefaceMs.size() << '\n';
	for (unsigned int i = 0; i < intefaceMs.size(); ++i)
	{
		out << "intefaceMs[" << i << "]\n";
		intefaceMs[i].Print(out);
	}
}

bool Subdomain_ElasticFractureModifier::Read_BulkInterfaces(istream& in, double directSpaceSizeModifier)
{
	string buf;
	//	immediately after each bulk and will be used between refined bulks in case of refinement of the edge includesInterfaceInfosFromrefinement
	bool includesInterfaceInfosFromrefinement = false;
	READ_NSTRING(in, buf, buf);
	if (buf == "includesInterfaceInfosFromrefinement")
	{
		READ_NBOOL(in, buf, includesInterfaceInfosFromrefinement);
		READ_NSTRING(in, buf, buf);
	}
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			THROW("istream should start with {");
		}
	}
	bool min_bulk_size = (directSpaceSizeModifier > 1e-40);
	unsigned int num_refinement;
	bool b_num_refinement = (directSpaceSizeModifier < -1);
	bool checkSize = (min_bulk_size || b_num_refinement);
	if (b_num_refinement)
		num_refinement = MAX((int)round(-directSpaceSizeModifier), 1);

	Bulk_Elastic_Modifier bem, bemNew;

	bool cont_reading = bem.Read_Bulk_Elastic_Modifier(in);
	double current_length, new_length;
	while (cont_reading)
	{
		Interface_Fracture_Modifier ifm0;
		if (includesInterfaceInfosFromrefinement)
		{
			cont_reading = ifm0.Read_Interface_Fracture_Modifier(in);
			if (!cont_reading)
			{
				THROW("It's supposed to read interface between the same bulk but it failed\n");
			}
		}
		else
			ifm0.interface_flag = 0;
		if (!checkSize)
			bulkMs.push_back(bem);
		else
		{
			current_length = bem.length;
			if (min_bulk_size)
				num_refinement = (int)ceil(current_length / directSpaceSizeModifier);
			if (num_refinement > 1)
			{
				new_length = current_length / num_refinement;
				bemNew = bem;
				bemNew.length = new_length;
				for (unsigned int i = 0; i < (num_refinement - 1); ++i)
				{
					bulkMs.push_back(bemNew);
					intefaceMs.push_back(ifm0);
				}
				bulkMs.push_back(bemNew);
			}
			else
				bulkMs.push_back(bem);
		}
		Interface_Fracture_Modifier ifm;
		cont_reading = ifm.Read_Interface_Fracture_Modifier(in);
		if (cont_reading)
		{
			intefaceMs.push_back(ifm);
			cont_reading = bem.Read_Bulk_Elastic_Modifier(in);
		}
	}
	return true;
}

void Subdomain_ElasticFractureModifier::Initialize_Subdomain_ElasticFractureModifier()
{
	if (b_serialNumber && (serialNumber >= 0))
	{
		toString(serialNumber, serialNumber_str);
		serialNumber_str = "_" + serialNumber_str;
	}
	if (interface_xs.size() > 1)
	{
		xm = interface_xs[0];
		xM = interface_xs[interface_xs.size() - 1];
		b_xm = true;
		b_xM = true;
	}
	if (bisrt == bisr_undecided)
	{
		if (directDataFileNameWOSerialExt != "none")
			bisrt = bisr_direct;
		else if (inhomogeneousConfigName != "none")
			bisrt = bisr_inhomogeneous;
		else 
			bisrt = bisr_homogeneous;
	}
	if (bisrt == bisr_direct)
		Initialize_Subdomain_ElasticFractureModifier_Direct();
	else if ((bisrt == bisr_inhomogeneous) || (bisrt == bisr_homogeneous))
		Initialize_Subdomain_ElasticFractureModifier__In_homogeneous();
}

void Subdomain_ElasticFractureModifier::Initialize_Subdomain_ElasticFractureModifier_Direct()
{
	// reading data
	string dataFileName;
	bool read2Start = false;
	if (directDataFileNameWOSerialExt == "infile")
	{
		dataFileName = configFileName;
		read2Start = true;
	}
	else
		dataFileName = directDataFileNameWOSerialExt + serialNumber_str + ".txt";
	fstream inData(dataFileName.c_str(), ios::in);
	if (!inData.is_open())
	{
		cout << "dataFileName\n" << dataFileName << '\n';
		THROW("cannot open file\n");
	}
	if (read2Start)
	{
		string buf = "";
		string str_subdomainNumber;
		toString(subdomainNumber, str_subdomainNumber);
		string check = "direction_subdomain_" + str_subdomainNumber;
		while (buf != check)
			READ_NSTRING(inData, buf, buf);
	}
	Read_BulkInterfaces(inData, directSpaceSizeModifier);

	num_bulk_one_seq = bulkMs.size();
	unsigned int interfaceSz = intefaceMs.size();
	if (interfaceSz == (num_bulk_one_seq - 1))
	{
		// adding a bonded interface at the end
		Interface_Fracture_Modifier intefaceM;
		intefaceMs.push_back(intefaceM);
	}
	else if (interfaceSz != num_bulk_one_seq)
	{
		cout << "interfaceSz\t" << interfaceSz << '\n';
		cout << "num_bulk_one_seq\t" << num_bulk_one_seq << '\n';
		THROW("Invalid size\t");
	}
	unsigned int sz_num_xs = interface_xs.size();
	if (sz_num_xs > 1)
	{
		if (sz_num_xs != (num_bulk_one_seq + 1))
		{
			THROW("(sz_num_xs != (num_bulk_one_seq + 1))\n");
		}
		length_one_seq = interface_xs[sz_num_xs - 1] - interface_xs[0];
		for (unsigned int i = 0; i < num_bulk_one_seq; ++i)
			bulkMs[i].length = interface_xs[i + 1] - interface_xs[i];
	}
	else
	{
		length_one_seq = 0.0;
		for (unsigned int i = 0; i < num_bulk_one_seq; ++i)
			length_one_seq += bulkMs[i].length;
		if (length_one_seq < 1e-7)
		{
			length_one_seq = xM - xm;
			double del = length_one_seq / num_bulk_one_seq;
			for (unsigned int i = 0; i < num_bulk_one_seq; ++i)
				bulkMs[i].length = del;
		}
	}
}

void Subdomain_ElasticFractureModifier::Initialize_Subdomain_ElasticFractureModifier__In_homogeneous()
{
	if (bisrt == bisr_inhomogeneous)
	{
		int *serialNumberPtr = NULL;
		if (b_serialNumber && (serialNumber >= 0))
			serialNumberPtr = &serialNumber;
		bool *isPeriodicPtr = NULL;
		if (b_isPeriodic)
			isPeriodicPtr = &isPeriodic;
		double* xMPtrIn = NULL; double* xmPtrIn = NULL;
		if (b_xM)
			xMPtrIn = &xM;
		if (b_xm)
			xmPtrIn = &xm;

		if (inhomogeneousConfigName == "infile")
			inhomogeneousConfigName = configFileName;
		efif.Read_ElasticFractureInhomogField(subdomainNumber, inhomogeneousConfigName, serialNumberPtr, isPeriodicPtr, xMPtrIn, xmPtrIn);
	}
	unsigned int sz_num_xs = interface_xs.size();
	if ((bisrt == bisr_homogeneous) || (efif.numVertices < 2)) // efif has no data!
	{
		bisrt = bisr_homogeneous;
		// field is homogeneous and basically need to make up the list of bulk and interfaces
		if (sz_num_xs < 2)
		{
			if ((b_xM == false) || (b_xm == false) || (num_bulk_one_seq <= 0))
			{
				THROW("no way to set interface_xs\n");
			}
			sz_num_xs = num_bulk_one_seq + 1;
			interface_xs.resize(sz_num_xs);
			length_one_seq = (xM - xm);
			double del = length_one_seq / num_bulk_one_seq;
			interface_xs[0] = xm;
			for (unsigned int i = 1; i <= num_bulk_one_seq; ++i)
				interface_xs[i] = interface_xs[i - 1] + del;
		}
		else
		{
			num_bulk_one_seq = sz_num_xs - 1;
			xm = interface_xs[0];
			xM = interface_xs[sz_num_xs - 1];
			b_xm = true;
			b_xM = true;
		}
		bulkMs.resize(num_bulk_one_seq);
		intefaceMs.resize(num_bulk_one_seq);
		for (unsigned int i = 0; i < num_bulk_one_seq; ++i)
		{
			Bulk_Elastic_Modifier* pemPtr = &bulkMs[i];
			pemPtr->bulk_flag = default_bulk_flag;
			pemPtr->length = interface_xs[i + 1] - interface_xs[i];
			pemPtr->b_modifies = false;
			pemPtr->CFactor = 1.0;
			pemPtr->rhoFactor = 1.0;
			pemPtr->dampingFactor = 1.0;

			Interface_Fracture_Modifier* ifmPtr = &intefaceMs[i];
			ifmPtr->interface_flag = default_interface_flag;
			ifmPtr->iniDamage = 0.0;
			ifmPtr->sigmaFactor = 1.0;
			ifmPtr->deltaFactor = 1.0;
		}
	}
	else
	{
		if (sz_num_xs < 2)
		{
			if (b_num_bulk_one_seq)
			{
				double xmDat = efif.xs[0], xMDat = efif.xs[efif.xs.size() - 1];
				if (b_xm)
				{
					if (xm > xmDat)
						xmDat = xm;
				}
				if (b_xM)
				{
					if (xM < xMDat)
						xMDat = xM;
				}
				if (num_bulk_one_seq < 0)
					num_bulk_one_seq *= -efif.numSegments;
				double del = (xMDat - xmDat) / num_bulk_one_seq;
				sz_num_xs = num_bulk_one_seq + 1;
				interface_xs.resize(sz_num_xs);
				for (unsigned int i = 0; i < sz_num_xs; ++i)
					interface_xs[i] = xmDat + i * del;
			}
			else
			{
				num_bulk_one_seq = efif.numSegments;
				sz_num_xs = num_bulk_one_seq + 1;
				interface_xs = efif.xs;
			}
		}
		num_bulk_one_seq = sz_num_xs - 1;
		bulkMs.resize(num_bulk_one_seq);
		intefaceMs.resize(num_bulk_one_seq);
		for (unsigned int i = 1; i <= num_bulk_one_seq; ++i)
		{
			double xInterface = interface_xs[i];
			ElasticFractureInhomogFactors efifrac;
			efif.getFactorsIniDamage_By_x(xInterface, efifrac);

			Interface_Fracture_Modifier* ifmPtr = &intefaceMs[i - 1];
			ifmPtr->interface_flag = default_interface_flag;
			ifmPtr->iniDamage = efifrac.iniDamage;
			ifmPtr->sigmaFactor = efifrac.sigmaFactor;
			ifmPtr->deltaFactor = efifrac.deltaFactor;

			// bulk
			double xBulk = 0.5 * (interface_xs[i] + interface_xs[i - 1]);
			ElasticFractureInhomogFactors efibulk;
			efif.getFactorsIniDamage_By_x(xBulk, efibulk);

			Bulk_Elastic_Modifier* pemPtr = &bulkMs[i - 1];
			pemPtr->bulk_flag = default_bulk_flag;
			pemPtr->length = interface_xs[i] - interface_xs[i - 1];
			pemPtr->CFactor = efibulk.CFactor;
			pemPtr->rhoFactor = efibulk.rhoFactor;
			pemPtr->dampingFactor = efibulk.dampingFactor;
			pemPtr->Initialize_Bulk_Elastic_Modifier();
		}
	}
}
