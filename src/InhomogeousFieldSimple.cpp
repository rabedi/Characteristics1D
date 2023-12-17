#include "InhomogeousFieldSimple.h"

string InhomogeneousSimple::root = "../InhomogeneousFiles/";

InhomogeneousSimple::InhomogeneousSimple()
{
	index_min_mappedValues_coarsened = -1;
	serialNo = 0;
	llc = -3.0;
	shape = 2; // triangular
	delta = 0.9;
	// 2, 4, ....
	resolutionFactor = 1;
	ssoFS = "min"; // min, mean_arithmetic, mean_harmonic, valStart
	baseResolution_p = 14;
	isUniform = ((fabs(shape - 0.0) < 1e-4) || (delta < 1e-3));
	inv_shape = 1.0;
	if (!isUniform)
		inv_shape = 1.0 / shape;
}

void InhomogeneousSimple::Compute_InhomogeneousSimple(int serialNoIn, double llcIn, double deltaIn, int shapeIn, double resolutionFactorIn, string ssofFSIn, int baseResolution_pIn)
{
	serialNo = serialNoIn;
	llc = llcIn;
	delta = deltaIn;
	shape = shapeIn;
	resolutionFactor = resolutionFactorIn;
	ssoFS = ssofFSIn;
	baseResolution_p = baseResolution_pIn;

	string s_serialNo;
	toStringTmp(serialNo, s_serialNo);
	isUniform = ((fabs(shape - 0.0) < 1e-4) || (delta < 1e-3));
	inv_shape = 1.0;
	if (!isUniform)
		inv_shape = 1.0 / shape;

	//! 1. reading the base file
	string fnAdded = "";
	if (baseResolution_pIn == 14)
	{
		string lc = "z";
		if (llc > -0.51)
			lc = "-0.5";
		else if (llc > -1.01)
			lc = "-1.0";
		else if (llc > -1.51)
			lc = "-1.5";
		else if (llc > -2.01)
			lc = "-2.0";
		else if (llc > -2.51)
			lc = "-2.5";
		else if (llc > -3.01)
			lc = "-3.0";
		else if (llc > -3.51)
			lc = "-3.5";
		else if (llc > -4.01)
			lc = "-4.0";
		fnAdded = "cl" + lc + "_np16385/initial_values_" + s_serialNo + ".txt";
	}
	else
	{
		THROW("Cannot open the file\n");
	}
	string fileName = InhomogeneousSimple::root + fnAdded;
	fstream in(fileName, ios::in);
	if (!in.is_open())
	{
		cout << "fileName\t" << fileName << '\n';
		THROW("File cannot be opened\n");
	}
	unsigned int sz;
	in >> sz;
	snValues.resize(sz);
	for (int i = 0; i < sz; ++i)
		in >> snValues[i];
	unsigned int numSegments = 0;
	if (sz > 0)
		numSegments = sz - 1;

	//! 2. map values
	mappedValues.resize(sz);
	for (int i = 0; i < sz; ++i)
		mappedValues[i] = GenSymTriangle_getInverseCDF(snValues[i]);

	//! 3 coarsening:
	if (resolutionFactor < 1.1)
	{
		mappedValues_coarsened = mappedValues;
		coarseningRat = 0;
	}
	else
	{
		coarseningRat = (int)round(log(resolutionFactor) / log(2));
		unsigned numSegmentsNew = numSegments / resolutionFactor;
		if (numSegmentsNew == 0)
			numSegmentsNew = 1;
		unsigned int numOldSegInNewSeg = numSegments / numSegmentsNew;
		vector<int> subsegmentSizes(numSegmentsNew), starPos(numSegmentsNew + 1);
		//vector<vector<double> > weights(numSegmentsNew);
		starPos[0] = 0;
		for (unsigned int i = 0; i < numSegmentsNew - 1; ++i)
		{
			subsegmentSizes[i] = numOldSegInNewSeg;
			starPos[i + 1] = starPos[i] + numOldSegInNewSeg;
		}
		starPos[numSegmentsNew] = numSegments;
		subsegmentSizes[numSegmentsNew - 1] = starPos[numSegmentsNew] - starPos[numSegmentsNew - 1];	//		numSegments - numOldSegInNewSeg * (numSegmentsNew - 1);
		numSegments = numSegmentsNew;
		int numVertices = numSegments + 1;
		mappedValues_coarsened.resize(numVertices);

		int indexMin;
		for (unsigned int i = 0; i < numSegments; ++i)
		{
			unsigned int st = starPos[i], en = starPos[i + 1], szz = en - st;
			vector<double> valsTmp(szz);
			for (unsigned int j = st; j < en; ++j)
				valsTmp[j - st] = mappedValues[j];
			mappedValues_coarsened[i] = getStatvalueSimple(valsTmp, ssoFS, indexMin);
		}
		//if (valsAtVertices == 1)
		{
			//if // (isPeriodic)
			mappedValues_coarsened[numSegments] = mappedValues_coarsened[0];
		}
	}
	sz_mappedValues_coarsened = mappedValues_coarsened.size();
	unsigned int num_interiorPts = sz_mappedValues_coarsened - 1;
	min_mappedValues_coarsened = 1e50;
	mean_mappedValues_coarsened = 0.0;
	index_min_mappedValues_coarsened = -1;
	double tmp;
	for (unsigned int i = 0; i < num_interiorPts; ++i)
	{
		tmp = mappedValues_coarsened[i];
		if (tmp < min_mappedValues_coarsened)
		{
			min_mappedValues_coarsened = tmp;
			index_min_mappedValues_coarsened = i;
		}
		mean_mappedValues_coarsened += tmp;
	}
	mean_mappedValues_coarsened /= (double)num_interiorPts;
}

double InhomogeneousSimple::GenSymTriangle_getInverseCDF(double sn)
{
	static double mean = 1.0;
	if (isUniform)
		return mean;
	double p = 0.5 * (1.0 + erf(sn / sqrt(2.0))); // sn_cdf
	if (p < 0.5)
		return mean * ((1.0 - delta) + delta * pow(2.0 * p, inv_shape));
	return mean * ((1.0 + delta) - delta * pow(2.0 * (1.0 - p), inv_shape));
}


double getStatvalueSimple(const vector<double>& vals, string ssoFS, int& index)
{
	unsigned int sz = vals.size();
	if (sz == 0)
		return 0.0;
	if (ssoFS == "min")
	{
		index = 0;
		double minV = vals[0];
		for (unsigned int i = 0; i < sz; ++i)
			if (vals[i] < minV)
			{
				minV = vals[i];
				index = i;
			}
			return minV;
	}
	if (ssoFS == "max")
	{
		int index = 0;
		double maxV = vals[0];
		for (unsigned int i = 0; i < sz; ++i)
			if (vals[i] > maxV)
			{
				maxV = vals[i];
				index = i;
			}
		return maxV;
	}
	if (ssoFS == "mean_arithmetic")
	{
		if (sz <= 1)
			return 0.0;
		double meanV = 0.0;
		for (unsigned int i = 0; i < sz; ++i)
			meanV += vals[i];
		meanV /= (double)sz;
		return meanV;
	}
	if (ssoFS == "mean_harmonic")
	{
		if (sz <= 1)
			return 0.0;
		double meanV = 0.0;
		for (unsigned int i = 0; i < sz; ++i)
			meanV += 1.0 / vals[i];
		meanV /= (double)(sz - 1.0);
		return 1.0 / meanV;
	}
	if (ssoFS == "mean_geometric")
	{
		double meanV = 1.0;
		for (unsigned int i = 0; i < sz; ++i)
			meanV *= vals[i];
		meanV = pow(meanV, 1.0 / (double)sz);
		return meanV;
	}
	if (ssoFS == "valStart")
		return vals[0];
	if (ssoFS == "valEnd")
		return vals[sz - 1];
	cout << ssoFS << '\n';
	THROW("Invalid sso\n");
}

void Test_InhomogeneousSimple()
{
	int serialNoIn = 0; // use 1000 - 1240
	double llcIn = -3.5; // -4.5 (white noise), -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5
	double deltaIn = 0.9;
	int shapeIn = 2; // 1, 1.5, 2, 3, 4
	double resolutionFactorIn = 64; // Use 1 in general
	string ssofFSIn = "min"; // mean_arithmetic, mean_harmonic, valStart
	int baseResolution_pIn = 14;

	InhomogeneousSimple is;
	is.Compute_InhomogeneousSimple(serialNoIn, llcIn, deltaIn, shapeIn, resolutionFactorIn, ssofFSIn, baseResolution_pIn);

	unsigned int sz;
	{
		fstream out("mesh_sn.txt", ios::out);
		sz = is.snValues.size();
		out << sz << '\n';
		for (int i = 0; i < sz; ++i)
			out << is.snValues[i] << '\n';
	}
	{
		fstream out("mesh_mapped.txt", ios::out);
		sz = is.mappedValues.size();
		out << sz << '\n';
		for (int i = 0; i < sz; ++i)
			out << is.mappedValues[i] << '\n';
	}
	{
		fstream out("mesh_mapped_coarsened.txt", ios::out);
		sz = is.mappedValues_coarsened.size();
		out << sz << '\n';
		for (int i = 0; i < sz; ++i)
			out << is.mappedValues_coarsened[i] << '\n';
	}
	cout << "sz_mappedValues_coarsened\t" << is.sz_mappedValues_coarsened << '\n';
	cout << "min_mappedValues_coarsened\t" << is.min_mappedValues_coarsened << '\n';
	cout << "mean_mappedValues_coarsened\t" << is.mean_mappedValues_coarsened << '\n';
	cout << "coarseningRat\t" << is.coarseningRat << '\n';
	cout << "index_min_mappedValues_coarsened\t" << is.index_min_mappedValues_coarsened << '\n';
}
