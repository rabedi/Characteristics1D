#include "globalFunctions.h"
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
using namespace std;

int SetNewtonCotes_Points_AndWeights(int numSpatialSubsegments_BulkInterfacePoints_Print, bool useRepeatedSimpsonRuleForHigherOrders, vector<double>& spatialIntegrationWeights, vector<double>& spatialIntegrationPoints)
{
	int numSpatialPointsPerSegment = numSpatialSubsegments_BulkInterfacePoints_Print + 1;
	spatialIntegrationWeights.resize(numSpatialPointsPerSegment);

	spatialIntegrationPoints.resize(numSpatialPointsPerSegment);
	double h = 1.0 / numSpatialSubsegments_BulkInterfacePoints_Print;
	for (int i = 0; i < numSpatialPointsPerSegment; ++i)
		spatialIntegrationPoints[i] = i * h;

	if (numSpatialPointsPerSegment == 2)
	{
		spatialIntegrationWeights[0] = 0.5;	spatialIntegrationWeights[1] = 0.5;
	}
	else if (numSpatialPointsPerSegment == 3)
	{
		spatialIntegrationWeights[0] = 1.0 / 6.0;	spatialIntegrationWeights[1] = 4.0 / 6.0;	spatialIntegrationWeights[2] = 1.0 / 6.0;
	}
	else
	{
		if (!useRepeatedSimpsonRuleForHigherOrders && (numSpatialPointsPerSegment <= 11))
		{
			if (numSpatialPointsPerSegment == 4)
			{
				double tmp = 1.0 / 8.0;
				spatialIntegrationWeights[0] = tmp;	spatialIntegrationWeights[1] = 3.0 * tmp;	spatialIntegrationWeights[2] = 3.0 * tmp;	spatialIntegrationWeights[3] = tmp;
			}
			else if (numSpatialPointsPerSegment == 5)
			{
				double tmp = 1.0 / 90.0;
				spatialIntegrationWeights[0] = 7.0 * tmp;		spatialIntegrationWeights[4] = 7.0 * tmp;
				spatialIntegrationWeights[1] = 32.0 * tmp;		spatialIntegrationWeights[3] = 32.0 * tmp;
				spatialIntegrationWeights[2] = 12.0 * tmp;
			}
			else if (numSpatialPointsPerSegment == 6)
			{
				double tmp = 1.0 / 288.0;
				spatialIntegrationWeights[0] = 19.0 * tmp;		spatialIntegrationWeights[5] = 19.0 * tmp;
				spatialIntegrationWeights[1] = 75.0 * tmp;		spatialIntegrationWeights[4] = 75.0 * tmp;
				spatialIntegrationWeights[2] = 50.0 * tmp;		spatialIntegrationWeights[3] = 50.0 * tmp;
			}
			else if (numSpatialPointsPerSegment == 7)
			{
				double tmp = 1.0 / 840.0;
				spatialIntegrationWeights[0] = 41.0 * tmp;		spatialIntegrationWeights[6] = 41.0 * tmp;
				spatialIntegrationWeights[1] = 216.0 * tmp;		spatialIntegrationWeights[5] = 216.0 * tmp;
				spatialIntegrationWeights[2] = 27.0 * tmp;		spatialIntegrationWeights[4] = 27.0 * tmp;
				spatialIntegrationWeights[3] = 272.0 * tmp;
			}
			else if (numSpatialPointsPerSegment == 8)
			{
				double tmp = 7.0 / 17280.0 / 7.0;
				spatialIntegrationWeights[0] = 751.0 * tmp;		spatialIntegrationWeights[7] = 751.0 * tmp;
				spatialIntegrationWeights[1] = 3577.0 * tmp;		spatialIntegrationWeights[6] = 3577.0 * tmp;
				spatialIntegrationWeights[2] = 1323.0 * tmp;		spatialIntegrationWeights[5] = 1323.0 * tmp;
				spatialIntegrationWeights[3] = 2989.0 * tmp;		spatialIntegrationWeights[4] = 2989.0 * tmp;
			}
			else if (numSpatialPointsPerSegment == 9)
			{
				double tmp = 4.0 / 14175.0 / 8.0;
				spatialIntegrationWeights[0] = 989.0 * tmp;		spatialIntegrationWeights[8] = 989.0 * tmp;
				spatialIntegrationWeights[1] = 5888.0 * tmp;		spatialIntegrationWeights[7] = 5888.0 * tmp;
				spatialIntegrationWeights[2] = -928.0 * tmp;		spatialIntegrationWeights[6] = -928.0 * tmp;
				spatialIntegrationWeights[3] = 10496.0 * tmp;		spatialIntegrationWeights[5] = 10496.0 * tmp;
				spatialIntegrationWeights[4] = -4540.0 * tmp;
			}
			else if (numSpatialPointsPerSegment == 10)
			{
				double tmp = 9.0 / 89600.0 / 9.0;
				spatialIntegrationWeights[0] = 2857.0 * tmp;		spatialIntegrationWeights[9] = 2857.0 * tmp;
				spatialIntegrationWeights[1] = 15741.0 * tmp;		spatialIntegrationWeights[8] = 15741.0 * tmp;
				spatialIntegrationWeights[2] = 1080.0 * tmp;		spatialIntegrationWeights[7] = 1080.0 * tmp;
				spatialIntegrationWeights[3] = 19344.0 * tmp;		spatialIntegrationWeights[6] = 19344.0 * tmp;
				spatialIntegrationWeights[4] = 5778.0 * tmp;		spatialIntegrationWeights[5] = 5778.0 * tmp;
			}
			else if (numSpatialPointsPerSegment == 11)
			{
				double tmp = 5.0 / 299376.0 / 10.0;
				spatialIntegrationWeights[0] = 16067.0 * tmp;		spatialIntegrationWeights[10] = 16067.0 * tmp;
				spatialIntegrationWeights[1] = 106300.0 * tmp;		spatialIntegrationWeights[9] = 106300.0 * tmp;
				spatialIntegrationWeights[2] = -48525.0 * tmp;		spatialIntegrationWeights[8] = -48525.0 * tmp;
				spatialIntegrationWeights[3] = 272400.0 * tmp;		spatialIntegrationWeights[7] = 272400.0 * tmp;
				spatialIntegrationWeights[4] = -260550.0 * tmp;		spatialIntegrationWeights[6] = -260550.0 * tmp;
				spatialIntegrationWeights[5] = 427368.0 * tmp;
			}
		}
		else
		{
			double H = 2.0 * h;
			double ovsh = 1.0 / 6.0 * H;
			spatialIntegrationWeights[0] = ovsh;
			int numPtsHalf = (numSpatialPointsPerSegment - 1) / 2;
			int cntr = 1;
			for (int i = 0; i < numPtsHalf; ++i)
			{
				spatialIntegrationWeights[cntr++] = 4.0 * ovsh;
				spatialIntegrationWeights[cntr++] = 2.0 * ovsh;
			}
			if (numSpatialPointsPerSegment % 2 == 0) // even segments
			{
				spatialIntegrationWeights[numSpatialSubsegments_BulkInterfacePoints_Print] = 0.5 * h;
				spatialIntegrationWeights[numSpatialSubsegments_BulkInterfacePoints_Print - 1] = ovsh + 0.5 * h;
			}
			else
				spatialIntegrationWeights[numSpatialSubsegments_BulkInterfacePoints_Print] = ovsh;
		}
	}
	return numSpatialPointsPerSegment;
}

double rlog10(double val)
{
	static double factor = 1.0 / log(10.0);
	if (val > 0.0)
		return log(val) * factor;
	return std::numeric_limits<double>::quiet_NaN();
}

void Test_SetNewtonCotes_Points_AndWeights()
{
	bool useRepeatedSimpsonRuleForHigherOrders = false;
	for (int numSpatialSubsegments_BulkInterfacePoints_Print = 1; numSpatialSubsegments_BulkInterfacePoints_Print < 20; ++numSpatialSubsegments_BulkInterfacePoints_Print)
	{ 
		vector<double> spatialIntegrationWeights, spatialIntegrationPoints;
		int numSpatialPointsPerSegment = SetNewtonCotes_Points_AndWeights(numSpatialSubsegments_BulkInterfacePoints_Print, useRepeatedSimpsonRuleForHigherOrders, spatialIntegrationWeights, spatialIntegrationPoints);
		double sm = 0.0;
		for (int i = 0; i < numSpatialPointsPerSegment; ++i)
			sm += spatialIntegrationWeights[i];
		if (fabs(sm - 1.0) > 1e-5)
		{
			for (int i = 0; i < numSpatialPointsPerSegment; ++i)
				cout << spatialIntegrationWeights[i] << '\t';
			cout << "numSpatialSubsegments_BulkInterfacePoints_Print\t" << numSpatialSubsegments_BulkInterfacePoints_Print << '\n';
			cout << "numSpatialPointsPerSegment\t" << numSpatialPointsPerSegment << '\n';
			cout << "\nsum\t" << sm << '\n';
			THROW("Sum is not one\t");
		}
	}
	cout << "success\n";
	getchar();
}

double computeRatio(double numerator, double denominator)
{
	static double max_ret = 1e40;
	static double zr = 1e10 * DBL_MIN, zrn = 100 * DBL_MIN;
	if (fabs(denominator) < zr)
	{
		if (fabs(numerator) < zr) //zrn)
		{
			return 0;
		}
		else
		{
			if (((numerator > 0) && (denominator > 0)) || ((numerator < 0) && (denominator < 0)))
				return max_ret;
			else
				return -max_ret;
		}
	}
	double ratio = numerator / denominator;

	if (ratio > max_ret)
		return max_ret;
	if (ratio < -max_ret)
		return -max_ret;

	return ratio;
}

bool DoublesAreEqual(double d1, double d2, double tol)
{
	double den = fabs(d1);
	if (fabs(d2) > den)	den = fabs(d2);
	if (den < 1.0) den = 1.0;
	return (fabs((d1 - d2) / tol) < den);
}

// returns size
unsigned int BreakString(const string& inString, vector<string>& parts)
{
	//#if 0
	std::stringstream ss(inString);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	unsigned int numF = vstrings.size();
	parts.resize(numF);
	for (unsigned int i = 0; i < numF; ++i)
		parts[i] = vstrings[i];
	return numF;
	//#endif
	//return 0;
}

unsigned int BreakStringBySeparator(const string& inString, vector<string>& parts, char separator)
{
	std::istringstream ss(inString);
	string token;
	parts.clear();
	while (std::getline(ss, token, ',')) {
		parts.push_back(token);
	}
	return parts.size();
}

std::string removeExtension(const std::string& filename)
{
	size_t lastdot = filename.find_last_of(".");
	if (lastdot == std::string::npos) return filename;
	return filename.substr(0, lastdot);
}
