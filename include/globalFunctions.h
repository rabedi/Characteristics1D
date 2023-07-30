#ifndef GLOBAL_FUNCTIONS__H
#define GLOBAL_FUNCTIONS__H

#include "globalMacros.h"
#include <float.h>
#include <iterator>

// return is numSpatialPointsPerSegment
// points are between 0 and 1
// useRepeatedSimpsonRuleForHigherOrders: for higher orders (numSpatialSubsegments_BulkInterfacePoints_Print >= 3), it uses repeated Simpson (+Trapezoidal rule)
int SetNewtonCotes_Points_AndWeights(int numSpatialSubsegments_BulkInterfacePoints_Print, bool useRepeatedSimpsonRuleForHigherOrders, vector<double>& spatialIntegrationWeights, vector<double>& spatialIntegrationPoints);

double computeRatio(double numerator, double denominator);
bool DoublesAreEqual(double d1, double d2, double tol = 1e-5);

// returns size
unsigned int BreakString(const string& inString, vector<string>& parts);
unsigned int BreakStringBySeparator(const string& inString, vector<string>& parts, char separator = ',');
std::string removeExtension(const std::string& filename);

void Test_SetNewtonCotes_Points_AndWeights();
#endif