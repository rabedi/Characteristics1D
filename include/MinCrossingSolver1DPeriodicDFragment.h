#ifndef MICRO_CROSSING_SOLVER_PERIODIC_D_FRAGMENT__H
#define MICRO_CROSSING_SOLVER_PERIODIC_D_FRAGMENT__H

#include "MinCrossingSolver1D.h"
#include "SLInterfaceCalculator.h"


class gFx2yPDF : public gFx2y
{
public:
	gFx2yPDF();
	virtual void InitializeValues(const map<string, string>& str_mapIn, const vector<unsigned int>& indices4ParasIn, const vector<double>& parasIn, double& xMin, double& xMax, double& tol_x, vector<double>& primary_xs, vector<double>& secondary_xs, int& num_y, vector<double>& tol_ys);
	virtual bool ComputeValue(genIndexVal& giv);
	virtual void Print_YHeader(ostream& out, int num_y = -1) const;
	SLFF_TSRType	tsrModel;
	// if true,		para0 =	log10(a)
	// else			para0 = log10(a_p), where a_p = a / 2
	// a is nondimensional strain rate. a_p is the version used by Zhu, Drugan, etc, mine is a
	bool para0_is_la;

	// writes log 10 files for a and l
	void Writel10_al_TSR_file(double l10l = 1.0, int index_primary_l = 0, int index_secondary_l = 0);
};

void MAIN_SolvePeriodicFragmentSize(string confSolve = "config/OneInterface/Periodic_segment/Periodic_segmen_confSolver1D.txt");

#endif

