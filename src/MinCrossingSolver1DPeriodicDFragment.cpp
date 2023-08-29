#include "MinCrossingSolver1DPeriodicDFragment.h"
#include "SL_OneInterfaceAllTimes.h"

gFx2yPDF::gFx2yPDF()
{
	tsrModel = tsr_Ortiz;
	para0_is_la = true;
}

void gFx2yPDF::InitializeValues(const map<string, string>& str_mapIn, const vector<unsigned int>& indices4ParasIn, const vector<double>& parasIn, double& xMin, double& xMax, double& tol_x, vector<double>& primary_xs, vector<double>& secondary_xs, int& num_y, vector<double>& tol_ys)
{
	indices4Paras = indices4ParasIn;
	str_map = str_mapIn;
	map<string, string>::iterator it;
	it = str_map.find("tsrModel");
	if ((it != str_map.end()) && (!string2Type(it->second, tsrModel)))
	{
		cout << "it->second\t" << it->second << '\n';
		THROW("Cannot find the appropriate cohesive model for solution\n");
	}
	it = str_map.find("para0_is_la");
	int i_para0_is_la = (int)para0_is_la;
	if ((it != str_map.end()) && (!fromString(it->second, i_para0_is_la)))
	{
		cout << "it->second\t" << it->second << '\n';
		THROW("Boolean not successfully read!\n");
	}
	para0_is_la = (bool)i_para0_is_la;

	paras = parasIn;
	Periodic1IntrFrag p1Conf;
	Writel10_al_TSR_file();
	p1Conf.Set_TSR_Model(tsrModel);
	p1Conf.Initialize_Periodic1IntrFrag(true);

	// log 10's of l's
	p1Conf.GetPrimaryFragmentSizes(primary_xs);
	for (unsigned int i = 0; i < primary_xs.size(); ++i)
		primary_xs[i] = log(primary_xs[i]) / log(10.0);

	double l10_l_dilute = primary_xs[fst_dilute];
	double l10_l_dilute_approx = primary_xs[fst_dilute_approx];
	double l10_l_Zhu6a = primary_xs[fst_zhu6a];
	double l10_l_Zhu6b = primary_xs[fst_zhu6b];
	double l10_l_Glenn = primary_xs[fst_Glenn];
	double l10_l_Grady = primary_xs[fst_Grady];

//	double xCenter = l10_l_dilute + l10_l_Zhu6a + l10_l_Zhu6b + l10_l_Glenn;
//	xCenter /= 4.0;
	double xm = l10_l_dilute;
	double xM = l10_l_dilute;
	xm = MIN(xm, l10_l_Zhu6a);
	xm = MIN(xm, l10_l_Zhu6b);
	xm = MIN(xm, l10_l_Glenn);

	xM = MAX(xM, l10_l_Zhu6a);
	xM = MAX(xM, l10_l_Zhu6b);
	xM = MAX(xM, l10_l_Glenn);
	xM += 1.5; //0.04;//*= 1.1;
	xm -= 0.04; //*= 0.9;
	unsigned int numSeg = 20; // 10;
	double dx = (xM - xm) / (double)numSeg;
	secondary_xs.resize(numSeg + 1);
	for (unsigned int i = 0; i <= numSeg; ++i)
		secondary_xs[i] = xm + i * dx;


	tol_x = 0.01; // start with this, this means the error in log is 0.01 -> relative error in x (fragment size) is 0.01
	tol_x = MIN(tol_x, dx * 0.03);
	xMin = xm - 0.5;
	xMax = xM + 0.5;

	num_y = 105;
}

bool gFx2yPDF::ComputeValue(genIndexVal& giv)
{
	string configName = "config/OneInterface/SampleConfig_Periodic_Segment.txt";
	double l10l = giv.x;
	int index_primary_l = giv.index_main;
	int index_secondary_l = giv.index_sec;
	Writel10_al_TSR_file(l10l, index_primary_l, index_secondary_l);
	MAIN_SL_OneInterfaceAllTimes_ONE_Interface(configName);

	int aIndex = 0;
	if (indices4Paras.size() > 0)
		aIndex = indices4Paras[0];

	int lIndex = index_primary_l;
	string aIndex_s, lIndex_s;
	toString(aIndex, aIndex_s);
	toString(lIndex, lIndex_s);
	string nameBase = "periodic1IntrFrag";
	string nameOne_a_SharedAll_l = nameBase + "_aI_" + aIndex_s;
	string ser;
	int aIndex_secondary = 0;
	toString(aIndex_secondary, ser);
	nameOne_a_SharedAll_l += "_";
	nameOne_a_SharedAll_l += ser;

	string nameOne_a_One_l = nameOne_a_SharedAll_l + "_lI_" + lIndex_s;
	toString(index_secondary_l, ser);
	nameOne_a_One_l += "_";
	nameOne_a_One_l += ser;
	string nameOneRun = nameOne_a_One_l + ".yvals";

	giv.ys.clear();
	fstream in(nameOneRun.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << nameOneRun << '\n';
		THROW("Cannot open output file\n");
	}
	double tmp;
	string str;
	while (true)
	{
		in >> str;
		if (!in.eof())
		{
			if (!fromString(str, tmp))
				tmp = std::numeric_limits<double>::quiet_NaN();
			giv.ys.push_back(tmp);
		}
		else
			break;
	}
	return true;
}

void gFx2yPDF::Print_YHeader(ostream& out, int num_y) const
{
	Periodic1IntrFrag::Output_Periodic1IntrFrag_Header(out);
}

void gFx2yPDF::Writel10_al_TSR_file(double l10l, int index_primary_l, int index_secondary_l)
{
	static double log10_2 = log(2.0) / log(10.0);
	double l10a = 2.0;
	int index_primary_a = 0, index_secondary_a = 0;
	if (paras.size() > 0)
	{
		l10a = paras[0];
		if (!para0_is_la)
			l10a += log10_2;
		if (indices4Paras.size() > 0)
			index_primary_a = indices4Paras[0];
	}
	fstream outla("l10a.txt", ios::out);
	outla << l10a << '\t';
	outla << index_primary_a << '\t';
	outla << index_secondary_a << '\n';

	fstream outll("l10l.txt", ios::out);
	outll << l10l << '\t';
	outll << index_primary_l << '\t';
	outll << index_secondary_l << '\n';

	fstream outtsr("tsr.txt", ios::out);
	outtsr << tsrModel;
}

void MAIN_SolvePeriodicFragmentSize(string confSolve)
{
	gFx2yPDF* functionIn = new gFx2yPDF();
	MAIN_Solver1D(confSolve, (gFx2y*)functionIn);
	delete functionIn;
	cout << "solution successful\n";
}
