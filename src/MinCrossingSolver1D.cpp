#include "MinCrossingSolver1D.h"
#include "commonMacros.h"
#include "commonFunctions.h"

bool genIndexVal_DB::g_genIndexVal_DB_equality_by_x = true;
double genIndexVal_DB::g_genIndexVal_DB_tol_x = 1e-6;

int Crossing_Sln_Mode(convergeT1D dat)
{
	if (dat == ct1d_yes)
		return 1;
	if (dat == ct1d_notCrossing)
		return 0;
	if (dat == ct1d_xmin_lim)
		return 0;
	if (dat == ct1d_xmax_lim)
		return 0;
	if (dat == ct1d_failedSln)
		return -2;
	if (dat == ct1d_nanSln)
		return -1;
	if (dat == ct1d_gettingSameXcrd)
		return 1;
	if (dat == ct1d_large_delx)
		return 0;
	if (dat == ct1d_large_dely)
		return 0;
	THROW("Invalid option");
}

int Extremum_Sln_Mode(convergeT1D dat) 
{
	if (dat == ct1d_yes)
		return 1;
	if (dat == ct1d_notCrossing)
		return 0;
	if (dat == ct1d_xmin_lim)
		return 1;
	if (dat == ct1d_xmax_lim)
		return 1;
	if (dat == ct1d_failedSln)
		return -2;
	if (dat == ct1d_nanSln)
		return -1;
	if (dat == ct1d_gettingSameXcrd)
		return 1;
	if (dat == ct1d_large_delx)
		return 0;
	if (dat == ct1d_large_dely)
		return 0;
	THROW("Invalid option");
}

int Sln_Mode(convergeT1D dat, bool isExtremum)
{
	if (isExtremum)
		return Extremum_Sln_Mode(dat);
	return Crossing_Sln_Mode(dat);
}

string getName(convergeT1D dat)
{
	if (dat == ct1d_yes)
		return "yes";
	if (dat == ct1d_notCrossing)
		return "notCrossing";
	if (dat == ct1d_xmin_lim)
		return "xmin_lim";
	if (dat == ct1d_xmax_lim)
		return "xmax_lim";
	if (dat == ct1d_failedSln)
		return "failedSln";
	if (dat == ct1d_nanSln)
		return "nanSln";
	if (dat == ct1d_gettingSameXcrd)
		return "gettingSameXcrd";
	if (dat == ct1d_large_delx)
		return "large_delx";
	if (dat == ct1d_large_dely)
		return "large_dely";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

bool name2Type(string& name, convergeT1D& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= convergeT1D_SIZE)
			THROW("too large of a number\n");
		typeVal = (convergeT1D)num;
		return true;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < convergeT1D_SIZE; ++i)
	{
		typeVal = (convergeT1D)i; // casting integer to convergeT1D, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return true;
	}
	return false;
}

//operator for output
ostream& operator<<(ostream& out, convergeT1D dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, convergeT1D& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	if (name2Type(name, dat) == false)
	{
		cout << name << '\n';
		THROW("Invalid name\n");
	}
	return in;
}

ostream& operator<<(ostream& out, const ConvergenceLog& dat)
{
	dat.Print(out, true, true, true, false, -1);
	return out;
}

void ConvergenceLog::Print(ostream& out, bool printHeader, bool printData, bool printLong, bool expandConvType, int isExtremumSln) const
{
	if (printHeader)
		out << "convType\t";
	if (printData)
		out << convType << "\t";
	if (expandConvType)
	{
		if (printHeader)
			out << "i_convType\t";
		if (printData)
			out << (int)convType << "\t";
		if (isExtremumSln != -1)
		{
			if (printHeader)
				out << "convType_valid\t";
			if (printData)
			{
				bool isExtremum = (bool)isExtremumSln;
				int slnMode = Sln_Mode(convType, isExtremum);
				out << slnMode << "\t";
			}
		}
	}
	if (printHeader)
		out << "iterCnt\t";
	if (printData)
		out << iterCnt << "\t";
	if (printHeader)
		out << "iterExtremumNotChanged\t";
	if (printData)
		out << iterExtremumNotChanged << "\t";
	if (printHeader)
		out << "delxs\t";
	if (printData)
	{
		unsigned int sz = delxs.size();
		if (printLong)
		{
			out << "sz\t" << sz << "\t";
			for (unsigned int i = 0; i < delxs.size(); ++i)
				out << delxs[i] << '\t';
		}
		else
		{
			double lastVal = -1;
			if (sz > 0)
				lastVal = delxs[sz - 1];
			out << lastVal << '\t';
		}
	}

	if (printHeader)
		out << "delys\t";
	if (printData)
	{
		unsigned int sz = delys.size();
		if (printLong)
		{
			out << "sz\t" << sz << "\t";
			for (unsigned int i = 0; i < delys.size(); ++i)
				out << delys[i] << '\t';
		}
		else
		{
			double lastVal = -1;
			if (sz > 0)
				lastVal = delys[sz - 1];
			out << lastVal << '\t';
		}
	}
}


ConvergenceLog::ConvergenceLog()
{
	convType = ct1d_yes;
	iterCnt = 0;
	iterExtremumNotChanged = 0;
}

void gFx2y::InitializeValues(Solver1D* conf, const map<string, string>& str_mapIn, const vector<unsigned int>& indices4ParasIn, const vector<double>& parasIn, double& xMin, double& xMax, double& tol_x, vector<double>& primary_xs, vector<double>& secondary_xs, int& num_y, vector<double>& tol_ys)
{
	double tol = 1e-2;
	str_map = str_mapIn;
	paras = parasIn;
	indices4Paras = indices4ParasIn;
	double x0 = 2.0;
	if (paras.size() > 0)
		x0 = paras[0];
	double del = 2.0 - x0;
	xMin = 0.0 + del;
	xMax = 6.0 + del;
	tol_x = tol;
	num_y = 2;
	tol_ys.resize(num_y);
	for (unsigned int i = 0; i < (unsigned int)num_y; ++i)
		tol_ys[i] = tol * 1e-3;

	primary_xs.clear();
	primary_xs.push_back(0.34 + del);
	primary_xs.push_back(1.21 + del);
	primary_xs.push_back(4.41 + del);

//	double xm = 0.2, xM = 5.7, del_x = 0.31; // 0.111; // 5.4
//	double xm = 0.2, xM = 6.7, del_x = 0.31; // 0.111; // 5.4
	double xm = 0.2, xM = 4.4, del_x = 0.31; // 0.111; // 5.4
//	double xm = 1.0, xM = 2.0, del_x = 0.31; // 0.111; // 5.4
	unsigned int sz = (unsigned int)(floor((xM - xm) / del_x + 1e-4));
	secondary_xs.resize(sz);
	for (unsigned int i = 0; i < sz; ++i)
		secondary_xs[i] = xm + del_x * i + +del;
}

bool gFx2y::ComputeValue(genIndexVal& giv)
{
	giv.ys.resize(2);
	double x0 = 2.0;
	if (paras.size() > 0)
		x0 = paras[0];
	double del = 2.0 - x0;

	double tmp = (giv.x - x0);
	giv.ys[0] = tmp * tmp  - 9.0;
	giv.ys[1] = sin(2.0 * giv.x);
	if (giv.index_main == 1)
	{
		int tmp = (int)floor(giv.x * 100);
		if (tmp % 2 == 0)
			giv.ys[1] = std::numeric_limits<double>::quiet_NaN();
	}
	return true;
}

void gFx2y::Print_Header(ostream& out, int num_y) const
{
	if (num_y <= 0)
		return;
	out << "indexMain\tindexSecondary\tx\t";
	Print_YHeader(out, num_y);
}

void gFx2y::Print_YHeader(ostream& out, int num_y) const
{
	out << "out0";
	for (unsigned int i = 1; i < (unsigned int)num_y; ++i)
		out << "\tout" << i;
}

sgnT Sgn(double val, double tol)
{
	if (val > tol)
		return sgn_p;
	if (val < tol)
		return sgn_m;
	return sgn_0;
}

genIndexVal::genIndexVal(int index_main_in, double x_in)
{
	index_main = index_main_in;
	x = x_in;
}

bool genIndexVal::operator<(const genIndexVal &other) const
{
	if (genIndexVal_DB::g_genIndexVal_DB_equality_by_x)
	{
		if (*this == other)
			return false;
		return (x < other.x);
	}
	else
	{
		if (index_main < other.index_main)
			return true;
		if (index_main > other.index_main)
			return false;
		if (index_sec <= other.index_sec)
			return true;
		return false;
	}
}
bool genIndexVal::operator==(const genIndexVal &other) const
{
	if (genIndexVal_DB::g_genIndexVal_DB_equality_by_x)
	{
		double del = fabs(x - other.x);
		return (del < genIndexVal_DB::g_genIndexVal_DB_tol_x);
#if 0
		double maxAbsVal = MAX(fabs(x), fabs(other.x));
		if (maxAbsVal < genIndexVal_DB::g_genIndexVal_DB_tol_x)
			return (del < genIndexVal_DB::g_genIndexVal_DB_tol_x * genIndexVal_DB::g_genIndexVal_DB_tol_x);
		return (del < genIndexVal_DB::g_genIndexVal_DB_tol_x * maxAbsVal);
#endif
	}
	return ((index_main == other.index_main) && (index_sec == other.index_sec));
}

istream& operator>>(istream& in, genIndexVal& dat)
{
	dat.Read_genIndexVal(in);
	return in;
}

ostream& operator<<(ostream& out, const genIndexVal& dat)
{
	dat.Write_genIndexVal(out);
	return out;
}

void genIndexVal::Write_genIndexVal(ostream& out) const
{
	out << "index_main\t" << index_main << '\t';
	out << "index_sec\t" << index_sec << '\t';
	out << "x\t" << x << '\t';
	out << "ys\tsz\t" << ys.size() << '\t';
	for (unsigned int i = 0; i < ys.size(); ++i)
		out << ys[i] << '\t';
	out << '\n';
}

void genIndexVal::Write_genIndexVal_JustVals(ostream& out) const
{
	out << index_main << '\t';
	out << index_sec << '\t';
	out << x << '\t';
	for (unsigned int i = 0; i < ys.size(); ++i)
		out << ys[i] << '\t';
	out << '\n';
}

bool genIndexVal::Read_genIndexVal(istream& in)
{
	string buf;
	in >> buf >> index_main;
	if (in.eof())
		return false;
	in >> buf >> index_sec;
	in >> buf >> x;
	unsigned int sz;
	in >> buf >> buf >> sz;
	ys.resize(sz);
	for (unsigned int i = 0; i < sz; ++i)
	{
		in >> buf;
		if (fromString(buf, ys[i]) == false)
			ys[i] = std::numeric_limits<double>::quiet_NaN();
	}
	return true;
}

void genIndexVal::MakeNaN()
{
	x = std::numeric_limits<double>::quiet_NaN();
	for (unsigned int i = 0; i < ys.size(); ++i)
		ys[i] = std::numeric_limits<double>::quiet_NaN();
 }

genIndexVal_DBH::genIndexVal_DBH()
{
	pos = -1;
	y_size = 0;
//	givPtr = NULL;
	found = false;
}

void genIndexVal_DB::Sort(bool equality_by_x, double tolx, vector<genIndexVal>* dbIn)
{
	bool bk_b = g_genIndexVal_DB_equality_by_x;
	g_genIndexVal_DB_equality_by_x = equality_by_x;
	if (tolx < 0.0)
		tolx = g_genIndexVal_DB_tol_x;
	double bktolx = g_genIndexVal_DB_tol_x;
	g_genIndexVal_DB_tol_x = tolx;

	if (dbIn == NULL)
		dbIn = &db;
	sort(dbIn->begin(), dbIn->end());

	g_genIndexVal_DB_equality_by_x = bk_b;
	g_genIndexVal_DB_tol_x = bktolx;
}

void genIndexVal_DB::CopySort(vector<genIndexVal>& dbOut, bool equality_by_x, double tolx)
{
	dbOut = db;
	Sort(equality_by_x, tolx, &dbOut);
}

void genIndexVal_DB::Write_genIndexVal_DB(ostream& out, bool print_dbAsIs, bool equality_by_x, double tolx)
{
	vector<genIndexVal>* db2Print = &db;
	vector<genIndexVal> dbOut;
	if (!print_dbAsIs)
	{
		CopySort(dbOut, equality_by_x, tolx);
		db2Print = &dbOut;
	}
	unsigned int sz = db2Print->size();
	for (unsigned int i = 0; i < sz; ++i)
		(*db2Print)[i].Write_genIndexVal(out);
}

void genIndexVal_DB::Write_genIndexVal_DB_WithHeader(gFx2y* functionIn, ostream& out, bool print_dbAsIs, bool equality_by_x, double tolx)
{
	vector<genIndexVal>* db2Print = &db;
	vector<genIndexVal> dbOut;
	if (!print_dbAsIs)
	{
		CopySort(dbOut, equality_by_x, tolx);
		db2Print = &dbOut;
	}
	// printing the header
	out << "index_main\t";
	out << "index_sec\t";
	out << "x\t";
	unsigned int num_y = 0;
	unsigned int sz = db2Print->size();
	if (sz > 0)
		num_y = db[0].ys.size();
	functionIn->Print_YHeader(out, num_y);
	out << '\n';
	for (unsigned int i = 0; i < sz; ++i)
		(*db2Print)[i].Write_genIndexVal_JustVals(out);
}

void genIndexVal_DB::Read_genIndexVal_DB(istream& in, bool read_dbAsIs, bool equality_by_x, double tolx)
{
	db.clear();
	while (true)
	{
		genIndexVal giv;
		if (giv.Read_genIndexVal(in) == false)
			break;
		db.push_back(giv);
	}
	if (!read_dbAsIs)
		Sort(equality_by_x, tolx);
}

genIndexVal_DBH genIndexVal_DB::Add_Pt(genIndexVal& addPt, bool checkIfPtExists, bool equality_by_x, double tolx)
{
	genIndexVal_DBH retVal;
	retVal.found = false;
	unsigned int sz = db.size();
	retVal.pos = sz;
	retVal.y_size = addPt.getSize_y();
	if (!checkIfPtExists)
	{
		db.push_back(addPt);
//		retVal.givPtr = &db[sz];
		return retVal;
	}
	// part 1, set sort_x and tol_x to vals set
	bool bk_b = g_genIndexVal_DB_equality_by_x;
	g_genIndexVal_DB_equality_by_x = equality_by_x;
	if (tolx < 0.0)
		tolx = g_genIndexVal_DB_tol_x;
	double bktolx = g_genIndexVal_DB_tol_x;
	g_genIndexVal_DB_tol_x = tolx;

	// part 2, do the actual calculation
	bool found = false;
	unsigned int i;
	for (i = 0; i < sz; ++i)
	{
		if (addPt == db[i])
		{
			found = true;
			break;
		}
	}
	if (!found)
	{
		db.push_back(addPt);
//		retVal.givPtr = &db[sz];
		return retVal;
	}
	retVal.found = true;
	// now the point exist, we want to use the one that has more y's computed
	retVal.pos = i;
	unsigned int existing_n_ys = db[i].getSize_y();
	if (existing_n_ys > retVal.y_size)
	{
		addPt = db[i];
		retVal.y_size = existing_n_ys;
//		retVal.givPtr = &db[i];
		return retVal;
	}
	// existing y size is smaller ...
	db[i] = addPt;
//	retVal.givPtr = &db[i];
	return retVal;

	// part 3, revert sort_x and tol_x to default vals
	g_genIndexVal_DB_equality_by_x = bk_b;
	g_genIndexVal_DB_tol_x = bktolx;
}

istream& operator>>(istream& in, Solver1D_1posConf& dat)
{
	dat.Read_Solver1D_1posConf(in);
	return in;
}

bool Solver1D_1posConf::Read_Solver1D_1posConf(istream& in)
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
		if (buf == "nameBase")
		{
			READ_NSTRING(in, buf, nameBase);
		}
		else if (buf == "isActive")
		{
			READ_NBOOL(in, buf, isActive);
		}
		else if (buf == "pos_y")
		{
			READ_NINTEGER(in, buf, pos_y);
		}
		else if (buf == "tol_y")
		{
			READ_NDOUBLE(in, buf, tol_y);
		}
		else if (buf == "isExtremum")
		{
			READ_NBOOL(in, buf, isExtremum);
		}
		else if (buf == "isMins")
		{
			vector<int> iv;
			ReadVectorInteger(in, iv);
			unsigned int sz = iv.size();
			isMins.resize(sz);
			for (unsigned int i = 0; i < sz; ++i)
				isMins[i] = (bool)iv[i];
		}
		else if (buf == "crossing_ys")
		{
			ReadVectorDouble(in, crossing_ys);
		}
		else if (buf == "x0s")
		{
			ReadVectorDouble(in, x0s);
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

Solver1D_1posConf::Solver1D_1posConf()
{
	nameBase = "none";
	pos_y = 0;
	isExtremum = true;
	tol_y = -1.0;
	isActive = true;
}


unsigned int Solver1D_1posConf::Initialize_Solver1D_1posConf()
{
	sz = 0;
	if (isExtremum)
	{
		sz = isMins.size();
		if (sz == 0)
		{
			isMins.push_back(true);
			sz = 1;
		}
	}
	else
	{
		sz = crossing_ys.size();
		if (sz == 0)
		{
			crossing_ys.push_back(0.0);
			sz = 1;
		}
	}
	string pos_y_s;
	toString(pos_y, pos_y_s);
	if ((nameBase == "none") || (nameBase == ""))
	{
		nameBase = "pos_y_" + pos_y_s;
		if (isExtremum)
			nameBase += "_extreme";
		else
			nameBase += "_crossing";
	}
	vec_name.resize(sz);
	for (unsigned int i = 0; i < sz; ++i)
		vec_name[i] = nameBase;
	if (isExtremum)
	{
		for (unsigned int i = 0; i < sz; ++i)
		{
			if (isMins[i])
				vec_name[i] += "_min";
			else
				vec_name[i] += "_max";
		}
	}
	else if (sz > 1)
	{
		for (unsigned int i = 0; i < sz; ++i)
		{
			string ser;
			toString(i, ser);
			vec_name[i] += "_I_";
			vec_name[i] += ser;
		}
	}
	vec_isExtremum.resize(sz);
	fill(vec_isExtremum.begin(), vec_isExtremum.end(), isExtremum);
	if (isExtremum)
	{
		crossing_ys.resize(sz);
		fill(crossing_ys.begin(), crossing_ys.end(), 0.0);
	}
	else
	{
		isMins.resize(sz);
		fill(isMins.begin(), isMins.end(), true);
	}
	tol_ys.resize(sz);
	fill(tol_ys.begin(), tol_ys.end(), tol_y);
	unsigned int sz_x0s = x0s.size();
	if (sz_x0s < sz)
	{
		for (unsigned int i = sz_x0s; i < sz; ++i)
			x0s.push_back(-1e15);
	}
	return sz;
}

Solver1D::Solver1D()
{
	startInds[0] = 0;
	startInds[1] = 0;
	startInds[2] = 0;
	startInds[3] = 0;
	baseName = "_Solver1D";
	read_Results = true;
	write_Results = true;
	equality_by_x = true;
	tol_x = 1e20 * DBL_MIN;
	num_y = 0;
	maxNumIter = 100;
	maxNumIterExtremumNotChanging = 8;
	del_secondary_x = -1;

	do_posConfs_first_notAddPtSolve = true;
	do_posConfs_AddPtSolve = true;
	do_posConfs_second_notAddPtSolve = true;
	b_print_PrimaryPoints = true;
}

Solver1D::~Solver1D()
{
	Store_pts_unsorted_SlnEnd();
}

void Solver1D::Read_Solver1D(istream& in)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "baseName")
		{
			READ_NSTRING(in, buf, baseName);
		}
		else if (buf == "b_print_PrimaryPoints")
		{
			READ_NBOOL(in, buf, b_print_PrimaryPoints);
		}
		else if (buf == "do_posConfs_first_notAddPtSolve")
		{
			READ_NBOOL(in, buf, do_posConfs_first_notAddPtSolve);
		}
		else if (buf == "do_posConfs_AddPtSolve")
		{
			READ_NBOOL(in, buf, do_posConfs_AddPtSolve);
		}
		else if (buf == "do_posConfs_second_notAddPtSolve")
		{
			READ_NBOOL(in, buf, do_posConfs_second_notAddPtSolve);
		}
		else if (buf == "maxNumIter")
		{
			READ_NINTEGER(in, buf, maxNumIter);
		}
		else if (buf == "maxNumIterExtremumNotChanging")
		{
			READ_NINTEGER(in, buf, maxNumIterExtremumNotChanging);
		}
		else if (buf == "read_Results")
		{
			READ_NBOOL(in, buf, read_Results);
		}
		else if (buf == "write_Results")
		{
			READ_NBOOL(in, buf, write_Results);
		}
		else if (buf == "equality_by_x")
		{
			READ_NBOOL(in, buf, equality_by_x);
		}
		else if (buf == "paras1D")
		{
			ReadVectorDouble(in, paras1D);
		}
		else if (buf == "indices4paras1D")
		{
			vector<int> v_tmpi;
			ReadVectorInteger(in, v_tmpi);
			unsigned int sz_v = v_tmpi.size();
			indices4paras1D.resize(sz_v);
			for (unsigned int i = 0; i < sz_v; ++i)
				indices4paras1D[i] = v_tmpi[i];
		}
		else if (buf == "str_map")
			ReadMapString2String(in, str_map);
		else if (buf == "posConfs2Solve")
		{
			posConfs2Solve.clear();
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				if (buf != "}")
					THROW("istream should start with {");
			}
			while (true)
			{
				Solver1D_1posConf spc;
				if (spc.Read_Solver1D_1posConf(in) == false)
					break;
				if (spc.isActive)
					posConfs2Solve.push_back(spc);
			}
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void Solver1D::InitializeFunction(gFx2y * functionIn)
{
	name_log = baseName;
	for (unsigned int i = 0; i < indices4paras1D.size(); ++i)
	{
		string ser;
		toString(indices4paras1D[i], ser);
		name_log += "_I";
		name_log += ser;
	}
	name_WOExt = name_log;
	string tmpFileName = name_log + ".err";
	out_error.open(tmpFileName.c_str(), ios::app);
	name_log_x_unsorted = name_log + "_i_sorted";
	name_log += "_x_sorted.log";

	function = functionIn;
	genIndexVal_DB::g_genIndexVal_DB_equality_by_x = equality_by_x;
	function->InitializeValues(this, str_map, indices4paras1D, paras1D, xMin, xMax, tol_x, primary_xs, secondary_xs, num_y, tol_ys);
	tol_x_eq_check = 1e-5 * tol_x;
	genIndexVal_DB::g_genIndexVal_DB_tol_x = tol_x_eq_check;

	if (read_Results)
	{
		bool read_dbAsIs = true, compute_ys_not_computed = true;
		Restore_pts(read_dbAsIs, compute_ys_not_computed);
	}
	bool checkIfPtExists = read_Results; // only when results are read it makes sense to see if the value is previously computed

	unsigned int sz_primary_xSet = primary_xs.size();
	unsigned int sz_secondary_xSet = secondary_xs.size();
	for (unsigned int i = 0; i < sz_primary_xSet; ++i)
	{
		double x = primary_xs[i];
		genIndexVal giv(0, x);
		Compute_Add_pt_value(giv, checkIfPtExists);
	}
	if (b_print_PrimaryPoints)
		Print_PrimaryPoints();
	for (unsigned int i = 0; i < sz_secondary_xSet; ++i)
	{
		double x = secondary_xs[i];
		genIndexVal giv(1, x);
		Compute_Add_pt_value(giv, checkIfPtExists);
	}
	if ((del_secondary_x < 0.0) && (sz_secondary_xSet > 1))
		del_secondary_x = secondary_xs[1] - secondary_xs[0];
}

ConvergenceLog Solver1D::Solve_x4_Crossing(genIndexVal& sln, double crossing_y, unsigned int y_pos, double x0, double tol_y)
{
	ConvergenceLog cl;
	if (tol_y < 0)
	{
		if (tol_ys.size() > y_pos)
			tol_y = tol_ys[y_pos];
		else
			tol_y = tol_x;
	}
	double tol_y_zero = MIN(0.01 * tol_y, 0.01 * tol_x);

	sgnT sgnjCrossing;
	unsigned int cntrAddedOnSides = 0;
	int jCrossing = -1, jCrossingm1 = -1;
	while (jCrossing < 0)
	{
		bool pointAddedOrCrossingHappened = Solve_x4_Crossing_Aux(cntrAddedOnSides++, crossing_y, y_pos, x0, tol_y_zero, jCrossingm1, jCrossing, sgnjCrossing, cl);
		if (cl.convType == ct1d_nanSln)
		{
			sln = pts.db[0];
			sln.MakeNaN();
			return cl;
		}
		if (!pointAddedOrCrossingHappened && (jCrossing < 0))
			return cl;
		if (jCrossing >= 0)
		{
			if (sgnjCrossing == sgn_0)
			{
				sln = pts.db[jCrossing];
				return cl;
			}
		}
	}

	double vL, vR, xL, xR;
	xL = pts.db[jCrossingm1].x;
	xR = pts.db[jCrossing].x;
	vL = pts.db[jCrossingm1].ys[y_pos] - crossing_y;
	vR = pts.db[jCrossing].ys[y_pos] - crossing_y;
	sgnT sgnL = (sgnT)(-1 * (sgnT)sgnjCrossing), sgnR = sgnjCrossing;
	double delx, yNew;

	cl.iterCnt = 0;
	bool checkIfPtExists = read_Results;
	
	double prev_xNew = 0.5 * (xL + xR);
	bool delx_good = true, dely_good = true;
	while (cl.iterCnt < maxNumIter)
	{
		double xNew = (vR * xL - vL * xR) / (vR - vL);
		sln.ys.clear();
		sln.x = xNew;
		sln.index_main = 3;
//		sln.index_sec = cl.iterCnt;
		genIndexVal_DBH retVal = Compute_Add_pt_value(sln, checkIfPtExists);
		if (isnan(sln.ys[y_pos]))
		{
			cl.convType = ct1d_nanSln;
			return cl;
		}
		if (retVal.pos < 0)
		{
			cl.convType = ct1d_failedSln;
			return cl;
		}
//		if (retVal.givPtr == NULL)
//			return false;
		yNew = sln.ys[y_pos] - crossing_y;
		sgnT sgn_y = Sgn(yNew, tol_y_zero);
		dely_good = (fabs(yNew) < tol_y);
		delx = fabs(prev_xNew - xNew);
		delx_good = (delx < tol_x);
		prev_xNew = xNew;
		cl.delxs.push_back(delx);
		cl.delys.push_back(yNew);

		if (sgn_y == sgn_0)
			break;
		if (retVal.found)
		{
			if (fabs(yNew) > 10.0 * tol_y)
			{
				//THROW("This is bug an algorithm and the code should return true here, perhaps need to tighten delx tolerance\n");
				cl.convType = ct1d_gettingSameXcrd;
				return cl;
			}
			break;
		}

		if (delx_good && dely_good)
			break;
		if (retVal.found)
		{
			cl.convType = ct1d_gettingSameXcrd;
			return cl;
//			THROW("This is bug an algorithm and the code should return true here\n");
//			return false;
		}
		if (sgn_y == sgnL)
		{
			xL = xNew;
			vL = yNew;
		}
		else
		{
			xR = xNew;
			vR = yNew;
		}
		++cl.iterCnt;
	}
	pts.Sort(true, tol_x_eq_check);
	if (cl.iterCnt >= maxNumIter)
	{
		if (!delx_good)
		{
			cl.convType = ct1d_large_delx;
			if ((!dely_good) && ((fabs(yNew) / tol_y) > (delx / tol_x)))
				cl.convType = ct1d_large_dely;
		}
		else
		{
			cl.convType = ct1d_large_dely;
			if (dely_good)
			{
				THROW("Not both convergence measures can be good");
			}
		}
		return cl;
	}
	return cl;
}

unsigned int Solver1D::Cross_Helper(double crossing_y, unsigned int y_pos, double tol_y_zero,
	vector<int>& poss, vector<sgnT>& signs, vector<cross_optimal_helper>& helpers)
{
	pts.Sort(true, tol_x_eq_check);
	unsigned int szPos = getValidValues(y_pos, poss);
	if (szPos < 2)
	{
		return szPos;
//		THROW("not enough points\n");
	}

	signs.resize(szPos);

	for (unsigned int ii = 0; ii < szPos; ++ii)
	{
		unsigned int i = poss[ii];
		signs[ii] = Sgn(pts.db[i].ys[y_pos] - crossing_y, tol_y_zero);
		if (signs[ii] == sgn_0)
		{
			cross_optimal_helper h;
			h.x = pts.db[i].x;
			h.xnext = h.x;
			h.i = i;
			h.ii = ii;
			h.signi = signs[ii];
			h.signinext = sgnT_none;
			helpers.push_back(h);
		}
	}
	// now see if there are other points after which the sign changes
	unsigned int szPosm1 = szPos - 1;
	sgnT signi, signinext;
	for (unsigned int ii = 0; ii < szPosm1; ++ii)
	{
		signi = signs[ii], signinext = signs[ii + 1];
		if ((signi == sgn_0) || (signinext == sgn_0) || (signi == signinext))
			continue;
		cross_optimal_helper h;
		h.i = poss[ii], h.ii = ii, h.inext = poss[ii + 1];
		h.signi = signi, h.signinext = signinext;
		h.x = pts.db[h.i].x;
		h.xnext = pts.db[h.inext].x;
		helpers.push_back(h);
	}
	return szPos;
}

void Solver1D::Print_PrimaryPoints()
{
	pts.Sort(false);
	for (unsigned int i = 0; i < pts.size(); ++i)
	{
		genIndexVal* giv = &pts.db[i];
		if (giv->index_main != 0)
			return;
		string ser;
		toString(giv->index_sec, ser);
		string fileName = baseName + "_PrimaryPt_" + ser + ".txt";
		fstream in(fileName.c_str(), ios::in);
		fstream out;
		if (!in.is_open())
		{
			out.open(fileName.c_str(), ios::app);
			PrintSln_Header(out);
		}
		else
		{
			in.close();
			out.open(fileName.c_str(), ios::app);
		}
		ConvergenceLog cl;
		cl.convType = ct1d_yes;
		PrintSln_Values(out, *giv, cl, false);
	}
}

bool Solver1D::Solve_x4_Crossing_Aux(unsigned int cntrAddedOnSides, double crossing_y, unsigned int y_pos, double x0, double tol_y_zero, int& jCrossingm1, int& jCrossing, sgnT& sgnjCrossing, ConvergenceLog& cl)
{
	vector<int> poss;
	vector<sgnT> signs;
	vector<cross_optimal_helper> helpers;

	unsigned int szPos = Cross_Helper(crossing_y, y_pos, tol_y_zero,
		poss, signs, helpers);
	if (szPos < 2)
	{
		cl.convType = ct1d_nanSln;
		return false;
	}

	unsigned int sz_helpers = helpers.size();
	if (sz_helpers == 1)
	{
		cross_optimal_helper h = helpers[0];
		if (h.signi == sgn_0)
		{
			sgnjCrossing = sgn_0;
			jCrossingm1 = h.i;
			jCrossing = jCrossingm1;
			return true;
		}
		sgnjCrossing = h.signinext;
		jCrossing = h.inext;
		jCrossingm1 = h.i;
		return true;
	}
	int state = 2; // 2: x0 between xMin and xMax;		 0, x0 < xMin;		1 x0 > xMax
	if (x0 < xMin)
	{
		x0 = xMin;
		state = 0;
	}
	if (x0 > xMax)
	{
		x0 = xMax;
		state = 1;
	}
	// now find whih the cloest point
	if (sz_helpers > 1)
	{
		unsigned int index_helper = 0;
		if (state == 1)
			index_helper = sz_helpers - 1;
		else if (state == 2)
		{
			map<double, int> dist2x0Map;
			double dist;
			for (unsigned int i = 0; i < sz_helpers; ++i)
			{
				dist = fabs(helpers[i].x - x0);
				dist2x0Map[dist] = i;
			}
			index_helper = dist2x0Map.begin()->second;
		}
		cross_optimal_helper h = helpers[index_helper];
		if (h.signi == sgn_0)
		{
			sgnjCrossing = sgn_0;
			jCrossingm1 = h.i;
			jCrossing = jCrossingm1;
			return true;
		}
		sgnjCrossing = h.signinext;
		jCrossing = h.inext;
		jCrossingm1 = h.i;
		return true;
	}
	/// now the size is zero and no crossing is found
	jCrossingm1 = -1;
	jCrossing = -1; // means crossing has not fonud
	unsigned int posL = poss[0], posR = poss[szPos - 1];
	double abs_yCrossing_L = fabs(pts.db[posL].ys[y_pos] - crossing_y);
	double abs_yCrossing_R = fabs(pts.db[posR].ys[y_pos] - crossing_y);
	double xL = pts.db[posL].x;
	double xR = pts.db[posR].x;
	if (del_secondary_x < 0.0)
	{
		cout << "del_secondary_x\t" << del_secondary_x << '\n';
		THROW("Need to have at least 2 secondary points or directly provide del_secondary_x\n");
	}
	double xLProposed = xL - del_secondary_x;
	double xRProposed = xR + del_secondary_x;
	bool canDecrease_xL = ((xLProposed >= xMin) && (posL == 0));
	bool canIncrease_xR = ((xRProposed <= xMax) && (posR == (pts.db.size() - 1)));
	int side = -1; // 0 -> L, 1 -> R
	if (canDecrease_xL)
	{
		side = 0;
		if ((canIncrease_xR) && (abs_yCrossing_L > abs_yCrossing_R))
			side = 1;
	}
	else
	{
		if (canIncrease_xR)
			side = 1;
	}
	if (side == -1)
	{
		if (abs_yCrossing_L <= abs_yCrossing_R)
			cl.convType = ct1d_xmin_lim;
		else
			cl.convType = ct1d_xmax_lim;
		return false;
	}
	genIndexVal giv;
	giv.index_main = 2;
//	giv.index_sec = cntrAddedOnSides;
	if (side == 0)
		giv.x = xLProposed;
	else
		giv.x = xRProposed;
	bool checkIfPtExists = true;
	genIndexVal_DBH retVal = Compute_Add_pt_value(giv, checkIfPtExists);
	if (isnan(giv.ys[y_pos]))
	{
		cl.convType = ct1d_nanSln;
		return false;
	}
	if (retVal.pos < 0)
	{
		cl.convType = ct1d_failedSln;
		return false;
	}
	return true;
}

ConvergenceLog Solver1D::Solve_x4_Crossing_NoPointAdded(genIndexVal& sln, double crossing_y, unsigned int y_pos, double x0, double tol_y)
{
	ConvergenceLog cl;

	if (tol_y < 0)
	{
		if (tol_ys.size() > y_pos)
			tol_y = tol_ys[y_pos];
		else
			tol_y = tol_x;
	}
	double tol_y_zero = MIN(0.01 * tol_y, 0.01 * tol_x);

	pts.Sort(true, tol_x_eq_check);
	vector<int> poss;
	unsigned int szPos = getValidValues(y_pos, poss);
	if (szPos < 2)
	{
		cl.convType = ct1d_nanSln;
		return cl;
//		THROW("not enough points\n");
	}

	vector<sgnT> signs(szPos);
	vector<double> vals(szPos);
	vector<cross_optimal_helper> helpers;

	for (unsigned int ii = 0; ii < szPos; ++ii)
	{
		unsigned int i = poss[ii];
		vals[ii] = pts.db[i].ys[y_pos] - crossing_y;
		signs[ii] = Sgn(vals[ii], tol_y_zero);
		if (signs[ii] == sgn_0)
		{
			cross_optimal_helper h;
			h.x = pts.db[i].x;
			h.xnext = h.x;
			h.i = i;
			h.ii = ii;
			h.signi = signs[ii];
			h.signinext = sgnT_none;
			helpers.push_back(h);
		}
	}
	// now see if there are other points after which the sign changes
	unsigned int szPosm1 = szPos - 1;
	sgnT signi, signinext;
	for (unsigned int ii = 0; ii < szPosm1; ++ii)
	{
		signi = signs[ii], signinext = signs[ii + 1];
		if ((signi == sgn_0) || (signinext == sgn_0) || (signi == signinext))
			continue;
		cross_optimal_helper h;
		h.i = poss[ii], h.ii = ii, h.inext = poss[ii + 1];
		h.signi = signi, h.signinext = signinext;
		h.x = pts.db[h.i].x;
		h.xnext = pts.db[h.inext].x;
		helpers.push_back(h);
	}
	unsigned int sz_helpers = helpers.size();
	cross_optimal_helper h;
	if (sz_helpers == 1)
		h = helpers[0];
	int state = 2; // 2: x0 between xMin and xMax;		 0, x0 < xMin;		1 x0 > xMax
	if (x0 < xMin)
	{
		x0 = xMin;
		state = 0;
	}
	if (x0 > xMax)
	{
		x0 = xMax;
		state = 1;
	}
	// now find whih the cloest point
	if (sz_helpers > 1)
	{
		unsigned int index_helper = 0;
		if (state == 1)
			index_helper = sz_helpers - 1;
		else if (state == 2)
		{
			map<double, int> dist2x0Map;
			double dist;
			for (unsigned int i = 0; i < sz_helpers; ++i)
			{
				dist = fabs(helpers[i].x - x0);
				dist2x0Map[dist] = i;
			}
			index_helper = dist2x0Map.begin()->second;
		}
		h = helpers[index_helper];
	}
	if (sz_helpers > 0) // crossing is found
	{
		if (h.signi == sgn_0)
		{
			sln = pts.db[h.i];
			cl.delxs.push_back(0.0);
			cl.delys.push_back(0.0);
			return cl;
		}
		double vali = pts.db[h.i].ys[y_pos] - crossing_y;
		double valinext = pts.db[h.inext].ys[y_pos] - crossing_y;
		unsigned int j = h.i;
		if (fabs(valinext) < fabs(vali))
			j = h.inext;
		double delx = h.xnext - h.x;
		double dely = pts.db[j].ys[y_pos];
		sln = pts.db[j];
		cl.delxs.push_back(delx);
		cl.delys.push_back(dely);
		return cl;
	}

	// no crossing is found, so find the minimum value
	double absValMin = fabs(vals[0]), tmp;
	unsigned int index_absValMin = 0;
	for (unsigned int ii = 0; ii < szPos; ++ii)
	{
		tmp = fabs(vals[ii]);
		if (tmp < absValMin)
		{
			absValMin = tmp;
			index_absValMin = ii;
		}
	}
	double delx = 0.0;
	double dely = vals[index_absValMin];
	if (index_absValMin == 0)
	{
		cl.convType = ct1d_xmin_lim;
		if (szPos > 1)
			delx = pts.db[poss[1]].x - pts.db[poss[0]].x;
	}
	else if (index_absValMin == (szPos - 1))
	{
		cl.convType = ct1d_xmax_lim;
		if (szPos > 1)
			delx = pts.db[poss[szPos - 1]].x - pts.db[poss[szPos - 2]].x;
	}
	else
	{
		cl.convType = ct1d_notCrossing;
		double delxm = pts.db[poss[index_absValMin]].x - pts.db[poss[index_absValMin - 1]].x;
		double delxp = pts.db[poss[index_absValMin + 1]].x - pts.db[poss[index_absValMin]].x;
		delx = MAX(delxm, delxp);
	}
	cl.delxs.push_back(delx);
	cl.delys.push_back(dely);
	sln = pts.db[poss[index_absValMin]];
	return cl;
}

ConvergenceLog Solver1D::Solve_x4_Extremum(genIndexVal& sln, bool isMin, unsigned int y_pos, double tol_y)
{
	ConvergenceLog cl;
	if (tol_y < 0)
	{
		if (tol_ys.size() > y_pos)
			tol_y = tol_ys[y_pos];
		else
			tol_y = tol_x;
	}
	double tol_x_denom = tol_x * tol_x;
	tol_x_denom *= (3.0 * tol_x_denom);
	unsigned int cntrAddedOnSides = 0;
	int jExtremum = -1, jExtremumm1 = -1, jExtremump1 = -1;
	while (jExtremum < 0)
	{
		bool pointAddedOrExtremum = Solve_x4_Extremum_Aux(cntrAddedOnSides++, isMin, y_pos, jExtremum, jExtremumm1, jExtremump1, cl);
		if (cl.convType == ct1d_nanSln)
		{
			sln = pts.db[0];
			sln.MakeNaN();
			return cl;
		}
		if ((cl.convType == ct1d_xmin_lim) || (cl.convType == ct1d_xmax_lim))
		{
			sln = pts.db[jExtremum];
			return cl;
		}
		if (!pointAddedOrExtremum && (jExtremum < 0))
			return cl;
	}

	double fact = 1.0;
	if (!isMin)
		fact = -1.0;

	double vL, vR, vM, xL, xR, xM;
	xL = pts.db[jExtremumm1].x;
	xM = pts.db[jExtremum].x;
	xR = pts.db[jExtremump1].x;

	vL = fact * pts.db[jExtremumm1].ys[y_pos];
	vM = fact * pts.db[jExtremum].ys[y_pos];
	vR = fact * pts.db[jExtremump1].ys[y_pos];

	cl.iterCnt = 0;
	bool checkIfPtExists = read_Results;

	double prev_xNew = 0.5 * (xL + xR);
	bool delx_good = true, dely_good = true;
	// xij = xj - xi, same for v
	double xNew, vNew = 0.0, delx, dely;
	double xLM, xLR, xMR;
	// use Largrange (FEM shape function) 2nd order interpolation between 3 points and find slope = 0 for extremum vlaue
	double denomL, denomM, denomR, dLMvR, dLRvM, dMRvL, denom;
	while (cl.iterCnt < maxNumIter)
	{
		xLM = xM - xL; xLR = xR - xL, xMR = xR - xM;
		denomL = xLM * xLR;	//	xML * xRL;
		denomM = -xLM * xMR;	//	xLM * xRM;
		denomR = xLR * xMR;
		dLMvR = denomL * denomM * vR;
		dLRvM = denomL * denomR * vM;
		dMRvL = denomM * denomR * vL;

		// v = vL (x - xM)(x - xR)/dL + vM (x - xL)(x - xR)/dM + vR (x - xL)(x - xM)/dR
		// derivative zero gives 
		// v' = vL (2x - (xM + xR))/dL + vM (2x - (xL + xR))/dM + vR (2x - (xL + xM))/dR = 0 -> gives xNew
		denom = dLMvR + dLRvM + dMRvL; // if zero it means all values are equal
		if (fabs(denom) < tol_x_denom)
		{
			sln = pts.db[jExtremum];
			pts.Sort(true, tol_x_eq_check);
			return cl;
		}
		xNew = 0.5 *  (dMRvL * (xM + xR) + dLRvM * (xL + xR) + dLMvR * (xL + xM)) / denom;
		if ((xL - xNew > tol_x_eq_check) || (xNew - xR > tol_x_eq_check))
		{
			cout << "xL\t" << xL << "\tvL\t" << vL << '\n';
			cout << "xM\t" << xM << "\tvM\t" << vM << '\n';
			cout << "xR\t" << xR << "\tvR\t" << vR << '\n';
			cout << "xNew\t" << xNew << "\tvNew\t" << vNew << '\n';
			THROW("Extreme point falls outside the interval, some possible error in calculations (can be finite precision)\n");
		}
		if (cl.iterCnt == 0)
		{
			if (fabs(xNew - xL) < tol_x_eq_check)
			{
				sln = pts.db[jExtremumm1];
				pts.Sort(true, tol_x_eq_check);
				return cl;
			}
			if (fabs(xNew - xR) < tol_x_eq_check)
			{
				sln = pts.db[jExtremump1];
				pts.Sort(true, tol_x_eq_check);
				return cl;
			}
			// now the point is "inside" the interval
			// check if it's within the tolerance range of the middle point
			if (fabs(xNew - xM) < tol_x_eq_check)
			{
				sln = pts.db[jExtremum];
				pts.Sort(true, tol_x_eq_check);
				return cl;
			}
		}
		// point "inside" and far from the middle point
		sln.ys.clear();
		sln.x = xNew;
		sln.index_main = 3;
//		sln.index_sec = cl.iterCnt;
		genIndexVal_DBH retVal = Compute_Add_pt_value(sln, checkIfPtExists);
		if (isnan(sln.ys[y_pos]))
		{
			cl.convType = ct1d_nanSln;
			pts.Sort(true, tol_x_eq_check);
			return cl;
		}
		if (retVal.pos < 0)
		{
			cl.convType = ct1d_failedSln;
			pts.Sort(true, tol_x_eq_check);
			return cl;
		}
		vNew = fact * sln.ys[y_pos];
		// if vNew is too close to vM, vNew is returned as the solution
		if (fabs(vNew - vM) < tol_y)
		{
			pts.Sort(true, tol_x_eq_check);
			return cl;
		}
		// two cases: whether vNew is smaller than vM or not
		delx = xNew - xM;
		dely = vNew - vM;
		delx_good = (fabs(delx) < tol_x);
		dely_good = (fabs(dely) < tol_y);

		if (vNew < vM) // vNew becomes the new minimum value
		{
			cl.iterExtremumNotChanged = 0;
			if (xNew < xM)	// on the left  side - L unchanged, New -> M, M -> R
			{
				xR = xM; vR = vM;
			}
			else			// on the right side - R unchanged, New -> M, M -> L
			{
				xL = xM; vL = vM;
			}
			xM = xNew; vM = vNew;
			cl.delxs.push_back(delx);
			cl.delys.push_back(dely);

			if (delx_good && dely_good)
			{
				pts.Sort(true, tol_x_eq_check);
				return cl;
			}
		}
		else // vM < vNew -> so M stays as the middle point, but L or R point change
		{

			if (++cl.iterExtremumNotChanged == maxNumIterExtremumNotChanging)
			{
				cl.delxs.push_back(0.0);
				cl.delys.push_back(0.0);
				pts.Sort(true, tol_x_eq_check);
				return cl;
			}
			// now shifting other points
			// delx, dely are really zero, 
			// but for convergence check purposes, they are calculated based on the change of L or R point
			if (xNew < xM)	// on the left  side - M, R unchanged, New -> L 
			{
				xL = xNew; vL = vNew;
			}
			else			// on the right side - M, L unchanged, New -> R 
			{
				xR = xNew; vR = vNew;
			}
			cl.delxs.push_back(delx);
			cl.delys.push_back(dely);

			if (delx_good && dely_good)
			{
				pts.Sort(true, tol_x_eq_check);
				return cl;
			}
		}
		++cl.iterCnt;
	}
	pts.Sort(true, tol_x_eq_check);
	if (cl.iterCnt >= maxNumIter)
	{
		if (!delx_good)
		{
			cl.convType = ct1d_large_delx;
			if ((!dely_good) && ((fabs(dely) / tol_y) > (delx / tol_x)))
				cl.convType = ct1d_large_dely;
		}
		else
		{
			cl.convType = ct1d_large_dely;
			if (dely_good)
			{
				THROW("Not both convergence measures can be good");
			}
		}
		return cl;
	}
	return cl;
}

bool Solver1D::Solve_x4_Extremum_Aux(unsigned int cntrAddedOnSides, bool isMin, unsigned int y_pos, int& jExtremum, int& jExtremumm1, int& jExtremump1, ConvergenceLog& cl)
{
	pts.Sort(true, tol_x_eq_check);
	vector<int> poss;
	unsigned int szPos = getValidValues(y_pos, poss);
	if (szPos < 2)
	{
		cl.convType = ct1d_nanSln;
		return false;
//		THROW("not enough points\n");
	}
	double fact = 1.0;
	if (!isMin)
		fact = -1.0;

	vector<double> vals(szPos);
	unsigned int posL = poss[0], posR = poss[szPos - 1];
	
	double tmp, minVal = fact * pts.db[posL].ys[y_pos];
	jExtremum = 0;
	int indexLast = posR;
	unsigned int jjExtremum = 0;
	for (unsigned int ii = 0; ii < szPos; ++ii)
	{
		unsigned int i = poss[ii];
		tmp = fact * pts.db[i].ys[y_pos];
		vals[ii] = tmp;
		if (tmp < minVal)
		{
			jjExtremum = ii;
			jExtremum = i;
			minVal = tmp;
		}
	}
	bool extremumAtL = (jjExtremum == 0);
	bool extremumAtR = (jjExtremum == (szPos - 1));
	if (!extremumAtL && !extremumAtR)
	{
		jExtremumm1 = poss[jjExtremum - 1];
		jExtremump1 = poss[jjExtremum + 1];
		return true;
	}

	jExtremum = -1;
	jExtremumm1 = -1;
	jExtremump1 = -1;
	double valL = fact * vals[0];
	double valR = fact * vals[szPos - 1];
	double xL = pts.db[posL].x;
	double xR = pts.db[posR].x;
	if (del_secondary_x < 0.0)
	{
		cout << "del_secondary_x\t" << del_secondary_x << '\n';
		THROW("Need to have at least 2 secondary points or directly provide del_secondary_x\n");
	}
	double xLProposed = xL - del_secondary_x;
	double xRProposed = xR + del_secondary_x;
	bool canDecrease_xL = ((xLProposed >= xMin) && (posL == 0));
	bool canIncrease_xR = ((xRProposed <= xMax) && (posR == (pts.db.size() - 1)));
	int side = -1; // 0 -> L, 1 -> R
	if (extremumAtL)
	{
		if (canDecrease_xL)
			side = 0;
		else if (canIncrease_xR)
			side = 1;
	}
	else if (extremumAtR)
	{
		if (canIncrease_xR)
			side = 1;
	}
	if (side == -1)
	{
		if (valR < valL)
		{
			jExtremum = posR;
			cl.convType = ct1d_xmax_lim;
		}
		else
		{
			jExtremum = posL;
			cl.convType = ct1d_xmin_lim;
		}
		return false;
	}
	genIndexVal giv;
	giv.index_main = 2;
//	giv.index_sec = cntrAddedOnSides;
	if (side == 0)
		giv.x = xLProposed;
	else
		giv.x = xRProposed;
	bool checkIfPtExists = true;
	genIndexVal_DBH retVal = Compute_Add_pt_value(giv, checkIfPtExists);
	if (isnan(giv.ys[y_pos]))
	{
		cl.convType = ct1d_nanSln;
		return false;
	}
	if (retVal.pos < 0)
	{
		cl.convType = ct1d_failedSln;
		return false;
	}
	return true;
}

ConvergenceLog Solver1D::Solve_x4_Extremum_NoPointAdded(genIndexVal& sln, bool isMin, unsigned int y_pos)
{
	ConvergenceLog cl;
	double fact = 1.0;
	if (!isMin)
		fact = -1.0;

	pts.Sort(true, tol_x_eq_check);
	vector<int> poss;
	unsigned int szPos = getValidValues(y_pos, poss);
	unsigned int posL = poss[0], posR = poss[szPos - 1];
	vector<double> vals(szPos);

	double valMin = fact * pts.db[posL].ys[y_pos];
	unsigned int index_ValMin, iindex_ValMin = 0;

	for (unsigned int ii = 0; ii < szPos; ++ii)
	{
		unsigned int i = poss[ii];
		double val = fact * pts.db[i].ys[y_pos];
		vals[ii] = val;
		if (val < valMin)
		{
			valMin = val;
			iindex_ValMin = ii;
		}
	}
	double delx = 0.0, dely = 0.0;
	index_ValMin = poss[iindex_ValMin];
	sln = pts.db[index_ValMin];

	if (iindex_ValMin == 0)
	{
		cl.convType = ct1d_xmin_lim;
		if (szPos > 1)
		{
			delx = pts.db[poss[1]].x - pts.db[poss[0]].x;
			dely = fabs(pts.db[poss[1]].ys[y_pos] - pts.db[poss[0]].ys[y_pos]);
		}
	}
	else if (iindex_ValMin == (szPos - 1))
	{
		cl.convType = ct1d_xmax_lim;
		if (szPos > 1)
		{
			unsigned int indr = posR, indl = poss[iindex_ValMin - 1];
			delx = pts.db[indr].x - pts.db[indl].x;
			dely = fabs(pts.db[indr].ys[y_pos] - pts.db[indl].ys[y_pos]);
		}
	}
	else // in the middle
	{
		cl.convType = ct1d_yes;
		unsigned int indr = index_ValMin, indl = poss[iindex_ValMin - 1];
		double delxm = pts.db[indr].x - pts.db[indl].x;
		double delym = fact * (pts.db[indr].ys[y_pos] - pts.db[indl].ys[y_pos]);
		indl = indr;
		indr = poss[iindex_ValMin + 1];
		double delxp = pts.db[indr].x - pts.db[indl].x;
		double delyp = pts.db[indr].ys[y_pos] - pts.db[indl].ys[y_pos];
		delx = MAX(delxm, delxp);
		dely = MAX(fabs(delym), fabs(delyp));
	}
	cl.delxs.push_back(delx);
	cl.delys.push_back(dely);
	return cl;
}

unsigned int Solver1D::getValidValues(unsigned int y_pos, vector<int>& poss)
{
	unsigned int sz = pts.size();
	poss.clear();
	double val;
	for (unsigned int i = 0; i < sz; ++i)
	{
		val = pts.db[i].ys[y_pos];
		if (!isnan(val))
			poss.push_back(i);
	}
	return poss.size();
}

void Solver1D::Store_pts(bool print_dbAsIs)
{
	fstream out(name_log.c_str(), ios::out);
	out << setprecision(22);
	pts.Write_genIndexVal_DB(out, print_dbAsIs, equality_by_x, tol_x_eq_check);
}

void Solver1D::Store_pts_unsorted_SlnEnd(int y_pos)
{
	string ser = "";
	if (y_pos >= 0)
	{
		toString(y_pos, ser);
		ser = "_ypos" + ser;
	}
	string fileName = name_log_x_unsorted + ser + ".log";
	fstream out(fileName.c_str(), ios::out);
	out << setprecision(22);
	pts.Write_genIndexVal_DB(out, false, false, tol_x_eq_check);
}

void Solver1D::Restore_pts(bool read_dbAsIs, bool compute_ys_not_computed)
{
	fstream in(name_log.c_str(), ios::in);
	if (!in.is_open())
		return;
	pts.Read_genIndexVal_DB(in, read_dbAsIs, equality_by_x, tol_x_eq_check);
	unsigned int sz = pts.size();
	for (unsigned int i = 0; i < sz; ++i)
	{
		unsigned int index_main = pts.db[i].index_main;
		unsigned int index_sec = pts.db[i].index_sec;
		startInds[index_main] = MAX((unsigned int)startInds[index_main], index_sec + 1);
	}
	if (!compute_ys_not_computed)
		return;
	for (unsigned int i = 0; i < sz; ++i)
	{
		if (pts.db[i].getSize_y() == 0)
			Compute_Add_pt_value(pts.db[i], true);
	}
}

genIndexVal_DBH Solver1D::Compute_Add_pt_value(genIndexVal& giv, bool checkIfPtExists)
{
	genIndexVal_DBH retVal = pts.Add_Pt(giv, checkIfPtExists, equality_by_x, tol_x_eq_check);
	unsigned int sz_y = retVal.y_size;
	if (sz_y > 0)
		return retVal;
	if (retVal.pos < 0)
	{
		out_error << "giv_not added\n" << giv << '\n';
		return retVal;
	}
	if (!retVal.found)
	{
		unsigned int index_main = giv.index_main;
//		unsigned int index_sec = giv.index_sec;
		unsigned int index_sec = startInds[index_main];
		++startInds[index_main];
		giv.index_sec = index_sec;
		pts.db[retVal.pos].index_sec = index_sec;
	}
	bool validSln = function->ComputeValue(pts.db[retVal.pos]);
	if (validSln == false)
	{
		out_error << "pts.db[retVal.pos]\n" << pts.db[retVal.pos] << '\n';
		// out_error << "value cannot be computed\n";
//		THROW("Value cannot be computed\n");
		retVal.pos = -1, retVal.y_size = -1;
		return retVal;
	}
	bool print_dbAsIs = false;
	Store_pts(print_dbAsIs);
	giv = pts.db[retVal.pos];
	return retVal;
}

void Solver1D::MAIN_ProcessConfigFile(gFx2y* functionIn, string confName)
{
	if (confName != "none")
	{
		fstream in(confName.c_str(), ios::in);
		if (!in.is_open())
		{
			cout << confName << '\n';
			THROW("Cannot open file\n");
		}
		Read_Solver1D(in);
	}
	SolveAllSolveMode_AfterReadingConfig(functionIn);
}

ConvergenceLog Solver1D::SolveYCrossingOrExtremum(genIndexVal& sln, gFx2y* functionIn, unsigned int y_pos, bool isExtremum, bool isMin, double crossing_y, bool addAdditionalPoints, double x0, double tol_y, bool doInitialization)
{
	if (isExtremum)
		return SolveYExtremum(functionIn, sln, isMin, y_pos, addAdditionalPoints, tol_y, doInitialization);
	return SolveYCrossing(functionIn, sln, crossing_y, y_pos, addAdditionalPoints, x0, tol_y, doInitialization);
}

ConvergenceLog Solver1D::SolveYCrossing(gFx2y* functionIn, genIndexVal& sln, double crossing_y, unsigned int y_pos, bool addAdditionalPoints, double x0, double tol_y, bool doInitialization)
{
	if (doInitialization)
		InitializeFunction(functionIn);
	ConvergenceLog cl;
	if (addAdditionalPoints)
		cl = Solve_x4_Crossing(sln, crossing_y, y_pos, x0, tol_y);
	else
		cl = Solve_x4_Crossing_NoPointAdded(sln, crossing_y, y_pos, x0, tol_y);
	if (doInitialization)
		Store_pts_unsorted_SlnEnd(y_pos);
	return cl;
}

ConvergenceLog Solver1D::SolveYExtremum(gFx2y* functionIn, genIndexVal& sln, bool isMin, unsigned int y_pos, bool addAdditionalPoints, double tol_y, bool doInitialization)
{
	if (doInitialization)
		InitializeFunction(functionIn);
	ConvergenceLog cl;
	if (addAdditionalPoints)
		cl = Solve_x4_Extremum(sln, isMin, y_pos, tol_y);
	else
		cl = Solve_x4_Extremum_NoPointAdded(sln, isMin, y_pos);
	if (doInitialization)
		Store_pts_unsorted_SlnEnd(y_pos);
	return cl;
}

unsigned int Solver1D::SolveYCrossings(gFx2y* functionIn, vector<genIndexVal>& slns, vector<ConvergenceLog>& cls, double crossing_y, unsigned int y_pos, bool addAdditionalPoints, double tol_y, bool doInitialization)
{
	if (doInitialization)
		InitializeFunction(functionIn);
	if (tol_y < 0)
	{
		if (tol_ys.size() > y_pos)
			tol_y = tol_ys[y_pos];
		else
			tol_y = tol_x;
	}
	double tol_y_zero = MIN(0.01 * tol_y, 0.01 * tol_x);
	vector<int> poss;
	vector<sgnT> signs;
	vector<cross_optimal_helper> helpers;
	unsigned int szPos = Cross_Helper(crossing_y, y_pos, tol_y_zero, poss, signs, helpers);
	if (szPos < 2)
		return 0;
	vector<double> x0s;
	unsigned numXs = helpers.size();
	for (unsigned int i = 0; i < numXs; ++i)
	{
		if (helpers[i].signi == sgn_0)
			x0s.push_back(helpers[i].x);
		else
			x0s.push_back(0.5 * (helpers[i].x + helpers[i].xnext));
	}
	slns.resize(numXs);
	cls.resize(numXs);
	for (unsigned int i = 0; i < numXs; ++i)
		cls[i] = SolveYCrossing(functionIn, slns[i], crossing_y, y_pos, addAdditionalPoints, x0s[i], tol_y, false);
	return numXs;
}

void Solver1D::Initialize_AfterReading()
{
	sz_indices4paras1D = paras1D.size();
	sz_posConfs2Solve = posConfs2Solve.size();
	for (unsigned int i = 0; i < sz_posConfs2Solve; ++i)
		posConfs2Solve[i].Initialize_Solver1D_1posConf();

	if (do_posConfs_AddPtSolve == false)
	{
		do_posConfs_second_notAddPtSolve = false;
		do_posConfs_first_notAddPtSolve = true;
	}
	sz_solve_mode = 0;
	solve_mode_addAdditionalPoints.clear();
	solve_mode_counter.clear();

	unsigned int cntr = 0;
	bool hasBeforeAfter = false;
	if (do_posConfs_first_notAddPtSolve)
	{
		++sz_solve_mode;
		solve_mode_counter.push_back(cntr++);
		if (do_posConfs_second_notAddPtSolve)
			hasBeforeAfter = true;
		solve_mode_addAdditionalPoints.push_back(false);
	}
	if (do_posConfs_AddPtSolve)
	{
		++sz_solve_mode;
		solve_mode_addAdditionalPoints.push_back(true);
		solve_mode_counter.push_back(cntr++);
	}
	if (do_posConfs_second_notAddPtSolve)
	{
		++sz_solve_mode;
		solve_mode_addAdditionalPoints.push_back(false);
		solve_mode_counter.push_back(cntr++);
	}
}

void Solver1D::SolveOneSolveMode_AfterInitialization(vector<genIndexVal>& slns, vector<ConvergenceLog>& cls, gFx2y* functionIn, bool addAdditionalPoints, int mode_counter)
{
	slns.clear();
	cls.clear();
	string addedName = "__SlnCntr_";
	string ser;
	toString(mode_counter, ser);
	addedName += ser;
	addedName += "_AddedPts";
	toString(addAdditionalPoints, ser);
	addedName += ser;
	string fileNameAllSlns = name_WOExt + addedName;
	string fileName_pt_x_sorted = fileNameAllSlns + "_pts_x_sorted.out";
	string fileName_pt_i_sorted = fileNameAllSlns + "_pts_i_sorted.out";
	fileNameAllSlns += ".txt";//"_slns.txt";

	fstream outSln(fileNameAllSlns.c_str(), ios::out);
	PrintSln_Header(outSln);

	for (unsigned int pci = 0; pci < sz_posConfs2Solve; ++pci)
	{
		cout << "ST: mode_counter\t" << mode_counter << "\tpci\t" << pci << "\ttsz_posConfs2Solve\t" << sz_posConfs2Solve << '\n';
		string fileName = baseName + addedName + "_solveNum_", ser;
		toString(pci, ser);
		fileName += ser;
		Solver1D_1posConf spc = posConfs2Solve[pci];
		for (unsigned int j = 0; j < spc.sz; ++j)
		{
			fstream out_pt_x(fileName_pt_x_sorted.c_str(), ios::out);
			fstream out_pt_i(fileName_pt_i_sorted.c_str(), ios::out);

			string nm = fileName;
			if (spc.sz > 1)
			{
				nm += "_subSolveIndex_";
				string ser2;
				toString(j, ser2);
				nm += ser2;
			}
			nm += ".txt";

			fstream in(nm.c_str(), ios::in);
			fstream out;
			if (!in.is_open())
			{
				out.open(nm.c_str(), ios::app);
				PrintSln_Header(out);
			}
			else
			{
				in.close();
				out.open(nm.c_str(), ios::app);
			}
			genIndexVal sln;
			ConvergenceLog cl = SolveYCrossingOrExtremum(sln, functionIn, spc.pos_y, spc.vec_isExtremum[j], spc.isMins[j], spc.crossing_ys[j], addAdditionalPoints, spc.x0s[j], spc.tol_ys[j], false);
			slns.push_back(sln);
			cls.push_back(cl);

			PrintSln_Values(out, sln, cl, spc.vec_isExtremum[j]);
			PrintSln_Values(outSln, sln, cl, spc.vec_isExtremum[j]);

			out_pt_x << setprecision(22);
			pts.Write_genIndexVal_DB_WithHeader(functionIn, out_pt_x, false, true, tol_x_eq_check);
			out_pt_i << setprecision(22);
			pts.Write_genIndexVal_DB_WithHeader(functionIn, out_pt_i, false, false, tol_x_eq_check);
		}
		cout << "EN: mode_counter\t" << mode_counter << "\tpci\t" << pci << "\ttsz_posConfs2Solve\t" << sz_posConfs2Solve << '\n';
	}
}

void Solver1D::SolveAllSolveMode_AfterReadingConfig(gFx2y* functionIn)
{
	Initialize_AfterReading();
	InitializeFunction(functionIn);
	for (unsigned int i = 0; i < sz_solve_mode; ++i)
	{
		cout << "ST slnMode\t" << i << "\tsz_solve_mode\t" << sz_solve_mode << '\n';
		vector<genIndexVal> slns;
		vector<ConvergenceLog> cls;
		SolveOneSolveMode_AfterInitialization(slns, cls, functionIn, solve_mode_addAdditionalPoints[i], solve_mode_counter[i]);
		cout << "EN slnMode\t" << i << "\tsz_solve_mode\t" << sz_solve_mode << '\n';
	}
}

void Solver1D::PrintSln_Header(ostream& out)
{
	for (unsigned int j = 0; j < sz_indices4paras1D; ++j)
		out << "solveIndex" << j << '\t';
	for (unsigned int j = 0; j < sz_indices4paras1D; ++j)
		out << "solvePara" << j << '\t';
	ConvergenceLog cl;
	cl.Print(out, true, false, false, true, 0);
	function->Print_Header(out, num_y);
	out << '\n';
}

void Solver1D::PrintSln_Values(ostream& out, const genIndexVal& sln, ConvergenceLog& cl, bool isExremum)
{
	for (unsigned int k = 0; k < sz_indices4paras1D; ++k)
		out << indices4paras1D[k] << '\t';
	for (unsigned int k = 0; k < sz_indices4paras1D; ++k)
		out << paras1D[k] << '\t';
	cl.Print(out, false, true, false, true, (int)isExremum);
	out << sln.index_main << '\t' << sln.index_sec << '\t' << sln.x;
	for (unsigned int k = 0; k < sln.ys.size(); ++k)
		out << '\t' << sln.ys[k];
	out << '\n';
}

void MAIN_Solver1D(string confName, gFx2y* functionIn)
{
	if (confName == "default")
		confName = "TestFiles/_confSolver1D.txt";
 	bool delFun = false;
	if (functionIn == NULL)
	{
		functionIn = new gFx2y();
		delFun = true;
	}
	Solver1D s1d;
	s1d.MAIN_ProcessConfigFile(functionIn, confName);

	if (delFun)
		delete functionIn;
}

void Test_SolveYCrossing()
{
	bool addAdditionalPoints = true;  //false;
	gFx2y* function = new gFx2y();
	Solver1D s1d;
	genIndexVal sln;
	double crossing_y = 0.0;
	unsigned int y_pos = 1;
	double tol_y = -1.0;
	double x0 = -100.0;
	if (false)
	{
		ConvergenceLog cl = s1d.SolveYCrossing(function, sln, crossing_y, y_pos, x0, tol_y, addAdditionalPoints);
		bool success = (cl.convType == ct1d_yes);
		cout << "success\t" << success;
		cout << "cl\n" << cl << '\n';
		cout << "sln\n";
		sln.Write_genIndexVal(cout);
	}
	else
	{
		vector<genIndexVal> slns;
		vector<ConvergenceLog> cls;
		unsigned int sz = s1d.SolveYCrossings(function, slns, cls, crossing_y, y_pos, addAdditionalPoints, tol_y, true);

		for (unsigned int i = 0; i < sz; ++i)
		{
			cout << "cl" << i << "\n" << cls[i] << '\n';
			cout << "sln" << i << "\n";
			slns[i].Write_genIndexVal(cout);
		}
	}
	delete function;
}

void Test_SolveYExtremum()
{
	bool addAdditionalPoints = false;  //true;
	gFx2y* function = new gFx2y();
	Solver1D s1d;
	genIndexVal sln;
	bool isMin = true;
	unsigned int y_pos = 1;
	double tol_y = -1.0;
	ConvergenceLog cl = s1d.SolveYExtremum(function, sln, isMin, y_pos, tol_y, addAdditionalPoints);
	bool success = (cl.convType == ct1d_yes);
	cout << "success\t" << success;
	cout << "cl\n" << cl << '\n';
	cout << "sln\n";
	sln.Write_genIndexVal(cout);
	delete function;
}
