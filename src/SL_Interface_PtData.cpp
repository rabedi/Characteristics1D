#include "SL_Interface_PtData.h"
#include "SLInterfaceFracturePF.h"
#include "Domain_AllInterfacesAllTimes.h"
#include "SLDescriptorData.h"
#include "globalFunctions.h"


/////////////////
string getName(OneSegmentPFT dat)
{
	if (dat == pft_time)
		return "time";
	if (dat == pft_numCyclesAfterCrackOpening)
		return "numCyclesAfterCrackOpening";
	if (dat == pft_idelu)
		return "delu";
	if (dat == pft_idelv)
		return "delv";
	if (dat == pft_isigma)
		return "sigma";
	if (dat == pft_ichar)
		return "char";
	if (dat == pft_iEne)
		return "iEnePL";
	if (dat == pft_iEner)
		return "Ener";
	if (dat == pft_iEned)
		return "Ened";
	if (dat == pft_ilog10Ened)
		return "log10Ened";
	if (dat == pft_iEned_PL)
		return "Ened_PL";
	if (dat == pft_ilog10Ened_PL)
		return "log10Ened_PL";
	if (dat == pft_ivsolid)
		return "vsolid";
	if (dat == pft_irelusolid)
		return "relusolid";
	if (dat == pft_bepsilon)
		return "bepsilon";
	if (dat == pft_bsigma)
		return "bsigma";
	if (dat == pft_bsigma_maxpw)
		return "bsigma_maxpw";
	if (dat == pft_bK_PL)
		return "K_PL";
	if (dat == pft_bU_PL)
		return "U_PL";
	if (dat == pft_bEneSource_PL)
		return "EneSource_PL";
	if (dat == pft_bEneN_PL)
		return "EneN_PL";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, OneSegmentPFT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= OneSegmentPFT_SIZE)
			THROW("too large of a number\n");
		typeVal = (OneSegmentPFT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < OneSegmentPFT_SIZE; ++i)
	{
		typeVal = (OneSegmentPFT)i; // casting integer to OneSegmentPFT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading OneSegmentPFT\n");
}

//operator for output
ostream& operator<<(ostream& out, OneSegmentPFT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, OneSegmentPFT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}


/////////////////
string getName(PITSS dat)
{
	if (dat == PITSS_none)
		return "none";
	if (dat == pit_sigma0)
		return "sigma0";
	if (dat == pit_sigmaMax)
		return "sigmaMax";
	if (dat == pit_timeDilute)
		return "timeDilute";
	if (dat == pit_timeInteraction)
		return "timeInteraction";
	if (dat == pit_timeSigmaAveMax)
		return "timeSigmaAveMax";
	if (dat == pit_timeSigmaPWMax)
		return "timeSigmaPWMax";
	if (dat == pit_vSolidNegative)
		return "vSolidNegative";
	if (dat == pit_reluSolid1)
		return "reluSolid1";
	if (dat == pit_reluSolid0)
		return "reluSolid0";
	if (dat == pit_reluSolidm1)
		return "reluSolidm1";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, PITSS& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= PITSS_SIZE)
			THROW("too large of a number\n");
		typeVal = (PITSS)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < PITSS_SIZE; ++i)
	{
		typeVal = (PITSS)i; // casting integer to PITSS, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading PITSS\n");
}

//operator for output
ostream& operator<<(ostream& out, PITSS dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, PITSS& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

/////////////////
string getName(FragmentSizeT dat)
{
	if (dat == fst_dilute)
		return "dilute";
	if (dat == fst_dilute_approx)
		return "dilute_approx";
	if (dat == fst_zhu6a)
		return "zhu6a";
	if (dat == fst_zhu6b)
		return "zhu6b";
	if (dat == fst_Glenn)
		return "Glenn";
	if (dat == fst_Grady)
		return "Grady";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, FragmentSizeT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= FragmentSizeT_SIZE)
			THROW("too large of a number\n");
		typeVal = (FragmentSizeT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < FragmentSizeT_SIZE; ++i)
	{
		typeVal = (FragmentSizeT)i; // casting integer to FragmentSizeT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading FragmentSizeT\n");
}

//operator for output
ostream& operator<<(ostream& out, FragmentSizeT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, FragmentSizeT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

SL_interface_Temp_PPtData::SL_interface_Temp_PPtData()
{
	interfacePropPtr = NULL;
	interface_damageLatestValue = 0.0;
	tsrStageLatestValue = TSR_loadingStages_none;
	dn = 0.0;
	dt = 0.0;
	dsep = 0.0;
	tau = 0.1;
	sigmaC = 0.1;
	deltaC = 0.1;

#if RING_PROBLEM
	v_r_latestValue = 0.0;
	v_r_source_latestValue = 0.0;
	sigma_theta_source_latestValue = 0.0;
#endif
}

SL_interface_Temp_PPtData::~SL_interface_Temp_PPtData()
{
	if (interfacePropPtr != NULL)
		delete interfacePropPtr;
}

void SL_interface_Temp_PPtData::SetZero()
{
	damageOffOnMix = sl_interfacial_damage_off;
	contactOffOnMix = sl_contact_off;
	slipOffOnMix = sl_slip_off;
	sl_side_temp_ptData[SDL].SetZero();
	sl_side_temp_ptData[SDR].SetZero();
	if (interfacePropPtr != NULL)
		interfacePropPtr->SetZero();
	interface_damageLatestValue = 0.0;
	tsrStageLatestValue = TSR_loadingStages_none;
	tau = 0.0, dn = 0.0, dt = 0.0, dsep = 0.0;
	sigmaC = 0.0, deltaC = 0.0;
#if RING_PROBLEM
	v_r_latestValue = 0.0;
	v_r_source_latestValue = 0.0;
	sigma_theta_source_latestValue = 0.0;
#endif
}

void SL_interface_Temp_PPtData::Read(istream& in)
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
		if (buf == "damageOffOnMix")
		{
			in >> damageOffOnMix;
		}
		else if (buf == "contactOffOnMix")
		{
			in >> contactOffOnMix;
		}
		else if (buf == "slipOffOnMix")
		{
			in >> slipOffOnMix;
		}
		else if (buf == "sl_side_temp_ptData_L")
		{
			sl_side_temp_ptData[SDL].Read(in);
		}
		else if (buf == "sl_side_temp_ptData_R")
		{
			sl_side_temp_ptData[SDR].Read(in);
		}
		else if (buf == "interfacePropPtr")
		{
			if (interfacePropPtr == NULL)
				interfacePropPtr = new SL_Interface_Temp_PtData_IntContFrac();
			interfacePropPtr->Read(in);
		}
		else if (buf == "interface_damageLatestValue")
		{
			READ_NDOUBLE(in, buf, interface_damageLatestValue);
		}
		else if (buf == "tsrStageLatestValue")
		{
			int tmpi;
			READ_NINTEGER(in, buf, tmpi);
			tsrStageLatestValue = (TSR_loadingStages)tmpi;
		}
		else if (buf == "tau")
		{
			READ_NDOUBLE(in, buf, tau);
		}
		else if (buf == "dn")
		{
			READ_NDOUBLE(in, buf, dn);
		}
		else if (buf == "dt")
		{
			READ_NDOUBLE(in, buf, dt);
		}
		else if (buf == "dsep")
		{
			READ_NDOUBLE(in, buf, dsep);
		}
		else if (buf == "sigmaC")
		{
			READ_NDOUBLE(in, buf, sigmaC);
		}
		else if (buf == "deltaC")
		{
			READ_NDOUBLE(in, buf, deltaC);
		}
#if RING_PROBLEM
		else if (buf == "v_r_latestValue")
		{
			READ_NDOUBLE(in, buf, v_r_latestValue);
		}
		else if (buf == "v_r_source_latestValue")
		{
			READ_NDOUBLE(in, buf, v_r_source_latestValue);
		}
		else if (buf == "sigma_theta_source_latestValue")
		{
			READ_NDOUBLE(in, buf, sigma_theta_source_latestValue);
		}
#endif
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void SL_interface_Temp_PPtData::Output_ScalarValues(ostream & out, IOF_type iot, double space_or_time)
{
	if (interfacePropPtr != NULL)
	{
		Output_ScalarValuesAux(out, iot, space_or_time, interfacePropPtr->interface_scalar_Vals, damageOffOnMix, contactOffOnMix, slipOffOnMix);
		return;
	}
	Output_ScalarValuesAux(out, iot, space_or_time, NULL, damageOffOnMix, contactOffOnMix, slipOffOnMix);
}

void SL_interface_Temp_PPtData::MakeReady_For_ContactFractureRuns(double beta_delU, double beta_traction, const VEC& sigmaI, double delu0Change)
{
	if (interfacePropPtr == NULL)
		interfacePropPtr = new SL_Interface_Temp_PtData_IntContFrac();
	VEC delU;
	for (int i = 0; i < DiM; ++i)
		delU[i] = sl_side_temp_ptData[SDR].u_downstream_latestValue[i] - sl_side_temp_ptData[SDL].u_downstream_latestValue[i];
	delU[0] += delu0Change;

	interfacePropPtr->del_u_nt_parts.nt_processed = false;
	interfacePropPtr->del_u_nt_parts.Compute_derivedValues_From_vec(delU, beta_delU);

	interfacePropPtr->sigma_I_nt_parts.nt_processed = false;
	interfacePropPtr->sigma_I_nt_parts.Compute_derivedValues_From_vec(sigmaI, beta_traction);

	// velocity term value is not computed as it may not be needed
	interfacePropPtr->del_v_nt_parts.nt_processed = false;
}


void SL_Interface_PtData_OneSide::Read(istream& in)
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
		if (buf == "v_downstream_final")
		{
			ReadV(v_downstream_final, in);
		}
		else if (buf == "sigma_downstream_final")
		{
			ReadV(sigma_downstream_final, in);
		}
		else if (buf == "u_downstream_final")
		{
			ReadV(u_downstream_final, in);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void SL_Interface_PtData_OneSide::SetZero()
{
	setValue(v_downstream_final, 0.0);
	setValue(sigma_downstream_final, 0.0);
	setValue(u_downstream_final, 0.0);
}

SL_interfacePPtData::SL_interfacePPtData()
{
	interface_damage_final = 0.0;
	interface_damage_source_final = 0.0;
	interface_time = 0.0;
	maxEffDelU = 0.0;
	tsrStage = TSR_loadingStages_none;
}

void SL_interfacePPtData::Read(istream& in)
{
	sl_side_ptData[SDL].SetZero();
	sl_side_ptData[SDR].SetZero();
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
		if (buf == "sl_side_ptData_L")
		{
			sl_side_ptData[SDL].Read(in);
		}
		else if (buf == "sl_side_ptData_R")
		{
			sl_side_ptData[SDR].Read(in);
		}
		else if (buf == "interface_damage_final")
		{
			READ_NDOUBLE(in, buf, interface_damage_final);
		}
		else if (buf == "interface_damage_source_final")
		{
			READ_NDOUBLE(in, buf, interface_damage_source_final);
		}
		else if (buf == "maxEffDelU")
		{
			READ_NDOUBLE(in, buf, maxEffDelU);
		}
		else if (buf == "tsrStage")
		{
			int tmpi;
			READ_NINTEGER(in, buf, tmpi);
			tsrStage = (TSR_loadingStages)tmpi;
		}
		else if (buf == "interface_time")
		{
			READ_NDOUBLE(in, buf, interface_time);
		}
#if RING_PROBLEM
		else if (buf == "v_r_final")
		{
			READ_NDOUBLE(in, buf, v_r_final);
		}
		else if (buf == "v_r_source_final")
		{
			READ_NDOUBLE(in, buf, v_r_source_final);
		}
		else if (buf == "sigma_theta_source_final")
		{
			READ_NDOUBLE(in, buf, sigma_theta_source_final);
		}
#endif
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void SL_interfacePPtData::Output_FinalSolution(ostream& out, IOF_type iot, double space_or_time, double x, double t, bool has_ring_opened1D_al)
{
	// time
	if (iot == iof_ascii)
	{
		out << setprecision(22);
		out << space_or_time << '\t';
		// stress
		for (int i = 0; i < DiM; ++i)
			out << sl_side_ptData[SDL].sigma_downstream_final[i] << '\t';
		for (int i = 0; i < DiM; ++i)
			out << sl_side_ptData[SDR].sigma_downstream_final[i] << '\t';
		if ((g_domain != NULL) && (!g_domain->b_ring_opened1D))
		{
			// velocities
			for (int i = 0; i < DiM; ++i)
				out << sl_side_ptData[SDL].v_downstream_final[i] << '\t';
			for (int i = 0; i < DiM; ++i)
				out << sl_side_ptData[SDR].v_downstream_final[i] << '\t';
			// displacements
			for (int i = 0; i < DiM; ++i)
				out << sl_side_ptData[SDL].u_downstream_final[i] << '\t';
			for (int i = 0; i < DiM; ++i)
				out << sl_side_ptData[SDR].u_downstream_final[i] << '\t';
		}
		else
		{
			double ax = g_SL_desc_data.load_parameters[0] * x, axt = ax * t;
			unsigned int i = 0;
			out << (sl_side_ptData[SDL].v_downstream_final[i] - ax) << '\t';
			if (!has_ring_opened1D_al)
				out << (sl_side_ptData[SDR].v_downstream_final[i] - ax) << '\t';
			else
				out << (sl_side_ptData[SDR].v_downstream_final[i] + ax) << '\t';
			out << (sl_side_ptData[SDL].u_downstream_final[i] - axt) << '\t';
			if (!has_ring_opened1D_al)
				out << (sl_side_ptData[SDR].u_downstream_final[i] - axt) << '\t';
			else
				out << (sl_side_ptData[SDR].u_downstream_final[i] + axt) << '\t';
		}
		out << interface_damage_final << '\t' << interface_damage_source_final << '\t';
		out << maxEffDelU;
#if TSR_STAGE_IO
		out << '\t' << tsrStage;
#endif
#if RING_PROBLEM
		out << '\t' << v_r_final << '\t' << v_r_source_final << '\t' << sigma_theta_source_final;
#endif
		out << '\n';
	}
	else if (iot == iof_binary)
	{
		out.write((char*)&space_or_time, sizeof(space_or_time));
		for (int i = 0; i < DiM; ++i)
			out.write((char*)&sl_side_ptData[SDL].sigma_downstream_final[i], sizeof(double));
		for (int i = 0; i < DiM; ++i)
			out.write((char*)&sl_side_ptData[SDR].sigma_downstream_final[i], sizeof(double));
		if (!g_domain->b_ring_opened1D)
		{
			// velocities
			for (int i = 0; i < DiM; ++i)
				out.write((char*)&sl_side_ptData[SDL].v_downstream_final[i], sizeof(double));
			for (int i = 0; i < DiM; ++i)
				out.write((char*)&sl_side_ptData[SDR].v_downstream_final[i], sizeof(double));
			// displacements
			for (int i = 0; i < DiM; ++i)
				out.write((char*)&sl_side_ptData[SDL].u_downstream_final[i], sizeof(double));
			for (int i = 0; i < DiM; ++i)
				out.write((char*)&sl_side_ptData[SDR].u_downstream_final[i], sizeof(double));
		}
		else
		{
			double ax = g_SL_desc_data.load_parameters[0] * x, axt = ax * t;
			unsigned int i = 0;
			double tmp = sl_side_ptData[SDL].v_downstream_final[i] - ax;
			out.write((char*)&tmp, sizeof(double));
			tmp = sl_side_ptData[SDR].v_downstream_final[i] - ax;
			out.write((char*)&tmp, sizeof(double));
			tmp = sl_side_ptData[SDL].u_downstream_final[i] - axt;
			out.write((char*)&tmp, sizeof(double));
			tmp = sl_side_ptData[SDR].u_downstream_final[i] - axt;
			out.write((char*)&tmp, sizeof(double));
		}
		out.write((char*)&interface_damage_final, sizeof(interface_damage_final));
		out.write((char*)&interface_damage_source_final, sizeof(interface_damage_source_final));
		out.write((char*)&maxEffDelU, sizeof(maxEffDelU));
#if TSR_STAGE_IO
		out.write((char*)&tsrStage, sizeof(tsrStage));
#endif
#if RING_PROBLEM
		out.write((char*)&v_r_final, sizeof(v_r_final));
		out.write((char*)&v_r_source_final, sizeof(v_r_source_final));
		out.write((char*)&sigma_theta_source_final, sizeof(sigma_theta_source_final));
#endif
	}
}

void SL_interfacePPtData::SetZero()
{
	sl_side_ptData[SDL].SetZero();
	sl_side_ptData[SDR].SetZero();
	interface_damage_final = 0.0;
	interface_damage_source_final = 0.0;
	maxEffDelU = 0.0;
	tsrStage = TSR_loadingStages_none;
#if RING_PROBLEM
	v_r_final = 0.0;
	v_r_source_final = 0.0;
	sigma_theta_source_final = 0.0;
#endif
}

void SL_interfacePPtData::get1DValues4Visualization(vector<double>& vals, double delC, InterfaceLocation1DT side, double t, double ax, unsigned int dir, bool periodicBoundaryPt)
{
	double u, v, s;
	double axt = ax * t;
	if (side == ilt_twoSided)
	{
		u = 0.5 * (sl_side_ptData[SDL].u_downstream_final[dir] + sl_side_ptData[SDR].u_downstream_final[dir]) - axt;
		v = 0.5 * (sl_side_ptData[SDL].v_downstream_final[dir] + sl_side_ptData[SDR].v_downstream_final[dir]) - ax;
		s = 0.5 * (sl_side_ptData[SDL].sigma_downstream_final[dir] + sl_side_ptData[SDR].sigma_downstream_final[dir]);
	}
	else if (side == ilt_left)
	{
		u = sl_side_ptData[SDL].u_downstream_final[dir] - axt;
		v = sl_side_ptData[SDL].v_downstream_final[dir] - ax;
		s = sl_side_ptData[SDL].sigma_downstream_final[dir];
	}
	else if (side == ilt_right)
	{
		u = sl_side_ptData[SDR].u_downstream_final[dir] - axt;
		v = sl_side_ptData[SDR].v_downstream_final[dir] - ax;
		s = sl_side_ptData[SDR].sigma_downstream_final[dir];
	}
	if (!g_domain->hasFracture)
	{
		vals.resize(3);
		vals[0] = u;
		vals[1] = v;
		vals[2] = s;
	}
	else
	{
		vals.resize(6);
		vals[0] = u;
		vals[1] = v;
		vals[2] = s;
		if (delC > 0)
		{
			vals[3] = interface_damage_final;
			double delu = sl_side_ptData[SDR].u_downstream_final[dir] - sl_side_ptData[SDL].u_downstream_final[dir];
			if (periodicBoundaryPt)
				delu += g_domain->ring_opened1D_al * t;
			vals[4] = MIN(maxEffDelU / delC, 1.0);
			vals[5] = delu / delC;
		}
		else
		{
			vals[3] = 0.0;
			vals[4] = 0.0;
			vals[5] = 0.0;
		}
	}
}

////////////////
SL_interfacePPtData_Time_Seq::SL_interfacePPtData_Time_Seq()
{
	for (int i = 0; i < SIZE_PT_TIME_SEQUENCE; ++i)
		timeSeqPtrs[i] = NULL;
	curPos = SIZE_PT_TIME_SEQUENCE - 1;
	sz = 0;
}

SL_interfacePPtData_Time_Seq::~SL_interfacePPtData_Time_Seq()
{
	for (int i = 0; i < SIZE_PT_TIME_SEQUENCE; ++i)
		if (timeSeqPtrs[i] != NULL)
			delete timeSeqPtrs[i];
}

SL_interfacePPtData_Time_Seq::SL_interfacePPtData_Time_Seq(const SL_interfacePPtData_Time_Seq& other)
{
	for (int i = 0; i < SIZE_PT_TIME_SEQUENCE; ++i)
		timeSeqPtrs[i] = NULL;
	(*this) = other;

}

SL_interfacePPtData_Time_Seq& SL_interfacePPtData_Time_Seq::operator=(const SL_interfacePPtData_Time_Seq& other)
{
	curPos = other.curPos;
	sz = other.sz;
	for (unsigned int i = 0; i < SIZE_PT_TIME_SEQUENCE; ++i)
	{
		if (timeSeqPtrs[i] != NULL)
			delete timeSeqPtrs[i];
		if (other.timeSeqPtrs[i] != NULL)
		{
			timeSeqPtrs[i] = new SL_interfacePPtData();
			*(timeSeqPtrs[i]) = *(other.timeSeqPtrs[i]);
		}
		else
			timeSeqPtrs[i] = NULL;
	}
	return *this;
}

void SL_interfacePPtData_Time_Seq::Read(istream& in)
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
		if (buf == "curPos")
		{
			READ_NINTEGER(in, buf, curPos);
		}
		if (buf == "sz")
		{
			READ_NINTEGER(in, buf, sz);
		}
		else if (buf == "timeSeqPtrs")
		{
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("istream should start with {");
			}
			READ_NSTRING(in, buf, buf);
			while (buf != "}")
			{
				unsigned int pos;
				if (fromString(buf, pos) == false)
				{
					cout << "buf\t" << buf << '\n';
					THROW("should have point ID\n");
					pos = pos % SIZE_PT_TIME_SEQUENCE;
					sz = MAX(sz, pos + 1);
					if (timeSeqPtrs[pos] != NULL)
						delete timeSeqPtrs[pos];
					timeSeqPtrs[pos] = new SL_interfacePPtData();
					timeSeqPtrs[pos]->Read(in);
				}
			}
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	if (curPos < 0)
		curPos = (sz - 1) % SIZE_PT_TIME_SEQUENCE;
}

void SL_interfacePPtData_Time_Seq::RemoveLastPoint()
{
	delete timeSeqPtrs[curPos];
	timeSeqPtrs[curPos] = NULL;
	--sz;
	if (curPos != 0)
		--curPos;
	else
		curPos = SIZE_PT_TIME_SEQUENCE - 1;
}

SL_interfacePPtData* SL_interfacePPtData_Time_Seq::AddNewPoint(double timeIn, unsigned int& pos)
{
	curPos = (curPos + 1) % SIZE_PT_TIME_SEQUENCE;
	pos = curPos;
	if (timeSeqPtrs[pos] == NULL)
		timeSeqPtrs[pos] = new SL_interfacePPtData();
	else
		timeSeqPtrs[pos]->SetZero();
	++sz;
	timeSeqPtrs[pos]->interface_time = timeIn;
	return timeSeqPtrs[pos];
}

SL_interfacePPtData* SL_interfacePPtData_Time_Seq::GetBackwardPosition(unsigned int relativeBackwardPos)
{
	int pos = curPos - (relativeBackwardPos + 1);
	if (pos < 0)
		pos += SIZE_PT_TIME_SEQUENCE;
	return timeSeqPtrs[pos];
}

SL_interfacePPtData* SL_interfacePPtData_Time_Seq::GetCurrentPosition()
{
	return  timeSeqPtrs[curPos];
}


bool SL_interfacePPtData_Time_Seq::Find_Pts_Around_Time(double timeIn, SL_interfacePPtData*& ptA, double& factor_ptA, SL_interfacePPtData*& ptB, double& factor_ptB)
{
	if (!g_slf_conf->between_steps_adaptivity)
	{
		double relTime = timeIn * g_slf_conf->inv_uniform_del_t;
		int relTimeInt = (int)floor(relTime + g_slf_conf->g_fp_tol);
		if (fabs(relTime - relTimeInt) < g_slf_conf->g_fp_tol)
		{
			ptA = timeSeqPtrs[relTimeInt % SIZE_PT_TIME_SEQUENCE];
			factor_ptA = 1.0;
			ptB = NULL;
			factor_ptB = 0.0;
			return (ptA != NULL);
		}
		ptA = timeSeqPtrs[relTimeInt % SIZE_PT_TIME_SEQUENCE];
		ptB = timeSeqPtrs[(relTimeInt + 1) % SIZE_PT_TIME_SEQUENCE];
		factor_ptB = relTime - relTimeInt;
		factor_ptA = 1.0 - factor_ptB;
		return ((ptA != NULL) && (ptB != NULL));
	}

	unsigned int ip, im;
	double diffp2Target, diffm2Target, del;
	for (unsigned int pt = 0; pt < SIZE_PT_TIME_SEQUENCE; ++pt)
	{
		if (pt < curPos)
		{
			ip = curPos - pt;
			im = ip - 1;
		}
		else if (pt > curPos)
		{
			ip = curPos - pt + SIZE_PT_TIME_SEQUENCE;
			im = ip - 1;
		}
		else
		{
			ip = 0;
			im = SIZE_PT_TIME_SEQUENCE - 1;
		}
		ptA = timeSeqPtrs[ip];
		ptB = timeSeqPtrs[im];
		if (ptA == NULL)
			continue;
		diffp2Target = ptA->interface_time - timeIn;
		if (fabs(diffp2Target) < g_slf_conf->g_fp_tol_time)
		{
			factor_ptA = 1.0;
			ptB = NULL;
			factor_ptB = 0.0;
			return true;
		}
		if (ptB == NULL)
			continue;
		diffm2Target = timeIn - ptB->interface_time;
		if ((diffp2Target > 0.0) && (diffm2Target > 0.0))
		{
			del = ptA->interface_time - ptB->interface_time;
			factor_ptA = diffm2Target / del;
			factor_ptB = 1.0 - factor_ptA;
			return true;
		}
	}
	return false;
}

bool SL_interfacePPtData_Time_Seq::Interpolate_Pt_Solution_At_time(double timeIn, SL_interfacePPtData*& ptSlnPtr, bool& ptSlnPtr_Deletable)
{
	SL_interfacePPtData *ptA, *ptB;
	double factor_ptA, factor_ptB;
	bool found = Find_Pts_Around_Time(timeIn, ptA, factor_ptA, ptB, factor_ptB);
	if (found == false)
	{
		ptSlnPtr = NULL;
		ptSlnPtr_Deletable = false;
		return false;
	}
	if (ptB == NULL)
	{
		ptSlnPtr = ptA;
		ptSlnPtr_Deletable = false;
		return true;
	}
	ptSlnPtr = new SL_interfacePPtData();
	ptSlnPtr_Deletable = true;
	ptSlnPtr->interface_time = timeIn;
	for (int eside = 0; eside < NUM_SIDES; ++eside)
	{
		LinearCombination2Vecs(ptA->sl_side_ptData[eside].sigma_downstream_final, ptB->sl_side_ptData[eside].sigma_downstream_final, factor_ptA, factor_ptB, ptSlnPtr->sl_side_ptData[eside].sigma_downstream_final);
		LinearCombination2Vecs(ptA->sl_side_ptData[eside].v_downstream_final, ptB->sl_side_ptData[eside].v_downstream_final, factor_ptA, factor_ptB, ptSlnPtr->sl_side_ptData[eside].v_downstream_final);
		LinearCombination2Vecs(ptA->sl_side_ptData[eside].u_downstream_final, ptB->sl_side_ptData[eside].u_downstream_final, factor_ptA, factor_ptB, ptSlnPtr->sl_side_ptData[eside].u_downstream_final);
	}
	ptSlnPtr->interface_damage_final = factor_ptA * ptA->interface_damage_final + factor_ptB * ptB->interface_damage_final;
	ptSlnPtr->interface_damage_source_final = factor_ptA * ptA->interface_damage_source_final + factor_ptB * ptB->interface_damage_source_final;
	ptSlnPtr->maxEffDelU = factor_ptA * ptA->maxEffDelU + factor_ptB * ptB->maxEffDelU;
	ptSlnPtr->tsrStage = ptB->tsrStage;

#if RING_PROBLEM
	ptSlnPtr->v_r_final = factor_ptA * ptA->v_r_final + factor_ptB * ptB->v_r_final;
	ptSlnPtr->v_r_source_final = factor_ptA * ptA->v_r_source_final + factor_ptB * ptB->v_r_source_final;
	ptSlnPtr->sigma_theta_source_final = factor_ptA * ptA->sigma_theta_source_final + factor_ptB * ptB->sigma_theta_source_final;
#endif
	return true;
}

SL_Interface_Temp_PtData_OneSide::SL_Interface_Temp_PtData_OneSide()
{
	SetZero();
}
void SL_Interface_Temp_PtData_OneSide::SetZero()
{
	setValue(v_downstream_latestValue, 0.0);
	setValue(sigma_downstream_latestValue, 0.0);
	setValue(u_downstream_latestValue, 0.0);
}

void SL_Interface_Temp_PtData_OneSide::Read(istream& in)
{
	SetZero();
	for (int i = 0; i < NUMSLRMN; ++i)
	{
		setValue(v_mode_Star[i], 0.0);
		setValue(sigma_mode_Star[i], 0.0);
	}
	setValue(v_Star, 0.0);
	setValue(sigma_Star, 0.0);

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
		if (buf == "v_downstream_latestValue")
		{
			ReadV(v_downstream_latestValue, in);
		}
		else if (buf == "sigma_downstream_latestValue")
		{
			ReadV(sigma_downstream_latestValue, in);
		}
		else if (buf == "u_downstream_latestValue")
		{
			ReadV(u_downstream_latestValue, in);
		}
		else if (buf == "v_mode_Star_I")
		{
			ReadV(v_mode_Star[rmode_stick], in);
		}
		else if (buf == "v_mode_Star_II")
		{
#if DiM2a3_F
			ReadV(v_mode_Star[rmode_slip], in);
#else
			VEC tmp;
			ReadV(tmp, in);
#endif
		}
		else if (buf == "v_mode_Star_III")
		{
			ReadV(v_mode_Star[rmode_sep], in);
		}
		else if (buf == "sigma_mode_Star_I")
		{
			ReadV(sigma_mode_Star[rmode_stick], in);
		}
		else if (buf == "sigma_mode_Star_II")
		{
#if DiM2a3_F
			ReadV(sigma_mode_Star[rmode_slip], in);
#else
			VEC tmp;
			ReadV(tmp, in);
#endif
		}
		else if (buf == "sigma_mode_Star_III")
		{
			ReadV(sigma_mode_Star[rmode_sep], in);
		}
		else if (buf == "v_Star")
		{
			ReadV(v_Star, in);
		}
		else if (buf == "sigma_Star")
		{
			ReadV(sigma_Star, in);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

SL_Interface_Temp_PtData_IntContFrac::SL_Interface_Temp_PtData_IntContFrac()
{
	for (int i = 0; i < s1i_SIZE; ++i)
		interface_scalar_Vals[i] = 0.0;
	interface_scalar_Vals[s1i_absStick] = 1.0;
	interface_scalar_Vals[s1i_contactRel] = 1.0;
	interface_scalar_Vals[s1i_stickRel] = 1.0;
}

void SL_Interface_Temp_PtData_IntContFrac::SetZero()
{
	for (int i = 0; i < s1i_SIZE; ++i)
		interface_scalar_Vals[i] = 0.0;
	sigma_I_nt_parts.SetZero();
	del_u_nt_parts.SetZero();
	del_v_nt_parts.SetZero();
}

void SL_Interface_Temp_PtData_IntContFrac::Read(istream& in)
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
		if (buf == "interface_scalar_Vals")
		{
			for (int i = 0; i < s1i_SIZE; ++i)
				READ_NDOUBLE(in, buf, interface_scalar_Vals[i]);
		}
		else if (buf == "sigma_I_nt_parts")
		{
			sigma_I_nt_parts.Read(in);
		}
		else if (buf == "del_u_nt_parts")
		{
			del_u_nt_parts.Read(in);
		}
		else if (buf == "del_v_nt_parts")
		{
			del_v_nt_parts.Read(in);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void Pt_Error_Data::PrintShort(ostream & out, double time, int itern)
{
	out << time << '\t';
	out << itern << '\t';
	out << norm_del_s << '\t';
	out << norm_del_v << '\t';
	out << del_damage << '\t';
	out << damage_source_tau << '\t';
	out << del_sep2cont_c;
}

void SL_Interface1DPtSeq_Short::AddPt(double time, double vL, double vR, double sigma, double delT2Keep)
{
	SL_Interface1DPt_Short pt;
	pt.sigma = sigma;
	pt.vL = vL;
	pt.vR = vR;
	pt.time = time;
	ptHistory.push_back(pt);
	if (delT2Keep < 0)
		return;
	double delTlimit = time - 1.00001 * delT2Keep;
	unsigned int sz = ptHistory.size() - 1;
	while (sz > 0)
	{
		if (ptHistory[0].time < delTlimit)
		{
			ptHistory.pop_front();
			--sz;
		}
		else
			return;
	}
}

bool SL_Interface1DPtSeq_Short::GetPt(double time, double& vL, double& vR, double& sigma) const
{
	unsigned int indexBase = 0, i;
	bool ti_equal_time = false;
	double time_i, time_ip1, factor_i, factor_ip1;
	int retStat = -1;
	bool found = false;
	unsigned int cntr = 0;
	while (!found & (cntr++ < 10))
	{
		retStat = Get_timeIndex(time, indexBase, i, ti_equal_time, time_i, time_ip1, factor_i, factor_ip1);
		if (retStat == 0)
		{
			found = true;
			break;
		}
		else if (retStat == -2) // out of bound
			return false;
	}
	if (found)
	{
		const SL_Interface1DPt_Short* ptr = &ptHistory[i];
		if (ti_equal_time)
		{
			vL = ptr->vL, vR = ptr->vR, sigma = ptr->sigma;
			return true;
		}
		const SL_Interface1DPt_Short* ptrN = &ptHistory[i + 1];
		vL = factor_i * ptr->vL + factor_ip1 * ptrN->vL;
		vR = factor_i * ptr->vR + factor_ip1 * ptrN->vR;
		sigma = factor_i * ptr->sigma + factor_ip1 * ptrN->sigma;
		return true;
	}
	return GetPtSlow(time, vL, vR, sigma);

}

bool SL_Interface1DPtSeq_Short::GetPtSlow(double time, double& vL, double& vR, double& sigma) const
{
	unsigned int sz = ptHistory.size();
	unsigned int i;
	bool found = false;
	for (i = 0; i < sz; ++i)
	{
		if (ptHistory[i].time > time)
		{
			found = true;
			break;
		}
	}
	if (!found)
		i = sz - 1;
	bool sameVal = DoublesAreEqual(ptHistory[i].time, time);
	if ((i == 0) || (!found))
	{
		if (!sameVal)
			return false;
	}
	const SL_Interface1DPt_Short* ptNext = &ptHistory[i];
	if (sameVal)
	{
		vL = ptNext->vL, vR = ptNext->vR, sigma = ptNext->sigma;
		return true;
	}
	const SL_Interface1DPt_Short* ptPrev = &ptHistory[i - 1];
	double factorPrev = (ptNext->time - time) / (ptNext->time - ptPrev->time);
	double factorNext = 1.0 - factorPrev;
	vL = ptNext->vL * factorNext + ptPrev->vL * factorPrev;
	vR = ptNext->vR * factorNext + ptPrev->vR * factorPrev;
	sigma = ptNext->sigma  * factorNext + ptPrev->sigma * factorPrev;
	return true;
}

int SL_Interface1DPtSeq_Short::GetPt(unsigned int index, double& time, double& vL, double& vR, double& sigma)
{
	if (index < 0)
		index = ptHistory.size() + index;

	SL_Interface1DPt_Short* ptr = &ptHistory[index];
	time = ptr->time;
	vL = ptr->vL;
	vR = ptr->vR;
	sigma = ptr->sigma;
	return index;
}

int SL_Interface1DPtSeq_Short::Get_timeIndex(double time, unsigned int& indexBase, unsigned int& i, bool& ti_equal_time, double& time_i, double& time_ip1, double& factor_i, double& factor_ip1) const
{
	double timeBase = ptHistory[indexBase].time;
	double delT;
	unsigned int lastInd = ptHistory.size() - 1;
	if (indexBase > 0)
		delT = timeBase - ptHistory[indexBase - 1].time;
	else if (indexBase < lastInd)
		delT = ptHistory[indexBase + 1].time - timeBase;
	else
		return -2; // out of bound

	static double relTol = 1e-5;
	double tolAbs = relTol * delT;
	double t0 = ptHistory[0].time - tolAbs;
	double tLast = ptHistory[lastInd].time + tolAbs;
	if ((time < t0) || (time > tLast)) // out of bound
		return -2;
	if (time > timeBase + tolAbs) // (time > timeBase)
		i = MAX(0, (int)floor(indexBase + (time - timeBase) / delT + relTol));
	else if (time < timeBase - tolAbs)
		i = MAX(0, (int)floor(indexBase - (timeBase - time) / delT + relTol));
	else
	{
		i = indexBase;
		ti_equal_time = true;
		time_i = timeBase;
		time_ip1 = timeBase;
		factor_i = 1.0;
		factor_ip1 = 0.0;
		return 0;
	}
	time_i = ptHistory[i].time;
	double timeChange = time - time_i;
	if (fabs(timeChange) < tolAbs)
	{
		ti_equal_time = true;
		time_ip1 = time_i;
		factor_i = 1.0;
		factor_ip1 = 0.0;
		return 0;
	}
	if (timeChange < 0.0) // time must be >= time_i for acceptable solution, so this signs an unacceptable solution
	{
		indexBase = i;
		return -1;
	}
	int ip1 = i + 1;
	if (ip1 > lastInd)
		return -2;
	ti_equal_time = false;
	time_ip1 = ptHistory[ip1].time;
	if (time - time_ip1 > tolAbs)
	{
		indexBase = ip1;
		return 1;
	}
	factor_i = MIN(MAX((time_ip1 - time) / (time_ip1 - time_i), 0.0), 1.0);
	factor_ip1 = 1.0 - factor_i;
	return 0;
}

SL_Interface1DPtSeq_Short* g_seq_short = NULL;

Pt_Error_Data::Pt_Error_Data()
{
	SetZero();
}

void Pt_Error_Data::SetZero()
{
	norm_del_v = 0.0;
	norm_del_s = 0.0;
	del_damage = 0.0;
	damage_source_tau = 0.0;
	del_sep2cont_c = 0.0;
}

void Pt_Error_Data::Read(istream& in)
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
		if (buf == "norm_del_v")
		{
			READ_NDOUBLE(in, buf, norm_del_v);
		}
		else if (buf == "norm_del_s")
		{
			READ_NDOUBLE(in, buf, norm_del_s);
		}
		else if (buf == "del_damage")
		{
			READ_NDOUBLE(in, buf, del_damage);
		}
		else if (buf == "damage_source_tau")
		{
			READ_NDOUBLE(in, buf, damage_source_tau);
		}
		else if (buf == "del_sep2cont_c")
		{
			READ_NDOUBLE(in, buf, del_sep2cont_c);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void Output_ScalarValuesAux(ostream& out, IOF_type iot, double space_or_time, double* interface_scalar_Vals, SLFF_InterfacialDamageModeType damageOffOnMix, SLFF_ContactType	contactOffOnMix, SLFF_SlipType slipOffOnMix)
{
	int sz = 4;
	if (interface_scalar_Vals != NULL)
		sz += s1i_SIZE;
	if (iot == iof_ascii)
	{
		out << sz << '\t';
		out << space_or_time << '\t' << (int)damageOffOnMix << '\t' << (int)contactOffOnMix << '\t' << (int)slipOffOnMix;
		if (interface_scalar_Vals != NULL)
			for (int i = 0; i < s1i_SIZE; ++i)
				out << '\t' << interface_scalar_Vals[i];
		out << '\n';
	}
	else if (iot == iof_binary)
	{
		out.write((char*)&sz, sizeof(sz));
		out.write((char*)&space_or_time, sizeof(space_or_time));
		out.write((char*)&damageOffOnMix, sizeof(damageOffOnMix));
		out.write((char*)&contactOffOnMix, sizeof(contactOffOnMix));
		out.write((char*)&slipOffOnMix, sizeof(slipOffOnMix));

		if (interface_scalar_Vals != NULL)
			for (int i = 0; i < s1i_SIZE; ++i)
				out.write((char*)&interface_scalar_Vals[i], sizeof(double));
	}
}
