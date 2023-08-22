#include "SLDescriptorData.h"
#include "LAfuncsFinalStep.h"
#include "SLBulk_Properties.h"
#include "SimpleConfigMaker.h"

SLDescriptorData g_SL_desc_data;


SLDescriptorData::SLDescriptorData()
{
	setValue(a_xt_prob, 0.0);
	ring_p0 = 1.0;
	load_number = AXT_LN;

	tdLoadComputer = NULL;
	tdLoadType = lmt_Incident;
	tdLoadSide = SDL;
	b_tdLoad_stressScale = false;
	b_tdLoad_velScale = false;
	sz_load_parameters = 0;

	tdLoad_Zi = 0.0, tdLoad_Zo = 0.0, tdLoad_ZEffective = 0.0;
	bndryLoad_inputEnergy = 0.0, tdLoad_stressScale = 0.0, tdLoad_velScale = 0.0, tdLoad_loadScale = 0.0, tdAmbientProjectileLength = 0.0;

	DiracLoadingPtr = NULL;
}

SLDescriptorData::~SLDescriptorData()
{
	if (tdLoadComputer != NULL)
		delete tdLoadComputer;
	if (DiracLoadingPtr != NULL)
		delete DiracLoadingPtr;
}

void SLDescriptorData::Read(istream & in)
{
	string key;	map<string, string>* mpPtr;
	double value = -1;
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
		else
		{
			while ((buf != "infile_icbc") && (!in.eof()))
				READ_NSTRING(in, buf, buf);
			if (in.eof())
				THROW("Reached end of file looking for infile_icbc\n");
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
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "load_number")
		{
			READ_NINTEGER(in, buf, load_number);
		}
		else if (buf == "DiracLoading")
		{
			DiracLoadingPtr = new DiracLoading();
			DiracLoadingPtr->Read_DiracLoading(in);
			if (!DiracLoadingPtr->IsActive())
			{
				delete DiracLoadingPtr;
				DiracLoadingPtr = NULL;
			}
		}
		else if (buf == "a_xt_prob")
		{
			ReadV(a_xt_prob, in);
		}
		else if (buf == "ring_p0")
		{
			READ_NDOUBLE(in, buf, ring_p0);
		}
		else if (buf == "load_parameters")
		{
			int sz;
			READ_NINTEGER(in, buf, sz);
			load_parameters.resize(sz);
			for (int i = 0; i < sz; ++i)
				READ_NDOUBLE(in, buf, load_parameters[i]);
		}
		else if (buf == "LoadComputer")
		{
			tdLoadComputer = new TD_FD();
			in >> *tdLoadComputer;
		}
		else if (buf == "LoadType")
		{
			in >> tdLoadType;
		}
		else if (buf == "LoadSide")
		{
			string side;
			READ_NSTRING(in, buf, side);
			to_lower(side);
			if ((side == "left") || (side == "l") || (side == "0"))
				tdLoadSide = SDL;
			else if ((side == "right") || (side == "r") || (side == "1"))
				tdLoadSide = SDR;
			else
			{
				cout << "side\n" << side << '\n';
				THROW("Invalid side symbol\n");
			}
		}
		else if (buf == "stressScale")
		{
			READ_NDOUBLE(in, buf, tdLoad_stressScale);
			key = "stressScale";
			if (Find_Version_Value(key, value, mpPtr) == true)
				tdLoad_stressScale = value;
			b_tdLoad_stressScale = true;
		}
		else if (buf == "velScale")
		{
			READ_NDOUBLE(in, buf, tdLoad_velScale);
			key = "velScale";
			if (Find_Version_Value(key, value, mpPtr) == true)
				tdLoad_velScale = value;
			b_tdLoad_velScale = true;
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	if (load_number == AXT_LN)
	{
		if (load_parameters.size() == 0)
		{
			load_parameters.resize(1);
			load_parameters[0] = a_xt_prob[0];
		}
		double value = 0.0; map<string, string>* mpPtr;
		if (Find_Version_Value("la", value, mpPtr))
		{
			load_parameters[0] = pow(10.0, value);
		}
		g_logout << "\tla\t" << value << "\tlrate\t" << load_parameters[0];
	}
	Finalize_SLDescriptorData();
}

void SLDescriptorData::Read(string configNameIn)
{
	fstream in(configNameIn.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "configNameIn\t" << configNameIn << '\n';
		THROW("Cannot open file\n");
	}
	Read(in);
}

void SLDescriptorData::Finalize_SLDescriptorData()
{
	int sz = load_parameters.size();
	if (load_number == AXT_LN)
	{
		if ((sz > 0) && (fabs(a_xt_prob[0]) < 1e-15))
		{
			unsigned int szz = MIN(sz, DiM);
			for (unsigned int i = 0; i < szz; ++i)
				a_xt_prob[i] = load_parameters[i];
		}
	}
	else if (load_number == LC_PIECE_WISE_LIN)
	{
		bndryLoad_inputEnergy = 0.0;
		unsigned int szTimes = load_parameters.size() / 2;
		if (szTimes <= 1)
			bndryLoad_inputEnergy = 1.0;
		else
		{
			for (int j = 0; j < (int)szTimes - 1; ++j)
			{
				int st = 2 * j;
				double tp = load_parameters[st];
				double lp = load_parameters[st + 1];
				double tn = load_parameters[st + 2];
				double ln = load_parameters[st + 3];
				bndryLoad_inputEnergy += 0.5 * (tn - tp) * (lp + ln);
			}
		}
		bndryLoad_inputEnergy = fabs(bndryLoad_inputEnergy); // this should really be scaled by the correct impedance, if it's Dirichlet or Neumann BC, but for the <E> = 1, <rho> = 1 problems considered this is fine
	}
	else if (load_number == LC_CYCLIC_TEN_COMP)
	{
		double period = 1.0;
		double relTensile = 0.25;
		double sigmaTensile = 0.1;
		double relRamp = 0.05;
		double sigmaCompressive = 0.0;
		bool sigmaCompRead = false;
		double additiveFactor = 0.0;

		if (sz > 0)
		{
			period = load_parameters[0];
			if (sz > 1)
			{
				relTensile = load_parameters[1];
				if (sz > 2)
				{
					sigmaTensile = load_parameters[2];
					if (sz > 3)
					{
						relRamp = load_parameters[3];
						if (sz > 4)
						{
							sigmaCompressive = load_parameters[4];
							sigmaCompRead = true;

							if (sz > 5)
								additiveFactor = load_parameters[5];
						}
					}
				}
			}
		}
		if (sigmaCompRead == false)
			sigmaCompressive = -sigmaTensile;
		load_parameters.resize(6);

		load_parameters[0] = period;
		load_parameters[1] = relTensile;
		load_parameters[2] = sigmaTensile;
		load_parameters[3] = relRamp;
		load_parameters[4] = sigmaCompressive;
		load_parameters[5] = additiveFactor;
	}
	sz_load_parameters = load_parameters.size();
}

bool SLDescriptorData::GetNonRingNonZeroTerm_SourceTerm(double x, double t, GID bulkFlag, VEC & source_v, VEC & source_sigma) const
{
//	if ((load_number == AXT_LN) || (tdLoadComputer != NULL))
		return false;
//	THROW("Non-axt problem is not implemented.\n");
}

double SLDescriptorData::GetRing_p(double x, double t)
{
	// later add option to have spacetime variable p
	return ring_p0;
}

void SLDescriptorData::Get_IC(SL_Bulk_Properties* bulk_Ptr, GID bulkFlag, double x, VEC& u_0, VEC& v_0, VEC& eps_0
#if RING_PROBLEM
	, double& v_r
#endif
) const
{
#if RING_PROBLEM
	v_r = 0.0;
#endif

	setValue(u_0, 0.0);
	setValue(eps_0, 0.0);
	if (load_number == AXT_LN)
	{
		for (int i = 0; i < DiM; ++i)
		{
			v_0[i] = a_xt_prob[i] * x;
			u_0[i] = v_0[i] * gt0;
			eps_0[i] = a_xt_prob[i] * gt0;
		}
	}
	else
	{
		for (int i = 0; i < DiM; ++i)
			v_0[i] = 0.0;
	}
	if (tdLoadComputer != NULL)
	{
		if (tdLoadType == lmt_Impact)
		{
			if (tdLoadSide == SDL)
			{
				if (bulkFlag == g_bulkFlags_IncImp[SDL])
					v_0[0] = -tdLoad_velScale;
			}
			else if (tdLoadSide == SDR)
			{
				if (bulkFlag == g_bulkFlags_IncImp[SDR])
					v_0[0] = tdLoad_velScale;
			}
		}
	}
}

void SLDescriptorData::GetBoundaryConditionValue_LeftRight(bool isLeft, SL_Elastic_InterfaceProperties* ts_bulkProps, double x, double t, VEC& BC_val, int timeIndex) const
{
	SL_Bulk_Properties *bulk_Ptr;
	if (isLeft)
		bulk_Ptr = ts_bulkProps->bulk_rightPtr;
	else
		bulk_Ptr = ts_bulkProps->bulk_leftPtr;
	if (load_number == AXT_LN)
	{
		if (ts_bulkProps->isIsoInterface)
		{
			for (int i = 0; i < DiM; ++i)
			{
				if (ts_bulkProps->directionalBCType[i] == bct_Dirichlet)
					BC_val[i] = a_xt_prob[i] * x;
				else if (ts_bulkProps->directionalBCType[i] == bct_Neumann)
					BC_val[i] = a_xt_prob[i] * (t + gt0) * bulk_Ptr->effective_CsIso[i];
				else if (ts_bulkProps->directionalBCType[i] == bct_Characteristics)
				{
					if (isLeft)
						BC_val[i] = a_xt_prob[i] * (t + gt0) * bulk_Ptr->effective_CsIso[i] - a_xt_prob[i] * x;
					else
						BC_val[i] = a_xt_prob[i] * (t + gt0) * bulk_Ptr->effective_CsIso[i] + a_xt_prob[i] * x;
				}
			}
		}
#if COMPUTE_ANISO_BULK
		else
		{
			BoundaryConditionT bct = ts_bulkProps->directionalBCType[0];
			if (bct == bct_Dirichlet)
				for (int i = 0; i < DiM; ++i)
					BC_val[i] = a_xt_prob[i] * x;
			else
			{
				VEC eps, sig;
				for (int i = 0; i < DiM; ++i)
					eps[i] = a_xt_prob[i] * t;
				ProductMatVec(bulk_Ptr->C, eps, sig);
				if (bct == bct_Neumann)
					CopyVec(sig, BC_val);
				else
				{
					VEC v, Fv, Esigma;
					for (int i = 0; i < DiM; ++i)
						v[i] = a_xt_prob[i] * x;

					ProductMatVec(bulk_Ptr->E, sig, Esigma);
					ProductMatVec(bulk_Ptr->F, v, Fv);
					if (isLeft)
						SubtractVec(Esigma, Fv, BC_val);
					else
						AddVec(Esigma, Fv, BC_val);
				}
			}
		}
#endif
		return;
	}
	setValue(BC_val, 0.0);
	if (tdLoadComputer != NULL)
	{
		bool isDirichlet = (tdLoadType == lmt_Dirichlet);
		bool isNeumann = (tdLoadType == lmt_Neumann);
		if (isDirichlet || isNeumann)
		{
			double value = tdLoad_loadScale * tdLoadComputer->ComputeSpecific_TD_Value(t);
			if (isDirichlet)
			{
				if (isLeft)
					value = -value;
			}
			BC_val[0] = value;
		}
		// nothing to be done for else as the BC is traction free for impact and 0 transmitting for incident case
		return;
	}
	if (DiracLoadingPtr != NULL)
	{
		BC_val[0] = 0.0;
		if (isLeft)
			BC_val[0] = DiracLoadingPtr->GetBaseValue(t, timeIndex);
		return;
	}
	if (load_number == SQUARE_PULSE)
	{
		if (isLeft)
		{
			double rampTime = 25, load = 0.1;
			unsigned int sz = load_parameters.size();
			if (sz > 0)
			{
				rampTime = load_parameters[0];
				if (sz > 1)
					load = load_parameters[1];
			}
			if (t < rampTime)
				BC_val[0] = load;
		}
		return;
	}
	if (load_number == LC_PIECE_WISE_LIN)
	{
		if (!isLeft)
		{
			BC_val[0] = 0.0;
			return;
		}
		double load = 1;
		if (sz_load_parameters >= 2)
		{
			double t0 = load_parameters[0], l0 = load_parameters[1];
			if (t < t0)
				load = 0.0;
			else if (sz_load_parameters == 2)
				load = load_parameters[1];
			else
			{
				unsigned int szd2 = sz_load_parameters / 2;
				unsigned int j;
				double tv;
				for (j = 0; j < szd2; ++j)
				{
					tv = load_parameters[2 * j];
					if (tv > t)
						break;
				}
				if (j == szd2)
					load = 0.0;
				else
				{
					double tp = load_parameters[2 * j - 2];
					double lp = load_parameters[2 * j - 1];
					double tn = load_parameters[2 * j];
					double ln = load_parameters[2 * j + 1];
					double Np = (tn - t) / (tn - tp), Nn = 1.0 - Np;
					load = lp * Np + ln * Nn;
				}
			}
		}
		BC_val[0] = load;
		return;
	}
	THROW("The option for other load numbers is not implemented\n");
}

void SLDescriptorData::Get1InterfaceImpingingCharacteristics(SL_Elastic_InterfaceProperties* ts_bulkProps, double t, VEC& wlSide_rGoing, VEC& wrSide_lGoing)
{
	if (load_number == LC_CYCLIC_TEN_COMP)
	{
		double period = load_parameters[0];
		double relTensile = load_parameters[1];
		double sigmaTensile = load_parameters[2];
		double relRamp = load_parameters[3];
		double sigmaCompressive = load_parameters[4];
		double additiveFactor = load_parameters[5];

		double relPer = t / period;
		double flr = floor(relPer);
		if (fabs(additiveFactor) > 0.0)
		{
			double mult = (flr + 1) * additiveFactor;
			sigmaTensile *= mult;
			sigmaCompressive *= mult;
		}
		relPer -= flr;
		if (relPer < 0.0)
			relPer = 0.0;
		else if (relPer > 1.0)
			relPer = 1.0;
		double w;
		if (relPer <= relRamp)
			w = relPer / relRamp * sigmaTensile;
		else if (relPer <= relTensile - relRamp)
			w = sigmaTensile;
		else if ((relPer > (relTensile - relRamp)) && (relPer <= (relTensile + relRamp)))
		{
			double rem = 0.5 * (relPer - (relTensile - relRamp)) / relRamp;
			w = sigmaCompressive * rem + sigmaTensile * (1.0 - rem);
		}
		else if (relPer <= (1.0 - relRamp))
			w = sigmaCompressive;
		else
		{
			double rem = (1.0 - relPer) / relRamp;
			w = sigmaCompressive * rem;
		}
		for (int i = 0; i < DiM; ++i)
		{
			wlSide_rGoing[i] = w;
			wrSide_lGoing[i] = w;
		}
		return;
	}

	if (load_number == AXT_LN)
	{
		if (ts_bulkProps->isIsoInterface)
		{
			for (int i = 0; i < DiM; ++i)
			{
				wlSide_rGoing[i] = a_xt_prob[i] * (t + gt0) * ts_bulkProps->bulk_leftPtr->effective_CsIso[i];
				wrSide_lGoing[i] = a_xt_prob[i] * (t + gt0) * ts_bulkProps->bulk_rightPtr->effective_CsIso[i];
			}
		}
#if COMPUTE_ANISO_BULK
		else
		{
			VEC eps, sigL, sigR;
			for (int i = 0; i < DiM; ++i)
				eps[i] = a_xt_prob[i] * t;
			ProductMatVec(ts_bulkProps->bulk_leftPtr->C, eps, sigL);
			ProductMatVec(ts_bulkProps->bulk_rightPtr->C, eps, sigR);
			ProductMatVec(ts_bulkProps->bulk_leftPtr->E, sigL, wlSide_rGoing);
			ProductMatVec(ts_bulkProps->bulk_rightPtr->E, sigR, wrSide_lGoing);
		}
#endif
		return;
	}
	THROW("The option for other load numbers is not implemented\n");
}

double SLDescriptorData::GetLoadingTimeScale() const
{
	if (load_number == AXT_LN)
		return 1.0 / load_parameters[0];
	if (tdLoadComputer != NULL)
		return tdLoadComputer->timeScale;
	if (load_number == LC_PIECE_WISE_LIN)
	{
		if (sz_load_parameters > 2)
			return load_parameters[2];
		return 1.0;
	}
	THROW("implement the option if needed\n");
	return 1.0;
}

void SLDescriptorData::Finalize_tdLoadParameters(double finalTime, double delT, double amientProjLength, double Ein, double rhoin, double Eout, double rhoout)
{
	delT *= 0.1;
	if (tdLoadComputer == NULL)
		return;
	tdLoad_Zi = sqrt(Ein * rhoin);
	if (Eout < 0)
		Eout = Ein;
	if (rhoout < 0)
		rhoout = rhoin;
	tdLoad_Zo = sqrt(Eout * rhoout);
	if (tdLoad_Zo < 0)
		tdLoad_Zo = tdLoad_Zi;
	tdAmbientProjectileLength = amientProjLength;
	bool isNeumannOrDirichlet = ((tdLoadType == lmt_Neumann) || (tdLoadType == lmt_Dirichlet));
	if (!isNeumannOrDirichlet)
	{
		if (tdLoadType == lmt_Incident)
		{
			if (b_tdLoad_stressScale)
			{
				tdLoad_velScale = tdLoad_stressScale / tdLoad_Zo;
				b_tdLoad_velScale = true;
			}
			else if (b_tdLoad_velScale)
			{
				tdLoad_stressScale  = tdLoad_velScale * tdLoad_Zo;
				b_tdLoad_stressScale = true;
			}
			else
			{
				THROW("stress or vel are not set\n");
			}
			tdLoad_loadScale = tdLoad_stressScale;
			tdLoad_ZEffective = (tdLoad_Zi + tdLoad_Zo);
			tdLoad_ZEffective = 0.25 * tdLoad_ZEffective * tdLoad_ZEffective / tdLoad_Zi;
			double integral_loadOut, integral_load2Out;
			tdLoadComputer->Compute_Inegral_load_load2(finalTime, delT, integral_loadOut, integral_load2Out);
			bndryLoad_inputEnergy = tdLoad_stressScale * tdLoad_stressScale * integral_load2Out / tdLoad_ZEffective;
		}
		else if (tdLoadType == lmt_Impact)
		{
			tdLoad_ZEffective = 2.0 * tdLoad_Zi * tdLoad_Zo / (tdLoad_Zi + tdLoad_Zo);

			if (b_tdLoad_stressScale)
			{
				tdLoad_velScale = tdLoad_stressScale / tdLoad_ZEffective;
				b_tdLoad_velScale = true;
			}
			else if (b_tdLoad_velScale)
			{
				tdLoad_stressScale = tdLoad_velScale * tdLoad_ZEffective;
				b_tdLoad_stressScale = true;
			}
			else
			{
				THROW("stress or vel are not set\n");
			}
			tdLoad_loadScale = tdLoad_velScale;
			bndryLoad_inputEnergy = 0.5 * rhoout * tdLoad_velScale * tdLoad_velScale * tdAmbientProjectileLength;
		}
		else
		{
			cout << "tdLoadType\t" << tdLoadType << '\n';
			THROW("Invalid tdLoadType\t");
		}
	}
	else
	{
		tdLoad_ZEffective = tdLoad_Zi;
		if (b_tdLoad_stressScale)
		{
			tdLoad_velScale = tdLoad_stressScale / tdLoad_ZEffective;
			b_tdLoad_velScale = true;
		}
		else if (b_tdLoad_velScale)
		{
			tdLoad_stressScale = tdLoad_velScale * tdLoad_ZEffective;
			b_tdLoad_stressScale = true;
		}
		else
		{
			THROW("stress or vel are not set\n");
		}
		if (tdLoadType == lmt_Neumann)
			tdLoad_loadScale = tdLoad_stressScale;
		else
			tdLoad_loadScale = tdLoad_velScale;
		double integral_loadOut, integral_load2Out;
		tdLoadComputer->Compute_Inegral_load_load2(finalTime, delT, integral_loadOut, integral_load2Out);
		bndryLoad_inputEnergy = tdLoad_stressScale * tdLoad_stressScale * integral_load2Out / tdLoad_ZEffective;
		tdAmbientProjectileLength = 0.0;
	}
	Print_tdLoad();
}

void SLDescriptorData::Print_tdLoad()
{
	string fileName = g_prefileName + "/" + "__tdLoad.txt";
	fstream out(fileName.c_str(), ios::out);

	out << "bndryLoad_inputEnergy\t" << bndryLoad_inputEnergy << '\n';
	out << "tdLoadType\t" << tdLoadType << '\n';
	out << "tdLoadSide\t" << tdLoadSide << '\n';
	out << "tdLoad_stressScale\t" << tdLoad_stressScale << '\n';
	out << "tdLoad_velScale\t" << tdLoad_velScale << '\n';
	out << "tdLoad_loadScale\t" << tdLoad_loadScale << '\n';
	out << "tdLoad_Zi\t" << tdLoad_Zi << '\n';
	out << "tdLoad_Zo\t" << tdLoad_Zo << '\n';
	out << "tdLoad_ZEffective\t" << tdLoad_ZEffective << '\n';
	out << "tdAmbientProjectileLength\t" << tdAmbientProjectileLength << '\n';
}

void SLDescriptorData::Read_tdLoad()
{
	string fileName = g_prefileName + "/" + "__tdLoad.txt";
	fstream in(fileName.c_str(), ios::in);
	if (!in.is_open())
		return;
	string buf;
	in >> buf >> bndryLoad_inputEnergy;
	in >> buf >> tdLoadType;
	in >> buf >> tdLoadSide;
	in >> buf >> tdLoad_stressScale;
	in >> buf >> tdLoad_velScale;
	in >> buf >> tdLoad_loadScale;
	in >> buf >> tdLoad_Zi;
	in >> buf >> tdLoad_Zo;
	in >> buf >> tdLoad_ZEffective;
	in >> buf >> tdAmbientProjectileLength;
}

void CyclicLoading::Initialize_FromParamaters(const vector<double>& paras)
{
	sigma0 = paras[0];
	T = paras[1];
	beta = -1.0;
	tr = -1.0;
	if (paras.size() > 2)
	{
		beta = paras[2];
		if (paras.size() > 3)
			tr = paras[3];
	}
	if (beta < 0)
		beta = 0.8104702832; // ~ .2579807036 PI // sinBeta = 0.7246113538
	double sinBeta = sin(beta);
	if (tr < 0)
		tr = sinBeta;
}

double CyclicLoading::getValue(double time)
{
	double tp = time / T * 2.0 * PI;
	double sinBeta = sin(beta);
	double w0p;
	if (tp < tr)
		w0p = -tp / tr * sinBeta;
	else
	{
		tp -= tr;
		w0p = -(sinBeta + sin(tp));
	}
	return w0p * sigma0;
}

DiracLoading::DiracLoading()
{
	maxTime = 1.0;
	numTimeSteps2Cover = 2;
	fld4Int1 = di_energy;
	Z = 1.0;
	canOverwiseZ = true;
	directionalBCType = bct_Neumann;
	timeStepping_delt = 9.765625000000000e-04; // 1/1024, 1024 mesh with CFL = 1
	totalFactor4Integral = 1.0;
	b_maxTimeRead = false;
	b_numTimeSteps2CoverRead = false;
	b_simple_parabola42TimeSteps = false;
	isActive = true;
}

void DiracLoading::Read_DiracLoading(istream& in)
{
	string key;	map<string, string>* mpPtr;
	double value = -1;
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "isActive")
		{
			READ_NBOOL(in, buf, isActive);
		}
		else if (buf == "maxTime")
		{
			double maxTimeIn;
			READ_NDOUBLE(in, buf, maxTimeIn);
			key = "Dirac_maxTime";
			if (Find_Version_Value(key, value, mpPtr) == true)
				maxTimeIn = value;
			SetDiractMaxTime(maxTimeIn);
		}
		else if (buf == "numTimeSteps2Cover")
		{
			READ_NINTEGER(in, buf, numTimeSteps2Cover);
			key = "Dirac_numTimeSteps2Cover";
			if (Find_Version_Value(key, value, mpPtr) == true)
				numTimeSteps2Cover = (int)round(value);
			b_numTimeSteps2CoverRead = true;
		}
		else if (buf == "fld4Int1")
		{
			in >> fld4Int1;
		}
		else if (buf == "directionalBCType")
		{
			in >> directionalBCType;
		}
		else if (buf == "Z")
		{
			READ_NDOUBLE(in, buf, Z);
		}
		else if (buf == "ZtoUse")
		{
			READ_NDOUBLE(in, buf, Z);
			canOverwiseZ = false;
		}
		else if (buf == "timeStepping_delt")
		{
			READ_NDOUBLE(in, buf, timeStepping_delt);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void DiracLoading::Write_DiracLoading(ostream& out)
{
	out << "{\n";
	out << "fld4Int1\t" << fld4Int1 << '\n';
	out << "directionalBCType\t" << directionalBCType << '\n';
	out << "timeStepping_delt\t" << timeStepping_delt << '\n';
	out << "maxTime\t" << maxTime << '\n';
	out << "numTimeSteps2Cover\t" << numTimeSteps2Cover << '\n';
	out << "Z\t" << Z << '\n';
	out << "}\n";
}

void DiracLoading::SetDiractMaxTime(double maxTimeIn)
{
	if (maxTimeIn < 0.0)
	{
		numTimeSteps2Cover = -(int)maxTimeIn;
		b_numTimeSteps2CoverRead = true;
	}
	else
	{
		maxTime = maxTimeIn;
		b_maxTimeRead = true;
	}
}

void DiracLoading::InitializeFromOutside(double delt, BoundaryConditionT bcTypeIn, double ZIn)
{
	timeStepping_delt = delt;
	directionalBCType = bcTypeIn;
	if (canOverwiseZ)
		Z = ZIn;
	Initialize_DiracLoading();
}

string getName(DiracIntegralT dat)
{
	if (dat == di_energy)
		return "energy";
	if (dat == di_stress)
		return "stress";
	if (dat == di_velocity)
		return "velocity";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, DiracIntegralT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= DiracIntegralT_SIZE)
			THROW("too large of a number\n");
		typeVal = (DiracIntegralT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < DiracIntegralT_SIZE; ++i)
	{
		typeVal = (DiracIntegralT)i; // casting integer to DiracIntegralT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading DiracIntegralT\n");
}

//operator for output
ostream& operator<<(ostream& out, DiracIntegralT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, DiracIntegralT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

ostream& operator<<(ostream& out, const twsvp& dat)
{
	out << dat.t << '\t';
	out << dat.w << '\t';
	out << dat.s << '\t';
	out << dat.v << '\t';
	out << dat.p;
	return out;
}

twsvp::twsvp()
{
	t = 0.0, w = 0.0, s = 0.0, v = 0.0, p = 0.0;
}

void DiracLoading::Initialize_DiracLoading()
{
	if (fld4Int1 == di_energy)
	{
		if (directionalBCType == bct_Dirichlet)
			physicalFactor4Int1 = 1.0 / sqrt(Z);
		else if (directionalBCType == bct_Neumann)
			physicalFactor4Int1 = sqrt(Z);
		else if (directionalBCType == bct_Characteristics)
			physicalFactor4Int1 = 2.0 * sqrt(Z);
	}
	else if (fld4Int1 == di_stress)
	{
		if (directionalBCType == bct_Dirichlet)
			physicalFactor4Int1 = 1.0 / Z;
		else if (directionalBCType == bct_Neumann)
			physicalFactor4Int1 = 1.0;
		else if (directionalBCType == bct_Characteristics)
			physicalFactor4Int1 = 2.0;
	}
	else if (fld4Int1 == di_velocity)
	{
		if (directionalBCType == bct_Dirichlet)
			physicalFactor4Int1 = 1.0;
		else if (directionalBCType == bct_Neumann)
			physicalFactor4Int1 = Z;
		else if (directionalBCType == bct_Characteristics)
			physicalFactor4Int1 = 2.0 * Z;
	}
	if (b_maxTimeRead)
		numTimeSteps2Cover = (int)floor(0.5 * maxTime / timeStepping_delt) * 2;
	if (numTimeSteps2Cover < 2)
		numTimeSteps2Cover = 2;
	maxTime = timeStepping_delt * numTimeSteps2Cover;
	b_simple_parabola42TimeSteps = (numTimeSteps2Cover == 2);
	if (b_simple_parabola42TimeSteps)
	{
		totalFactor4Integral = 0.75 / timeStepping_delt;
		if (fld4Int1 == di_energy)
			totalFactor4Integral = sqrt(totalFactor4Integral);
		totalFactor4Integral *= physicalFactor4Int1;
		cout << "this is too sharp of a time step. Increase Dirac delt by at least a factor of 2 (now it's 2x delt -> should be at least 4x delt\n";
		cout << "maxTime at least should be\t" << 4.0 * timeStepping_delt;
		cout << "or maxTime entered = " << -4.0 << '\n';
#if DB_STRICT_EXIT
		THROW("Code exits because of the error above\n");
#endif
	}
	else
	{
		double dydt;
		double delt3rd = timeStepping_delt / 3.0;
		double integralVal = 0.0;
		offset = 0.5 * maxTime;
		for (unsigned int i = 0; i < (unsigned int)numTimeSteps2Cover / 2; ++i)
		{
			double t0 = i * 2 * timeStepping_delt;
			double t1 = t0 + timeStepping_delt;
			double t2 = t1 + timeStepping_delt;
			double v0 = getInfinitlySmoothDeltaDirac(t0 - offset, dydt, offset, false);
			double v1 = getInfinitlySmoothDeltaDirac(t1 - offset, dydt, offset, false);
			double v2 = getInfinitlySmoothDeltaDirac(t2 - offset, dydt, offset, false);
			// end points end up being computed twice, but it's OK
			integralVal += delt3rd * (v0 + 4.0 * v1 + v2);
		}
		if (fld4Int1 == di_energy)
			totalFactor4Integral = physicalFactor4Int1 / sqrt(integralVal);
		else
			totalFactor4Integral = physicalFactor4Int1 / integralVal;
	}
	inputEnergy_Dirac = 1.0;
	if (fld4Int1 == di_stress)
	{
		if (b_simple_parabola42TimeSteps)
			inputEnergy_Dirac = 0.75 / timeStepping_delt / Z;
		else
			inputEnergy_Dirac = 1.0 / timeStepping_delt / Z; // rought estimate
	}
	else if (fld4Int1 == di_velocity)
	{
		if (b_simple_parabola42TimeSteps)
			inputEnergy_Dirac = 0.75 * Z / timeStepping_delt;
		else
			inputEnergy_Dirac = Z / timeStepping_delt; // rought estimate
	}
}

double DiracLoading::GetBaseValue(double time, int timeIndex)
{
	if (b_simple_parabola42TimeSteps)
	{
		if (timeIndex == -1)
			timeIndex = (int)round(time / timeStepping_delt);
		if (timeIndex == 1)
			return totalFactor4Integral;
		return 0.0;
	}
	static double dydt;
	double fn = getInfinitlySmoothDeltaDirac(time - offset, dydt, offset, false);
	if (fld4Int1 == di_energy)
		fn = sqrt(fn);
	return totalFactor4Integral * fn;
}

void DiracLoading::Get_svp_Values(double& s, double& v, double& p, double time, int timeIndex)
{
	double baseVal = GetBaseValue(time, timeIndex);
	if (directionalBCType == bct_Neumann)
	{
		s = baseVal;
		v = s / Z;
	}
	else if (directionalBCType == bct_Dirichlet)
	{
		v = baseVal;
		s = v * Z;
	}
	else if (directionalBCType == bct_Characteristics)
	{
		s = 0.5 * baseVal;
		v = s / Z;
	}
	p = s * v;
}

void DiracLoading::CalculateValuesAndIntegrals(vector<twsvp>& vals, twsvp& integrals)
{
	unsigned int numSteps = (unsigned int)round(maxTime / timeStepping_delt);
	unsigned int szVals = numSteps + 1;
	vals.resize(szVals);
	double intWeightBase = timeStepping_delt / 3.0;
	double int1 = 0.0;
	for (unsigned int i = 0; i < szVals; ++i)
	{
		twsvp* valPtr = &vals[i];
		valPtr->t = i * timeStepping_delt;
		Get_svp_Values(valPtr->s, valPtr->v, valPtr->p, valPtr->t, i);
		valPtr->w = 4.0 * intWeightBase;
		if (i % 2 == 0)
		{
			if ((i == 0) || (i == numSteps))
				valPtr->w = intWeightBase;
			else
				valPtr->w = 2.0 * intWeightBase;
		}
		int1 += valPtr->w;
		integrals.s += valPtr->w * valPtr->s;
		integrals.v += valPtr->w * valPtr->v;
		integrals.p += valPtr->w * valPtr->p;
	}
	if (fabs(int1 - maxTime) > 1e-4 * maxTime)
	{
		cout << "int1\t" << int1 << '\n';
		cout << "maxTime\t" << maxTime << '\n';
		cout << "error in integration\t";
		THROW("int1 not equal to maxTime\n");
	}
}

bool DiracLoading::DiracLoadingValid(ostream& out, bool printVals)
{
	vector<twsvp> vals; twsvp integrals;
	CalculateValuesAndIntegrals(vals, integrals);
	if (printVals)
	{
		out << "integrals\n" << integrals << '\n';
		unsigned int sz = vals.size();
		out << "vals\tsz\t" << sz << '\n';
		for (unsigned int i = 0; i < sz; ++i)
			out << vals[i] << '\n';
	}
	double intVal = 0.0;
	if (fld4Int1 == di_energy)
		intVal = integrals.p;
	else if (fld4Int1 == di_stress)
		intVal = integrals.s;
	else if (fld4Int1 == di_velocity)
		intVal = integrals.v;
	return (fabs(intVal - 1.0) < 1e-5);
}

void Test_DiracLoading()
{
	double timeStepping_delt = 0.01;
//	double maxTime = 10.0;
	vector <double> maxTimes; // < 0, num Steps, > 0 actual times
	if (true) // use number of steps
	{
		maxTimes.push_back(-2);
		maxTimes.push_back(-4);
		maxTimes.push_back(-6);
		maxTimes.push_back(-10);
	}
	else
	{
		maxTimes.push_back(2.5);
		maxTimes.push_back(3.2);
		maxTimes.push_back(6.0);
		maxTimes.push_back(10.2);
		for (unsigned int i = 0; i < maxTimes.size(); ++i)
			maxTimes[i] *= timeStepping_delt;
	}

	vector<DiracIntegralT> fld4Int1s;
	fld4Int1s.push_back(di_energy);
	fld4Int1s.push_back(di_stress);
	fld4Int1s.push_back(di_velocity);

	double Z = 4.0;

	vector<BoundaryConditionT> directionalBCTypes;
	directionalBCTypes.push_back(bct_Neumann);
	directionalBCTypes.push_back(bct_Dirichlet);
	directionalBCTypes.push_back(bct_Characteristics);

	unsigned int sz_fld4Int1 = fld4Int1s.size();
	unsigned int sz_directionalBCType = directionalBCTypes.size();
	unsigned int sz_maxTime = maxTimes.size();

	fstream out("TestFiles/TestDirac.txt", ios::out);
	unsigned int cntr = 0;
	for (unsigned int fi = 0; fi < sz_fld4Int1; ++fi)
	{
		DiracIntegralT fld4Int1 = fld4Int1s[fi];
		for (unsigned int bi = 0; bi < sz_fld4Int1; ++bi)
		{
			BoundaryConditionT directionalBCType = directionalBCTypes[bi];
			for (unsigned int mi = 0; mi < sz_maxTime; ++mi)
			{
				double maxTime = maxTimes[mi];
				DiracLoading dl;
				dl.timeStepping_delt = timeStepping_delt;
				dl.SetDiractMaxTime(maxTime);
				dl.directionalBCType = directionalBCType;
				dl.fld4Int1 = fld4Int1;
				dl.Z = Z;
				dl.Initialize_DiracLoading();

				out << "cntr_" << cntr << "st_fi_" << fld4Int1 << "_bi_" << directionalBCType << "_mi_" << maxTime << "\n";
				bool valid = dl.DiracLoadingValid(out, true);
				out << "valid" << valid << "\tcntr_" << cntr++ << "en_fi_" << fld4Int1 << "_bi_" << directionalBCType << "_mi_" << maxTime << "\n\n";
			}
		}
	}
}
