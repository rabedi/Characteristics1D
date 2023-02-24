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


	tdLoad_Zi = 0.0, tdLoad_Zo = 0.0, tdLoad_ZEffective = 0.0;
	tdLoad_inputEnergy = 0.0, tdLoad_stressScale = 0.0, tdLoad_velScale = 0.0, tdLoad_loadScale = 0.0, tdAmbientProjectileLength = 0.0;
}

SLDescriptorData::~SLDescriptorData()
{
	if (tdLoadComputer != NULL)
		delete tdLoadComputer;
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
		double value; map<string, string>* mpPtr;
		if (Find_Version_Value("la", value, mpPtr))
			load_parameters[0] = pow(10.0, value);
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
}

bool SLDescriptorData::GetNonRingNonZeroTerm_SourceTerm(double x, double t, GID bulkFlag, VEC & source_v, VEC & source_sigma) const
{
	if ((load_number == AXT_LN) || (tdLoadComputer != NULL))
		return false;
	THROW("Non-axt problem is not implemented.\n");
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
			v_0[i] = a_xt_prob[i] * x;
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

void SLDescriptorData::GetBoundaryConditionValue_LeftRight(bool isLeft, SL_Elastic_InterfaceProperties* ts_bulkProps, double x, double t, VEC& BC_val) const
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
					BC_val[i] = a_xt_prob[i] * t * bulk_Ptr->effective_CsIso[i];
				else if (ts_bulkProps->directionalBCType[i] == bct_Characteristics)
				{
					if (isLeft)
						BC_val[i] = a_xt_prob[i] * t * bulk_Ptr->effective_CsIso[i] - a_xt_prob[i] * x;
					else
						BC_val[i] = a_xt_prob[i] * t * bulk_Ptr->effective_CsIso[i] + a_xt_prob[i] * x;
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
				wlSide_rGoing[i] = a_xt_prob[i] * t * ts_bulkProps->bulk_leftPtr->effective_CsIso[i];
				wrSide_lGoing[i] = a_xt_prob[i] * t * ts_bulkProps->bulk_rightPtr->effective_CsIso[i];
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
			tdLoad_inputEnergy = tdLoad_stressScale * tdLoad_stressScale * integral_load2Out / tdLoad_ZEffective;
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
			tdLoad_inputEnergy = 0.5 * rhoout * tdLoad_velScale * tdLoad_velScale * tdAmbientProjectileLength;
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
		tdLoad_inputEnergy = tdLoad_stressScale * tdLoad_stressScale * integral_load2Out / tdLoad_ZEffective;
		tdAmbientProjectileLength = 0.0;
	}
	Print_tdLoad();
}

void SLDescriptorData::Print_tdLoad()
{
	string fileName = g_prefileName + "/" + "__tdLoad.txt";
	fstream out(fileName.c_str(), ios::out);

	out << "tdLoad_inputEnergy\t" << tdLoad_inputEnergy << '\n';
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
	in >> buf >> tdLoad_inputEnergy;
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
