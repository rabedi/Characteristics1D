#include "SL_OneInterfaceAllTimes.h"
#include "SLDescriptorData.h"
#include "Domain_AllInterfacesAllTimes.h"

SL_OneInterfaceAllTimes::SL_OneInterfaceAllTimes()
{
	interface_x = 0.0;
	interface_pos = 0;
	interface_flag = 1;
	ts_bulkProps = NULL;
	interfacePFs = NULL;

	ts_bulkProps_Deletable = true;
	interfacePFs_Deletable = true;
	sigmaCFactor = 1.0;
	deltaCFactor = 1.0;
	iniDamage = 0.0;
	outScalars = NULL;
	outFinalSln = NULL;
	outAdaptivity = NULL;
	outIterationConv = NULL;
	b_single_interfaceProblem = true;
	for (int eside = 0; eside < NUM_SIDES; ++eside)
	{
		sides_timeSeqData[eside] = NULL;
		sides_x[eside] = 0.0;
		sides_delx[eside] = 0.0;
		setValue(sides_delts[eside], 0.0);
	}
	sz_subDomainNos = 0;
	iofFinalSolution = iof_ascii;
	iofScalarVals = iof_ascii;
	incidentSide = ilt_noSided;
	only1DOrtizModel = false;
	has_ring_opened1D_al = false;
}

SL_OneInterfaceAllTimes::~SL_OneInterfaceAllTimes()
{
	if (ts_bulkProps_Deletable && (ts_bulkProps != NULL))
		delete ts_bulkProps;
	if (interfacePFs_Deletable && (interfacePFs != NULL))
		delete interfacePFs;

	if (outScalars != NULL)
		delete outScalars;
	if (outFinalSln != NULL)
		delete outFinalSln;
	if (outAdaptivity != NULL)
		delete outAdaptivity;
	if (outIterationConv != NULL)
		delete outIterationConv;
}

SL_OneInterfaceAllTimes::SL_OneInterfaceAllTimes(const SL_OneInterfaceAllTimes& other)
{
	ts_bulkProps = NULL;
	interfacePFs = NULL;
	ts_bulkProps_Deletable = true;
	interfacePFs_Deletable = true;
	outScalars = NULL;
	outFinalSln = NULL;
	outAdaptivity = NULL;
	outIterationConv = NULL;
	(*this) = other;
}

SL_OneInterfaceAllTimes& SL_OneInterfaceAllTimes::operator=(const SL_OneInterfaceAllTimes& other)
{
	interface_pos = other.interface_pos;
	interface_x = other.interface_x;
	interface_flag = other.interface_flag;
	b_single_interfaceProblem = other.b_single_interfaceProblem;
	sigmaCFactor = other.sigmaCFactor;
	deltaCFactor = other.deltaCFactor;
	min_delT = other.min_delT;
	max_delT = other.max_delT;
	iofFinalSolution = other.iofFinalSolution;
	iofScalarVals = other.iofScalarVals;

	if (ts_bulkProps != NULL)
	{
		if (ts_bulkProps_Deletable)
			delete ts_bulkProps;
		ts_bulkProps = NULL;
	}
	if (interfacePFs != NULL)
	{
		if (interfacePFs_Deletable)
			delete interfacePFs;
		interfacePFs = NULL;
	}
	ts_bulkProps_Deletable = other.ts_bulkProps_Deletable;
	interfacePFs_Deletable = other.interfacePFs_Deletable;
	if (other.ts_bulkProps != NULL)
	{
		if (ts_bulkProps_Deletable)
		{
			ts_bulkProps = new SL_Elastic_InterfaceProperties();
			*ts_bulkProps = *other.ts_bulkProps;
		}
		else
			ts_bulkProps = other.ts_bulkProps;
	}
	if (other.interfacePFs != NULL)
	{
		if (interfacePFs_Deletable)
		{
			interfacePFs = new SL_Interface_Fracture_PF();
			*interfacePFs = *other.interfacePFs;
		}
		else
			interfacePFs = other.interfacePFs;
	}
	for (unsigned int eside = 0; eside < NUM_SIDES; ++eside)
	{
		sides_x[eside] = other.sides_x[eside];
		sides_delx[eside] = other.sides_delx[eside];
		CopyVec(other.sides_delts[eside], sides_delts[eside]);

		// this operation will be done at the end ... to ensure all have the right pointers
		sides_timeSeqData[eside] = other.sides_timeSeqData[eside];
#if 0
		if (other.sides_timeSeqData[eside] != NULL)
		{
			sides_timeSeqData[eside] = new SL_interfacePPtData_Time_Seq();
			*(sides_timeSeqData[eside]) = *(other.sides_timeSeqData[eside]);
		}
		else
			sides_timeSeqData[eside] = NULL;
#endif
	}
	timeSeqData = other.timeSeqData;
	if (outScalars != NULL)
	{
		delete outScalars;
		outScalars = NULL;
	}
	if (outFinalSln != NULL)
	{
		delete outFinalSln;
		outFinalSln = NULL;
	}
	if (outAdaptivity != NULL)
	{
		delete outAdaptivity;
		outAdaptivity = NULL;
	}
	if (outIterationConv != NULL)
	{
		delete outIterationConv;
		outIterationConv = NULL;
	}

	if (outFinalSln != NULL)
		Open_fixed_x_files_SL_OneInterfaceAllTimes(iofFinalSolution, iofScalarVals);
	return *this;
}

void SL_OneInterfaceAllTimes::Set_EF_Properties(SL_Interface_Fracture_PF* interfacePFsIn, SL_Elastic_InterfaceProperties* ts_bulkPropsIn)
{
	interfacePFs = interfacePFsIn;
	interfacePFs_Deletable = false;
	if (ts_bulkPropsIn != NULL)
	{
		ts_bulkProps = ts_bulkPropsIn;
		ts_bulkProps_Deletable = false;
	}
}

double SL_OneInterfaceAllTimes::getDeltaC() const
{
	if (interfacePFs != NULL)
		return interfacePFs->deltaC * deltaCFactor;
	return -1.0;
}

double SL_OneInterfaceAllTimes::getSigmaC() const
{
	if (interfacePFs != NULL)
		return interfacePFs->gen_strength * sigmaCFactor;
	return -1.0;
}

void SL_OneInterfaceAllTimes::Set_left_right_TimeSequenceData(InterfaceLocation1DT interfaceLoc, BoundaryConditionT directionalBCType[],
	SL_Bulk_Properties* bulkLeft, SL_Bulk_Properties* bulkRight,
	bool ts_bulkPropsDeletableIn, SL_Elastic_InterfaceProperties*& ts_bulkPropsInOut,
	SL_OneInterfaceAllTimes* left_allTimes, SL_OneInterfaceAllTimes* right_allTimes, double* periodic_totalLength)
{
	if (!ts_bulkPropsDeletableIn) // means that both left and right bulk are homogeneous, so the interface property is storable and reusable (in output interfacePFsInOut)
	{
		if (ts_bulkPropsInOut == NULL) // this thing is not set yet, forming it
		{
			ts_bulkPropsInOut = new SL_Elastic_InterfaceProperties();
			ts_bulkPropsInOut->bulk_leftPtr = bulkLeft;
			ts_bulkPropsInOut->bulk_leftPtr_deletable = false;

			ts_bulkPropsInOut->bulk_rightPtr = bulkRight;
			ts_bulkPropsInOut->bulk_rightPtr_deletable = false;
			ts_bulkPropsInOut->FormMaps();
		}
		// now set interior elastic interface equal to this and mark it as undeletable
		ts_bulkProps = ts_bulkPropsInOut;
	}
	else // either boundary of the domain OR inhomogenous problems
	{
		ts_bulkProps = new SL_Elastic_InterfaceProperties();
		ts_bulkProps->bulk_leftPtr = bulkLeft;
		ts_bulkProps->bulk_leftPtr_deletable = false;

		ts_bulkProps->bulk_rightPtr = bulkRight;
		ts_bulkProps->bulk_rightPtr_deletable = false;
		if (directionalBCType != NULL)
		{
			for (unsigned int i = 0; i < DiM; ++i)
				ts_bulkProps->directionalBCType[i] = directionalBCType[i];
		}
		else
		{
			for (unsigned int i = 0; i < DiM; ++i)
				ts_bulkProps->directionalBCType[i] = bct_Unspecified;
		}
		ts_bulkProps->FormMaps();
	}
	ts_bulkProps_Deletable = ts_bulkPropsDeletableIn;

	// now going over sides form time sequence, x, delx, delt
	min_delT = DBL_MAX;
	max_delT = 0.0;

	if ((left_allTimes != NULL) && (ts_bulkProps->bulk_leftPtr != NULL))
	{
		sides_timeSeqData[SDL] = &left_allTimes->timeSeqData;
		b_single_interfaceProblem = false;
		sides_x[SDL] = left_allTimes->interface_x;
		sides_delx[SDL] = interface_x - left_allTimes->interface_x;
		if ((periodic_totalLength != NULL) && (sides_delx[SDL] < 0))
			sides_delx[SDL] += *periodic_totalLength;
		for (int i = 0; i < DiM; ++i)
			sides_delts[SDL][i] = sides_delx[SDL] / ts_bulkProps->bulk_leftPtr->ws[i];
		min_delT = MIN(min_delT, sides_delts[SDL][0]);
		max_delT = MAX(max_delT, sides_delts[SDL][0]);
	}
	if ((right_allTimes != NULL) && (ts_bulkProps->bulk_rightPtr != NULL))
	{
		sides_timeSeqData[SDR] = &right_allTimes->timeSeqData;
		b_single_interfaceProblem = false;
		sides_x[SDR] = right_allTimes->interface_x;
		sides_delx[SDR] = right_allTimes->interface_x - interface_x;
		if ((periodic_totalLength != NULL) && (sides_delx[SDR] < 0))
			sides_delx[SDR] += *periodic_totalLength;
		for (int i = 0; i < DiM; ++i)
			sides_delts[SDR][i] = sides_delx[SDR] / ts_bulkProps->bulk_rightPtr->ws[i];
		min_delT = MIN(min_delT, sides_delts[SDR][0]);
		max_delT = MAX(max_delT, sides_delts[SDR][0]);
	}
}

void SL_OneInterfaceAllTimes::Open_fixed_x_files_SL_OneInterfaceAllTimes(IOF_type iofFinalSolutionIn, IOF_type iofScalarValsIn, int subdomainLeft, int subdomainRight)
{
	iofFinalSolution = iofFinalSolutionIn;
	iofScalarVals = iofScalarValsIn;

//	do_mkdir(nameBase.c_str());
//	nameBase += "/";
	bool b_directlyCalculateBondedInFinalSln = ((ts_bulkProps->interfaceLoc != ilt_twoSided) || (!interfacePFs->HasAnyContactFriction()));
#if HAS_SOURCE
	b_directlyCalculateBondedInFinalSln = false;
#endif
	string ext, outName;
	bool genFile = getExt(iofFinalSolution, ext);
	if (genFile)
	{
		string specificName = "InterfaceRawFinalSln";
		if (subdomainLeft < 0)
			GetSubdomainIndexed_SpaceIndexed_FileName(outName, interface_pos, specificName, ext);
		else
			GetSubdomainIndexed_SpaceIndexed_FileName(outName, subdomainLeft, specificName, ext, subdomainRight);
		if (iofFinalSolution == iof_ascii)
			outFinalSln = new fstream(outName.c_str(), ios::out);
		else
			outFinalSln = new fstream(outName.c_str(), ios::out | ios::binary);
	}
	if (!b_directlyCalculateBondedInFinalSln)
	{
		if (g_slf_conf->print_Scalars)
		{
			bool genFile = getExt(iofScalarVals, ext);
			if (genFile)
			{
				string specificName = "InterfaceRawScalars";
				GetSubdomainIndexed_SpaceIndexed_FileName(outName, interface_pos, specificName, ext);
				if (iofScalarVals == iof_ascii)
					outScalars = new fstream(outName.c_str(), ios::out);
				else
					outScalars = new fstream(outName.c_str(), ios::out | ios::binary);
			}
		}
		if (g_slf_conf->print_outIterationConv)
		{
			string specificName = "InterationConv";
			GetSubdomainIndexed_SpaceIndexed_FileName(outName, interface_pos, specificName, ext);
			outIterationConv = new fstream(outName.c_str(), ios::out);
		}
	}
	if (g_slf_conf->print_Adaptivity)
	{
		string specificName = "Adaptivity";
		GetSubdomainIndexed_SpaceIndexed_FileName(outName, interface_pos, specificName, ext);
		outAdaptivity = new fstream(outName.c_str(), ios::out);
	}
}

void SL_OneInterfaceAllTimes::Compute_DownStream_Characteristics_wo_in_situ(VEC& wlSide_rGoing_WO_in_situ, VEC& wrSide_lGoing_WO_in_situ,
	double currentTime, int currentTimeIndex, const SL_interface_Temp_PPtData* current_ptData)
{
	// much simpler way of calculating impinging characteristics
	if (b_single_interfaceProblem)
	{
		g_SL_desc_data.Get1InterfaceImpingingCharacteristics(ts_bulkProps, currentTime, wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ);
		return;
	}
	setValue(wlSide_rGoing_WO_in_situ, 0.0);
	setValue(wrSide_lGoing_WO_in_situ, 0.0);
	VEC *wDownstream;
	vector<int> esideMissing_wCals;

	for (int eside = 0; eside < NUM_SIDES; ++eside)
	{
#if RING_PROBLEM
		double v_r;
#endif
		bool is_wr;
		/////////////////////////////////////////////////////////////////
		// A. Upstream values
		SL_Bulk_Properties* bulk_Ptr;
		if (sides_timeSeqData[eside] == NULL)
		{
			esideMissing_wCals.push_back(eside);
			continue;
		}
		int epside;
		double* ring_opened1D_alPtr = NULL;
		if (eside == SDL)
		{
			bulk_Ptr = ts_bulkProps->bulk_leftPtr;
			wDownstream = &wlSide_rGoing_WO_in_situ;
			epside = SDR;
			is_wr = true;
		}
		else
		{
			bulk_Ptr = ts_bulkProps->bulk_rightPtr;
			wDownstream = &wrSide_lGoing_WO_in_situ;
			epside = SDL;
			is_wr = false;
			if (has_ring_opened1D_al)
				ring_opened1D_alPtr = &g_domain->ring_opened1D_al;
		}

#if HAVE_SOURCE
#if RING_PROBLEM
		double rho, E;
		rho = bulk_Ptr->rho;
		E = bulk_Ptr->E_iso;
#endif
#endif
		AllPoints_Data_Upstream_w_lOr_r upstream_pts;
		for (int i = 0; i < DiM; ++i)
		{
			double tAbsolute = currentTime - sides_delts[eside][i];
			if (tAbsolute > 0) // hit the neighbor time sequence
			{
				upstream_pts.up_Pts[i].delT = sides_delts[eside][i];
				SL_interfacePPtData* ptSlnPtr;
				bool ptSlnPtr_Deletable;
				bool found = sides_timeSeqData[eside]->Interpolate_Pt_Solution_At_time(tAbsolute, ptSlnPtr, ptSlnPtr_Deletable);
				if (found == false)
				{
					THROW("Cannot find the point\n");
				}

				CopyVec(ptSlnPtr->sl_side_ptData[epside].sigma_downstream_final, upstream_pts.up_Pts[i].sigma);
				CopyVec(ptSlnPtr->sl_side_ptData[epside].v_downstream_final, upstream_pts.up_Pts[i].v);
#if HAVE_SOURCE
				g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(sides_x[eside], tAbsolute, bulk_Ptr->flag, upstream_pts.up_Pts[i].source_v, upstream_pts.up_Pts[i].source_sigma);
#if RING_PROBLEM
				v_r = ptSlnPtr->v_r_final;
#endif
#endif
				if (ptSlnPtr_Deletable)
					delete ptSlnPtr;
			}
			else // characteristics hits IC
			{
				double delT = currentTime;
				double del_x = bulk_Ptr->ws[i] * delT;
				double x;
				if (eside == SDL)
				{
					x = interface_x - del_x;
					if (g_domain->isPeriodic && (x < g_domain->x_min))
						x += g_domain->L;

//					double ratioFar = del_x / sides_delx[SDL];
//					double omratioFar = 1.0 - ratioFar;
				}
				else
				{
					x = interface_x + del_x;
					if (g_domain->isPeriodic && (x > g_domain->x_max))
						x -= g_domain->L;

//					double ratioFar = del_x / sides_delx[SDR];
//					double omratioFar = 1.0 - ratioFar;
				}
				upstream_pts.up_Pts[i].delT = delT;

				// getting ICs
				VEC u_0, eps_0;
				g_SL_desc_data.Get_IC(bulk_Ptr, bulk_Ptr->flag, x, u_0, upstream_pts.up_Pts[i].v, eps_0
#if RING_PROBLEM
					, v_r
#endif
				);
				bulk_Ptr->Compute_Stress_from_Strain(eps_0, upstream_pts.up_Pts[i].sigma);

#if HAVE_SOURCE
				g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(x, 0.0, bulk_Ptr->flag, upstream_pts.up_Pts[i].source_v, upstream_pts.up_Pts[i].source_sigma);
#endif
			}
#if HAVE_SOURCE
#if RING_PROBLEM
			double sigma_theta_source_final = bulk_Ptr->E_iso * v_r / g_domain->ring_R; // E v_r / R (same on both sides)
			upstream_pts.up_Pts[i].source_sigma[0] += sigma_theta_source_final;
#endif
#endif
		}
		/////////////////////////////////////////////////////////////////
		// B. Downstream values, parts needed
		OnePoint_Data_Downstream_w_lOr_r downstream_pt;
		if (current_ptData != NULL)
		{
			CopyVec(current_ptData->sl_side_temp_ptData[eside].v_downstream_latestValue, downstream_pt.last_estimate_v);
			CopyVec(current_ptData->sl_side_temp_ptData[eside].sigma_downstream_latestValue, downstream_pt.last_estimate_sigma);
#if RING_PROBLEM
			v_r = current_ptData->v_r_latestValue;
#endif
		}
		else
		{
			// get the latest permanent solution
			SL_interfacePPtData* lastPt = timeSeqData.GetBackwardPosition(0);
			CopyVec(lastPt->sl_side_ptData[eside].v_downstream_final, downstream_pt.last_estimate_v);
			CopyVec(lastPt->sl_side_ptData[eside].sigma_downstream_final, downstream_pt.last_estimate_sigma);
#if RING_PROBLEM
			v_r = lastPt->v_r_final;
#endif
		}
#if HAVE_SOURCE
		g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(interface_x, currentTime, bulk_Ptr->flag, downstream_pt.source_v, downstream_pt.source_sigma);
#if RING_PROBLEM
		double sigma_theta_source_final = bulk_Ptr->E_iso * v_r / g_domain->ring_R; // E v_r / R (same on both sides)
		downstream_pt.source_sigma[0] += sigma_theta_source_final;
#endif
#endif
		bulk_Ptr->Compute_Downstream_characteristic_From_Upstream_qs_and_Forces(upstream_pts, downstream_pt, is_wr, interface_x, false, ring_opened1D_alPtr);
		CopyVec(downstream_pt.wDownstream, *wDownstream);
	}

	if (ts_bulkProps->interfaceLoc != ilt_twoSided)
	{
		unsigned int szMissingCalc = esideMissing_wCals.size();
		if (szMissingCalc == 1)
		{
			if (ts_bulkProps->isIsoInterface == false)
				THROW("(anti)symmetric BC, Dirichlet, Neumann, etc are only supported for isotropic bulk interfaces\n");

			unsigned int miss_side = esideMissing_wCals[0];
			BoundaryConditionT directionalBCt = ts_bulkProps->directionalBCType[0];
			if ((directionalBCt == bct_Symmetric) || (directionalBCt == bct_AntiSymmetric))
			{
				if (miss_side == 0)
				{
					for (unsigned int d = 0; d < DiM; ++d)
					{
						directionalBCt = ts_bulkProps->directionalBCType[d];
						if (directionalBCt == bct_Symmetric)
							wlSide_rGoing_WO_in_situ[d] = wrSide_lGoing_WO_in_situ[d];
						else if (directionalBCt == bct_AntiSymmetric)
							wlSide_rGoing_WO_in_situ[d] = -wrSide_lGoing_WO_in_situ[d];
						else
						{
							cout << "directionalBCt\t" << directionalBCt << '\n';
							THROW("Invalid directionalBCt\n");
						}
					}
				}
				else // if (miss_side == 1)
				{
					for (unsigned int d = 0; d < DiM; ++d)
					{
						directionalBCt = ts_bulkProps->directionalBCType[d];
						if (directionalBCt == bct_Symmetric)
							wrSide_lGoing_WO_in_situ[d] = wlSide_rGoing_WO_in_situ[d];
						else if (directionalBCt == bct_AntiSymmetric)
							wrSide_lGoing_WO_in_situ[d] = -wlSide_rGoing_WO_in_situ[d];
						else
						{
							cout << "directionalBCt\t" << directionalBCt << '\n';
							THROW("Invalid directionalBCt\n");
						}
					}
				}
			}
			else // normal Dirichlet, Neumann BC
			{
				// wlSide_rGoing_WO_in_situ and wrSide_lGoing_WO_in_situ are in fact "loads" and in function SL_Elastic_InterfaceProperties::Compute_sigmaI_vI_from_wsAndBC__BoundaryCase they are interpreted correctly
				if (miss_side == 0)
					g_SL_desc_data.GetBoundaryConditionValue_LeftRight(true, ts_bulkProps, interface_x, currentTime, wlSide_rGoing_WO_in_situ, currentTimeIndex);
				else //if (miss_side == 1)
					g_SL_desc_data.GetBoundaryConditionValue_LeftRight(false, ts_bulkProps, interface_x, currentTime, wrSide_lGoing_WO_in_situ, currentTimeIndex);
			}
		}
		else if (szMissingCalc != 0)
		{
			cout << szMissingCalc << '\n';
			THROW("Invalid domain problem for fracture problem\n");
		}
	}
}

double SL_OneInterfaceAllTimes::getEffectiveDamage_4_InterfacialDamage_TSRs(const SL_interfacePPtData& pt)
{
	if ((interfacePFs != NULL) && (interfacePFs->damageOffOnMix == sl_interfacial_damage_on) && (interfacePFs->damageTractionModel == sl_interfacial_damage_traction_TSR))
		return pt.maxEffDelU / getDeltaC();
	return pt.interface_damage_final;
}

void SL_OneInterfaceAllTimes::Set_ring_opened1D_left_side_jump_handling_true()
{
	has_ring_opened1D_al = true;
}

void SL_OneInterfaceAllTimes::Set1DOrtizType()
{
#if DiM1
	only1DOrtizModel = (interfacePFs != NULL) && (interfacePFs->IsSimpleOrtizModel());
#else
	only1DOrtizModel = false;
#endif
}

void SL_OneInterfaceAllTimes::Print(ostream& out) const
{
	out << "interface_flag\t" << interface_flag;
	out << "\tinterface_pos\t" << interface_pos;
	out << "\tinterface_x\t" << interface_x;
	out << "\tmin_delT\t" << min_delT;
	out << "\tmax_delT\t" << max_delT;
	out << "\tsigmaCFactor\t" << sigmaCFactor;
	out << "\tdeltaCFactor\t" << deltaCFactor;
	out << "\tiniDamage\t" << iniDamage;
	out << "\tonly1DOrtizModel\t" << only1DOrtizModel;
	out << "\tb_single_interfaceProblem\t" << b_single_interfaceProblem;

	out << "\nincidentSide\t" << incidentSide;
	out << "\tsz_subDomainNos\t" << sz_subDomainNos;
	out << "\tsubDomainNos";	for (unsigned int pos = 0; pos < subDomainNos.size(); ++pos) out << "\t" << subDomainNos[pos]; out << '\n';
	out << "\relPos_wrt_subDomainStartPoints";	for (unsigned int pos = 0; pos < relPos_wrt_subDomainStartPoints.size(); ++pos) out << "\t" << relPos_wrt_subDomainStartPoints[pos]; out << '\n';

	out << "\ntimeSeqData\t" << &timeSeqData;
	for (unsigned int sidei = 0; sidei < NUM_SIDES; ++sidei)
	{
		out << "\nside" << sidei;
		out << "\tsides_timeSeqData\t" << sides_timeSeqData[sidei];
		out << "\tsides_bulk_index\t" << sides_bulk_index[sidei];
		out << "\tsides_x\t" << sides_x[sidei];
		out << "\tsides_delx\t" << sides_delx[sidei];
		out << "\tsides_delts\t" << sides_delts[sidei][0];
	}

	out << "\ninterfacePFs_Deletable\t" << interfacePFs_Deletable;
	out << "\tinterfacePFs\t" << interfacePFs;
	if (interfacePFs != NULL)
	{
		out << "\n";
		interfacePFs->PrintShort(out);
	}
	out << "\n";

	out << "\nts_bulkProps_Deletable\t" << ts_bulkProps_Deletable;
	out << "\ts_bulkProps\t" << ts_bulkProps;
	if (ts_bulkProps != NULL)
	{
		out << "\n";
		ts_bulkProps->Print(out, true);
	}
	out << "\n";
}

void SL_OneInterfaceAllTimes::InitialStep_Use_IC(SLInterfaceCalculator& slic)
{
	// u, v, sig must be continuous across the interface. // eps may not be (if stiffness changes)
	VEC u_0, v_0, eps_0;
	VEC sig_0;
	double damage_source_0 = 0.0;
	double maxEffDelU = 0.0;
#if RING_PROBLEM
	double v_r;
#endif
	if (ts_bulkProps->bulk_leftPtr != NULL)
	{
		g_SL_desc_data.Get_IC(ts_bulkProps->bulk_leftPtr, ts_bulkProps->bulk_leftPtr->flag, interface_x, u_0, v_0, eps_0
#if RING_PROBLEM
			, v_r
#endif
		);
		ts_bulkProps->bulk_leftPtr->Compute_Stress_from_Strain(eps_0, sig_0);
		// note: no check is done whether stress is continuous from two side calculation. This is left to the user
		if (ts_bulkProps->bulk_rightPtr != NULL)
		{
			if ((outScalars != NULL) && (interfacePFs->HasAnyContactFriction()))
			{
				double scalar_Vals[s1i_SIZE];
				interfacePFs->GetIC_ScalarValues(sig_0, sigmaCFactor, deltaCFactor, iniDamage, maxEffDelU, damage_source_0, scalar_Vals);
				Output_ScalarValuesAux(*outScalars, iofScalarVals, 0.0, scalar_Vals, interfacePFs->damageOffOnMix, interfacePFs->contactOffOnMix, interfacePFs->slipOffOnMix);

				SL_interface_Temp_PPtData* tPtSlns = new SL_interface_Temp_PPtData();
				slic.Set_tPtSlns(tPtSlns, true);
				slic.tPtSlns->interfacePropPtr = new SL_Interface_Temp_PtData_IntContFrac();
				for (unsigned int i = 0; i < s1i_SIZE; ++i)
					slic.tPtSlns->interfacePropPtr->interface_scalar_Vals[i] = scalar_Vals[i];
				slic.tPtSlns->damageOffOnMix = interfacePFs->damageOffOnMix;
				slic.tPtSlns->contactOffOnMix = interfacePFs->contactOffOnMix;
				slic.tPtSlns->slipOffOnMix = interfacePFs->slipOffOnMix;
			}
		}
	}
	else
	{
		g_SL_desc_data.Get_IC(ts_bulkProps->bulk_rightPtr, ts_bulkProps->bulk_rightPtr->flag, interface_x, u_0, v_0, eps_0
#if RING_PROBLEM
			, v_r
#endif
		);
		ts_bulkProps->bulk_rightPtr->Compute_Stress_from_Strain(eps_0, sig_0);
	}

	// getting a new point
	unsigned int pos;
	SL_interfacePPtData* ppt = timeSeqData.AddNewPoint(0.0, pos);
	ppt->maxEffDelU = maxEffDelU;
	ppt->tsrStage = TSR_loadingStages_none;
	if ((only1DOrtizModel) || ((interfacePFs != NULL) && (interfacePFs->damageOffOnMix == sl_interfacial_damage_on)))
		iniDamage = 1.0;
	ppt->interface_damage_final = iniDamage;
	ppt->interface_damage_source_final = damage_source_0;
	ppt->interface_time = 0.0;

//	if (ts_bulkProps->bulk_leftPtr != NULL)
	{
		CopyVec(u_0, ppt->sl_side_ptData[SDL].u_downstream_final);
		CopyVec(v_0, ppt->sl_side_ptData[SDL].v_downstream_final);
		CopyVec(sig_0, ppt->sl_side_ptData[SDL].sigma_downstream_final);
	}
//	if (ts_bulkProps->bulk_rightPtr != NULL)
	{
		CopyVec(u_0, ppt->sl_side_ptData[SDR].u_downstream_final);
		CopyVec(v_0, ppt->sl_side_ptData[SDR].v_downstream_final);
		CopyVec(sig_0, ppt->sl_side_ptData[SDR].sigma_downstream_final);
	}
	if (has_ring_opened1D_al)
		ppt->sl_side_ptData[SDR].v_downstream_final[0] -= g_domain->ring_opened1D_al;
#if RING_PROBLEM
	double invrho, sigma_theta, E;
	sigma_theta = sig_0[0];
	invrho = 0.5 * (1.0 / ts_bulkProps->bulk_leftPtr->rho + 1.0 / ts_bulkProps->bulk_rightPtr->rho);
	E = 0.5 * (ts_bulkProps->bulk_leftPtr->E_iso + ts_bulkProps->bulk_rightPtr->E_iso);
	ppt->v_r_final = v_r;
	double ring_p = g_SL_desc_data.GetRing_p(interface_x, 0.0);
	ppt->v_r_source_final = ring_p - sigma_theta * invrho / g_domain->ring_R; // p - sigma_theta/rho/R
	ppt->sigma_theta_source_final = E * v_r / g_domain->ring_R; // E v_r / R (same on both sides)
#endif
	// handling incident interface
	if (incidentSide != ilt_noSided) 
	{
		double s, stressBar = g_SL_desc_data.tdLoad_stressScale * g_SL_desc_data.tdLoadComputer->ComputeSpecific_TD_Value(0.0);
		ppt->sl_side_ptData[SDL].SetZero();
		ppt->sl_side_ptData[SDR].SetZero();
		double Zl = ts_bulkProps->bulk_leftPtr->c_rhos[0], Zr = ts_bulkProps->bulk_rightPtr->c_rhos[0];
		double ZsumInv = 1.0 / (Zl + Zr);
		if (g_SL_desc_data.tdLoadSide == SDL)
		{
			s = (Zr - Zl) * ZsumInv * stressBar;
			ppt->sl_side_ptData[SDL].sigma_downstream_final[0] = s;
			ppt->sl_side_ptData[SDL].v_downstream_final[0] = s / Zl;
			s = 2.0 * Zr * ZsumInv * stressBar;
			ppt->sl_side_ptData[SDR].sigma_downstream_final[0] = s;
			ppt->sl_side_ptData[SDR].v_downstream_final[0] = -s / Zr;
		}
		else
		{
			s = (Zl - Zr) * ZsumInv * stressBar;
			ppt->sl_side_ptData[SDR].sigma_downstream_final[0] = s;
			ppt->sl_side_ptData[SDR].v_downstream_final[0] = -s / Zr;
			s = 2.0 * Zl * ZsumInv * stressBar;
			ppt->sl_side_ptData[SDL].sigma_downstream_final[0] = s;
			ppt->sl_side_ptData[SDL].v_downstream_final[0] = s / Zl;
		}
	}
	slic.Set_pPtSlns(ppt, false);
	if (outFinalSln != NULL)
		ppt->Output_FinalSolution(*outFinalSln, iofFinalSolution, 0.0, interface_x, 0.0, has_ring_opened1D_al);
}

AdaptivityS SL_OneInterfaceAllTimes::NonInitialStep(double deltaT, bool &accept_point, double& maxTime, int currentTimeIndex, SLInterfaceCalculator& slic)
{
	double maxTimeBeforeThisPoint = timeSeqData.GetMaxTime();
	double currentTime = maxTimeBeforeThisPoint + deltaT;
	// getting a new point

	unsigned int pos;
	SL_interfacePPtData* ppt = timeSeqData.AddNewPoint(currentTime, pos);

	slic.current_timeIndex = currentTimeIndex;
	slic.oneIntAlltimesPtr = this;
	// adding elastic and fracture pointers
	slic.Initialize_setFromOutside_ElsticFractureProperties(ts_bulkProps, interfacePFs);
	// adding permanent point
	slic.Initialize_setFromOutside_Current_Step_EmptyPermanentPt_Storage(ppt); // , interface_x);
	// adding the sequence of earlier solutions
	slic.Initialize_setFromOutside_Earlier_Steps_NonEmptyPermanentPts_Storage(timeSeqData);

	// computing downstream characteristics and passing on to calculator.
	// note: for nonzero source problems characteristics are updated in SLInterfaceCalculator::Main_Compute_OnePoint
	VEC wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ;
	SL_interface_Temp_PPtData* current_ptData = NULL;
	bool compute_w_from_IC = true;
	double t0;
	if (g_slf_conf->isPeriodic_1Fragment)
	{
		double ws = slic.ts_bulkProps->bulk_leftPtr->cd_iso;
		double L = g_slf_conf->periodic_1Fragment_size;
		t0 = L / ws;
		if (currentTime > t0)
		{
			compute_w_from_IC = false;
			double timePrev = currentTime - t0;
			double vL, vR, sigma;
			bool found = g_seq_short.GetPt(timePrev, vL, vR, sigma);
			if (!found)
			{
				THROW("Previous point cannot be found\n");
			}
			double Z = slic.ts_bulkProps->bulk_leftPtr->c_rhos[0];
			double a = g_SL_desc_data.a_xt_prob[0];
			double Zla = Z * L * a;
			wlSide_rGoing_WO_in_situ[0] = sigma - Z * vR + Zla;
			wrSide_lGoing_WO_in_situ[0] = sigma + Z * vL + Zla;
		}
	}
	if (compute_w_from_IC)
		Compute_DownStream_Characteristics_wo_in_situ(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ,
		currentTime, slic.current_timeIndex, current_ptData);
	slic.Initialize_setFromOutside_Incoming_Characteristics_etc(&wlSide_rGoing_WO_in_situ, &wrSide_lGoing_WO_in_situ);
	// setting sigmaC and deltaC factors
	slic.Initialize_setFromOutside_Fracture_InhomogeneousFactors(sigmaCFactor, deltaCFactor);

	// adding incident part
	slic.incidentSide = incidentSide;
	if (incidentSide != ilt_noSided)
	{
		setValue(slic.incientIncomingStress, 0.0);
		slic.incientIncomingStress[0] = g_SL_desc_data.tdLoad_stressScale * g_SL_desc_data.tdLoadComputer->ComputeSpecific_TD_Value(currentTime);
	}

	// now doing the adaptivity
	AdaptivityS as = slic.Main_Compute_OnePoint(accept_point, iofFinalSolution, iofScalarVals, currentTime, outScalars, outFinalSln, outAdaptivity, outIterationConv);
	// for refine flag, the point is removed
	if (!accept_point)
	{
		timeSeqData.RemoveLastPoint();
		maxTime = maxTimeBeforeThisPoint;
	}
	else
	{
		maxTime = currentTime;
		if (g_slf_conf->isPeriodic_1Fragment)
		{
			double sigma = slic.pPtSlns->sl_side_ptData[SDL].sigma_downstream_final[0];
			double vL = slic.pPtSlns->sl_side_ptData[SDL].v_downstream_final[0], vR = slic.pPtSlns->sl_side_ptData[SDR].v_downstream_final[0];
			g_seq_short.AddPt(maxTime, vL, vR, sigma, 2.0 * t0);
			double delv = vR - vL;
			double uL = slic.pPtSlns->sl_side_ptData[SDL].u_downstream_final[0], uR = slic.pPtSlns->sl_side_ptData[SDR].u_downstream_final[0];
			double delu = uR - uL;
			PITSS terminateState = per_if.UpdateStats(maxTime, delu, delv, sigma);
		}
	}
	return as;
}
int SL_OneInterfaceAllTimes::Main_One_InterfaceProblem()
{
	int terminateFlag = Main_One_InterfaceProblem_Aux();
	if (per_if.isActive)
	{
		string fileName = per_if.nameOne_a_SharedAll_l +  "_laStudies.txt";
		fstream in(fileName.c_str(), ios::in);
		bool fileExists = in.is_open();
		if (fileExists)
			in.close();
		fstream out(fileName.c_str(), ios::app);
		if (!fileExists)
			per_if.Output_Periodic1IntrFrag_Header(out);
		per_if.Output_Periodic1IntrFrag(out);
	}
	return terminateFlag;
}

int SL_OneInterfaceAllTimes::Main_One_InterfaceProblem_Aux()
{
	SLInterfaceCalculator slic;
	InitialStep_Use_IC(slic);
	bool cont_time_steps = true;
	iofFinalSolution = iof_ascii;
	iofScalarVals = iof_ascii;
	AdaptivityS as;
	as.Update_AdaptFlag(a_unassigned);
	as.a_delt = g_slf_conf->uniform_del_t;
	bool accept_point;
	double maxTime;
	Set1DOrtizType();
	if (g_slf_conf->isPeriodic_1Fragment)
		g_seq_short.AddPt(0.0, 0.0, 0.0, 1.0, 2.0 * per_if.aInv);

	long cntr = 0;
	while (cntr < 1000000000000)
	{
		if (cntr++ % 1000 == 0)
			cout << cntr << '\n';
		SLInterfaceCalculator slic;
		as = NonInitialStep(as.a_delt, accept_point, maxTime, cntr, slic);
		if (per_if.isActive)
		{
			double maxDelU = slic.pPtSlns->maxEffDelU;
			double delta_u = slic.pPtSlns->sl_side_ptData[SDR].u_downstream_final[0] - slic.pPtSlns->sl_side_ptData[SDL].u_downstream_final[0];
			bool ready2Return = false;
#if 0
			if (maxTime >= per_if.t_dilute_zeroStressCheck)
			{
				if (per_if.diluteFractureModel)
				{
					per_if.terminateState = pit_sigma0;
					return 1;
				}
				cout << "delta_u\t" << delta_u << '\n';
				if (per_if.delu4Sigma0 - delta_u < 1e-6)
				{
					per_if.terminateState = pit_sigma0;
					return 1;
				}
			}
			else 
#endif
			if (per_if.terminateState != PITSS_none)
				return 1;
		}
		AdaptivityF a_flag = as.get_a_flag();
		if ((maxTime >= g_slf_conf->terminate_run_target_time) || (a_flag == a_terminate_run_correctly))
		{
			if (!per_if.isActive)
				return 1;
		}
		if (a_flag == a_terminate_run_prematurely)
			return 0;
	}
	return -1;
}

void MAIN_SL_OneInterfaceAllTimes_ONE_Interface(string configNameIn)
{
	string icbc_configName;	// source term, ic, bc, ...
	string slfg_configName; // SLFractureGlobal_Configuration, g_slf_conf
	string pf_configName;	// SL_Interface_Fracture_PF
	string el_configName;   // SL_Elastic_InterfaceProperties
	fstream in(configNameIn.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "config file name:\t" << configNameIn << '\n';
		THROW("file cannot be opened\n");
	}
	string buf;
	READ_NSTRING(in, buf, icbc_configName);
	READ_NSTRING(in, buf, slfg_configName);
	READ_NSTRING(in, buf, pf_configName);
	READ_NSTRING(in, buf, el_configName);

	if (icbc_configName != "default")
		g_SL_desc_data.Read(icbc_configName);

	fstream inb(pf_configName.c_str(), ios::in);
	if (!inb.is_open())
	{
		cout << "pf_configName\t" << pf_configName << '\n';
		THROW("File cannot be opened\n");
	}
	SL_Interface_Fracture_PF interfacePFs;
	interfacePFs.Read_SL_Interface_Fracture_PF(inb, 1);
	inb.close();

	if (slfg_configName != "default")
		g_slf_conf->Read(slfg_configName);
	else
		g_slf_conf->Initialize_SLFractureGlobal_Configuration_After_Reading();

	if (g_slf_conf->isPeriodic_1Fragment)
	{
		per_if.isActive = true;
		per_if.l = g_slf_conf->periodic_1Fragment_size;
		per_if.a = g_SL_desc_data.a_xt_prob[0];
		per_if.isExtrinsic = IsExtrinsic(interfacePFs.tsrModel);
		per_if.energyScale = GetEnergyConstantFactor(interfacePFs.tsrModel);
		if (!per_if.isExtrinsic)
			per_if.t0 = 0.0;
		else
			per_if.t0 = 1.0 / per_if.a;

		gt0 = per_if.t0;
		per_if.sigma_t0 = per_if.t0 * per_if.a;
		if (g_slf_conf->uniform_del_t < 0)
			per_if.relTol = -g_slf_conf->uniform_del_t;
		per_if.Initialize_Periodic1IntrFrag();
		g_slf_conf->uniform_del_t = per_if.suggeste_delt;
		g_slf_conf->terminate_run_target_time = per_if.suggested_final_time4ZeroSigmaC;
		g_slf_conf->between_steps_adaptivity = false;
	}
	inb.open(el_configName.c_str(), ios::in);
	if (!inb.is_open())
	{
		cout << "el_configName\t" << el_configName << '\n';
		THROW("File cannot be opened\n");
	}
	SL_Elastic_InterfaceProperties ts_bulkProps;
	ts_bulkProps.Read(inb);
	inb.close();


	SL_OneInterfaceAllTimes oiat;
	oiat.interfacePFs = &interfacePFs;
	oiat.Set_EF_Properties(&interfacePFs, &ts_bulkProps);
	oiat.Open_fixed_x_files_SL_OneInterfaceAllTimes(iof_ascii, iof_ascii);

	int successMode = oiat.Main_One_InterfaceProblem();
	cout << "successMode\t" << successMode << '\n';
}

