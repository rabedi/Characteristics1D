#include "SLInterfaceCalculator.h"
#include "SLDescriptorData.h"
#include "SL_OneInterfaceAllTimes.h"
#include "Domain_AllInterfacesAllTimes.h"
#include "globalFunctions.h"
#include <limits>
#include <cmath>
using namespace std;

Periodic1IntrFrag* per_if = NULL;

SLInterfaceCalculator::SLInterfaceCalculator()
{
	pPtSlns = NULL;
	tPtSlns = NULL;
	interfacePFs = NULL;
	sigmaCFactor = 1.0;
	deltaCFactor = 1.0;
	dcont = 0.0;
//	x = 0.0;
	itern = 0;
	current_delT = 0.0;
	current_timeIndex = -1;
	b_directlyCalculateBondedInFinalSln = false;
	ts_bulkProps = NULL;
	interfacePFs = NULL;
	tPtSlns_Deletable = true;
	pPtSlns_Deletable = true;
	ts_bulkProps_Deletable = true;
	interfacePFs_Deletable = true;
	ptTimeSequenceSlns_Deletable = false;
	oneIntAlltimesPtr = NULL;

	terminateConditionMet = false;
	ptTimeSequenceSlns = NULL;
	ortiz1DTSRptr = NULL;
	incidentSide = ilt_noSided;
}

SLInterfaceCalculator::~SLInterfaceCalculator()
{
	if ((tPtSlns != NULL) && (tPtSlns_Deletable == true))
		delete tPtSlns;
	if ((pPtSlns != NULL) && (pPtSlns_Deletable == true))
		delete pPtSlns;
	if ((ts_bulkProps != NULL) && (ts_bulkProps_Deletable == true))
		delete ts_bulkProps;
	if ((interfacePFs != NULL) && (interfacePFs_Deletable == true))
		delete interfacePFs;
	if (ptTimeSequenceSlns_Deletable)
		delete ptTimeSequenceSlns;
	if (ortiz1DTSRptr != NULL)
		delete 	ortiz1DTSRptr;
}

SLInterfaceCalculator::SLInterfaceCalculator(const SLInterfaceCalculator& other)
{
	pPtSlns = NULL;
	tPtSlns = NULL;
	interfacePFs = NULL;
	ts_bulkProps = NULL;
	interfacePFs = NULL;
	tPtSlns_Deletable = true;
	pPtSlns_Deletable = true;
	ts_bulkProps_Deletable = true;
	interfacePFs_Deletable = true;
	ptTimeSequenceSlns_Deletable = false;
	oneIntAlltimesPtr = NULL;
	ptTimeSequenceSlns = NULL;
	(*this) = other;
}

SLInterfaceCalculator& SLInterfaceCalculator::operator=(const SLInterfaceCalculator& other)
{
	itern = other.itern;
	b_directlyCalculateBondedInFinalSln = other.b_directlyCalculateBondedInFinalSln;
	sigmaCFactor = other.sigmaCFactor;
	deltaCFactor = other.deltaCFactor;
	dcont = other.dcont;
	current_delT = other.current_delT;
	current_timeIndex = other.current_timeIndex;
	within_step_convergence_check = other.within_step_convergence_check;
	between_steps_timestep_check = other.between_steps_timestep_check;
	between_steps_adaptivity_s = other.between_steps_adaptivity_s;

	terminateConditionMet = other.terminateConditionMet;

	CopyVec(other.wlSide_rGoing, wlSide_rGoing);
	CopyVec(other.wrSide_lGoing, wrSide_lGoing);

	CopyVec(other.wlSide_rGoing_WO_in_situ, wlSide_rGoing_WO_in_situ);
	CopyVec(other.wrSide_lGoing_WO_in_situ, wrSide_lGoing_WO_in_situ);

	////
	if (tPtSlns != NULL)
	{
		if (tPtSlns_Deletable)
			delete tPtSlns;
	}
	tPtSlns = NULL;
	tPtSlns_Deletable = other.tPtSlns_Deletable;
	if (other.tPtSlns != NULL)
	{
		if (tPtSlns_Deletable == true)
		{
			tPtSlns = new SL_interface_Temp_PPtData();
			*tPtSlns = *other.tPtSlns;
		}
		else
			tPtSlns = other.tPtSlns;
	}

	////
	if (pPtSlns != NULL)
	{
		if (pPtSlns_Deletable)
			delete pPtSlns;
	}
	pPtSlns = NULL;
	pPtSlns_Deletable = other.pPtSlns_Deletable;
	if (other.pPtSlns != NULL)
	{
		if (pPtSlns_Deletable == true)
		{
			pPtSlns = new SL_interfacePPtData();
			*pPtSlns = *other.pPtSlns;
		}
		else
			pPtSlns = other.pPtSlns;
	}

	////
	if (ts_bulkProps != NULL)
	{
		if (ts_bulkProps_Deletable)
			delete ts_bulkProps;
	}
	ts_bulkProps = NULL;
	ts_bulkProps_Deletable = other.ts_bulkProps_Deletable;
	if (other.ts_bulkProps != NULL)
	{
		if (ts_bulkProps_Deletable == true)
		{
			ts_bulkProps = new SL_Elastic_InterfaceProperties();
			*ts_bulkProps = *other.ts_bulkProps;
		}
		else
			ts_bulkProps = other.ts_bulkProps;
	}

	////
	if (interfacePFs != NULL)
	{
		if (interfacePFs_Deletable)
			delete interfacePFs;
	}
	interfacePFs = NULL;
	interfacePFs_Deletable = other.interfacePFs_Deletable;
	if (other.interfacePFs != NULL)
	{
		if (interfacePFs_Deletable == true)
		{
			interfacePFs = new SL_Interface_Fracture_PF();
			*interfacePFs = *other.interfacePFs;
		}
		else
			interfacePFs = other.interfacePFs;
	}

	////
	if (ptTimeSequenceSlns != NULL)
	{
		if (ptTimeSequenceSlns_Deletable)
			delete ptTimeSequenceSlns;
	}
	ptTimeSequenceSlns = NULL;
	ptTimeSequenceSlns_Deletable = other.ptTimeSequenceSlns_Deletable;
	if (other.ptTimeSequenceSlns != NULL)
	{
		if (ptTimeSequenceSlns_Deletable == true)
		{
			ptTimeSequenceSlns = new SL_interfacePPtData_Time_Seq();
			*ptTimeSequenceSlns = *other.ptTimeSequenceSlns;
		}
		else
			ptTimeSequenceSlns = other.ptTimeSequenceSlns;
	}

	oneIntAlltimesPtr = other.oneIntAlltimesPtr;
	return *this;
}

void SLInterfaceCalculator::Compute_Riemann_I_Bonded()
{
	if (!b_directlyCalculateBondedInFinalSln)
	{
		VEC* vIPtr = &tPtSlns->sl_side_temp_ptData[SDL].v_mode_Star[rmode_stick];
		VEC* sIPtr = &tPtSlns->sl_side_temp_ptData[SDL].sigma_mode_Star[rmode_stick];
		ts_bulkProps->Compute_sigmaI_vI_from_ws(wlSide_rGoing, wrSide_lGoing, *vIPtr, *sIPtr);
		CopyVec(*vIPtr, tPtSlns->sl_side_temp_ptData[SDR].v_mode_Star[rmode_stick]);
		CopyVec(*sIPtr, tPtSlns->sl_side_temp_ptData[SDR].sigma_mode_Star[rmode_stick]);
	}
	else // skip temp-data altogether
	{
		VEC* vIPtr = &pPtSlns->sl_side_ptData[SDL].v_downstream_final;
		VEC* sIPtr = &pPtSlns->sl_side_ptData[SDL].sigma_downstream_final;

		if (ts_bulkProps->interfaceLoc == ilt_twoSided)
			ts_bulkProps->Compute_sigmaI_vI_from_ws(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, *vIPtr, *sIPtr);
		else
			ts_bulkProps->Compute_sigmaI_vI_from_wsAndBC__BoundaryCase(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, *vIPtr, *sIPtr);
		CopyVec(*vIPtr, pPtSlns->sl_side_ptData[SDR].v_downstream_final);
		CopyVec(*sIPtr, pPtSlns->sl_side_ptData[SDR].sigma_downstream_final);
	}
}

SLFF_ContactType SLInterfaceCalculator::Compute_Riemann_III_Separation(SLFF_InterfacialDamageModeType damageOffOnMix, double& relCont)
{
	double DDot, DTarget, DSrc, DSrc_Allix, DSrc_RBH, effS;
	VEC sigmaIII;
	// computing commong traction
	SLFF_ContactType contactOffOnMix = interfacePFs->contactOffOnMix;
	Compute_Riemann_III_Separation_Damage_Rate_Aux(damageOffOnMix, contactOffOnMix, tPtSlns->interfacePropPtr->sigma_I_nt_parts, tPtSlns->interfacePropPtr->del_u_nt_parts,
		tPtSlns->sl_side_temp_ptData[SDL].v_downstream_latestValue,
		tPtSlns->sl_side_temp_ptData[SDR].v_downstream_latestValue,
		pPtSlns->maxEffDelU, pPtSlns->tsrStage, tPtSlns->tsrStageLatestValue, tPtSlns->interface_damageLatestValue,
		tPtSlns->sigmaC, tPtSlns->deltaC, tPtSlns->tau, interfacePFs,
		tPtSlns->dn, tPtSlns->dt,
		tPtSlns->interfacePropPtr->del_v_nt_parts, 
		DDot, effS, DTarget, DSrc_Allix, DSrc_RBH, sigmaIII, terminateConditionMet);
	// computing velocities for two sides
	if (contactOffOnMix == sl_contact_on)
	{
		relCont = 1.0;
		return contactOffOnMix;
	}

	SetStressComputeVelocityStarFromStressStar(sigmaIII, rmode_sep);
	// velocities
	VEC *vIIIL = &tPtSlns->sl_side_temp_ptData[SDL].v_mode_Star[rmode_sep],
		*vIIIR = &tPtSlns->sl_side_temp_ptData[SDR].v_mode_Star[rmode_sep];
	double v_IIIL0, v_IIIR0, delvn_sepVal, dpn_separation, dps_separation;
	v_IIIL0 = (*vIIIL)[0];
	v_IIIR0 = (*vIIIR)[0];
	delvn_sepVal = v_IIIR0 - v_IIIL0;

	dpn_separation = delvn_sepVal * sigmaIII[0];
	dps_separation = 0.0;
#if DiM2a3_F
	dps_separation = ((*vIIIR)[1] - (*vIIIL)[1]) * sigmaIII[1];
#if DiM3
	dps_separation += ((*vIIIR)[2] - (*vIIIL)[2]) * sigmaIII[2];
#endif
#endif
	DSrc = DSrc_Allix + DSrc_RBH;

	// setting scalar values
	SL_Interface_Temp_PtData_IntContFrac* interfacePropPtr = tPtSlns->interfacePropPtr;
	interfacePropPtr->interface_scalar_Vals[s1i_effectiveStrs] = effS;
	interfacePropPtr->interface_scalar_Vals[s1i_Dtarget] = DTarget;
	interfacePropPtr->interface_scalar_Vals[s1i_DsrcAllix] = DSrc_Allix;
	interfacePropPtr->interface_scalar_Vals[s1i_DsrcRBHV] = DSrc_RBH;
	interfacePropPtr->interface_scalar_Vals[s1i_Dsrc] = DSrc;
	interfacePropPtr->interface_scalar_Vals[s1i_dpn_separation] = dpn_separation;
	interfacePropPtr->interface_scalar_Vals[s1i_dps_separation] = dps_separation;

	if (contactOffOnMix == sl_contact_off)
	{
		relCont = 0.0;
		return contactOffOnMix;
	}
	double over_bar_delvnVal = 0.0;

	if (delvn_sepVal > over_bar_delvnVal)
	{
		contactOffOnMix = sl_contact_off;
		relCont = 0.0;
		return contactOffOnMix;
	}
	// now velocity is driving the two sides together and we need to look at the velocity value
	double delun_Val = tPtSlns->interfacePropPtr->del_u_nt_parts.vec_n;
	if (delun_Val <= dcont)
	{
		contactOffOnMix = sl_contact_on;
		relCont = 1.0;
		return contactOffOnMix;
	}
	if (delun_Val >= tPtSlns->dsep)
	{
		contactOffOnMix = sl_contact_off;
		relCont = 0.0;
		return contactOffOnMix;
	}
	// it already is in this state but just to be sure
	SLFF_ContactSeparationRegularizationType sep2cont_regType = interfacePFs->sep2cont_regType;
	contactOffOnMix = sl_contact_mixed;
	double factor = 1.0 / (tPtSlns->dsep - dcont);
	double cRaw = (tPtSlns->dsep - delun_Val) * factor;
	double cComputed, tmpVal;
	if (sep2cont_regType == csreg_sqrt1)
	{
		tmpVal = sqrt(1.0 - cRaw);
		cComputed = 1 - tmpVal;
	}
	else if (sep2cont_regType == csreg_linear)
		cComputed = cRaw;
	else if (sep2cont_regType == csreg_sqrt2)
	{
		tmpVal = sqrt(1.0 - cRaw * cRaw);
		cComputed = 1 - tmpVal;
	}
	relCont = cComputed;
	return contactOffOnMix;
}

onOffT SLInterfaceCalculator::Compute_Riemann_II_Slip(double& relStick, SLFF_SlipType& slipOffOnMix)
{
#if DiM1
	relStick = 1.0;
	return allOffV;
#endif
	slipOffOnMix = interfacePFs->slipOffOnMix;
	if (interfacePFs->slipOffOnMix == sl_slip_off)
	{
		relStick = 1.0;
		return allOffV;
	}
	double k, thetaTauFriction, delvIIAbs, thetaDelv, powerIDissipationSlip;
	VEC sigmaII;
	onOffT retVal = Compute_Riemann_II_Slip_Aux(interfacePFs->slipOffOnMix,
		tPtSlns->interfacePropPtr->sigma_I_nt_parts, tPtSlns->interfacePropPtr->del_u_nt_parts,
		interfacePFs, ts_bulkProps, 
		sigmaII, k, thetaTauFriction, delvIIAbs, thetaDelv, powerIDissipationSlip);

	if (retVal == allOnV)
	{
		relStick = 0.0;
		slipOffOnMix = sl_slip_on;
		// computing velocities for two sides
		SetStressComputeVelocityStarFromStressStar(sigmaII, rmode_slip);
	}
	else if (retVal == allOffV)
	{
		relStick = 1.0;
		slipOffOnMix = sl_slip_off;
		// solution is equal to stick
		for (int sidei = 0; sidei < NUM_SIDES; ++sidei)
		{
			CopyVec(tPtSlns->sl_side_temp_ptData[sidei].sigma_mode_Star[rmode_stick], tPtSlns->sl_side_temp_ptData[sidei].sigma_mode_Star[rmode_slip]);
			CopyVec(tPtSlns->sl_side_temp_ptData[sidei].v_mode_Star[rmode_stick], tPtSlns->sl_side_temp_ptData[sidei].v_mode_Star[rmode_slip]);
		}
	}
	else
		THROW("Cannot have partial stick/slip\n");

	// setting scalar values
	SL_Interface_Temp_PtData_IntContFrac* interfacePropPtr = tPtSlns->interfacePropPtr;
	interfacePropPtr->interface_scalar_Vals[s1i_fricCoef] = k; // interfacePFs->Get_Const_friction_Coef();
	interfacePropPtr->interface_scalar_Vals[s1i_slip_dvelTheta] = thetaDelv;
	interfacePropPtr->interface_scalar_Vals[s1i_slip_tauTheta] = thetaTauFriction;
	interfacePropPtr->interface_scalar_Vals[s1i_tau_minus_delv_Theta] = thetaTauFriction - thetaDelv;
	interfacePropPtr->interface_scalar_Vals[s1i_magnitude_delv] = delvIIAbs;
	interfacePropPtr->interface_scalar_Vals[s1i_dps_slip] = powerIDissipationSlip;
	return retVal;
}

void SLInterfaceCalculator::Read(istream& in)
{
	setValue(wlSide_rGoing_WO_in_situ, 0.0);
	setValue(wrSide_lGoing_WO_in_situ, 0.0);
	VEC ul, ur, vl, vr, sl, sr;
	bool u_read = false, v_read = false, s_read = false;

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
		if (buf == "itern")
		{
			READ_NINTEGER(in, buf, itern);
		}
		else if (buf == "u")
		{
			ReadV(ul, in);
			ReadV(ur, in);
			u_read = true;
		}
		else if (buf == "v")
		{
			ReadV(vl, in);
			ReadV(vr, in);
			v_read = true;
		}
		else if (buf == "s")
		{
			ReadV(sl, in);
			ReadV(sr, in);
			s_read = true;
		}
		else if (buf == "wlSide_rGoing_WO_in_situ")
		{
			ReadV(wlSide_rGoing_WO_in_situ, in);
		}
		else if (buf == "wrSide_lGoing_WO_in_situ")
		{
			ReadV(wrSide_lGoing_WO_in_situ, in);
		}
		else if (buf == "ts_bulkProps")
		{
			if ((ts_bulkProps != NULL) && (ts_bulkProps_Deletable == true))
				delete ts_bulkProps;
			ts_bulkProps = new SL_Elastic_InterfaceProperties();
			ts_bulkProps_Deletable = true;
			ts_bulkProps->Read(in);
		}
		else if (buf == "interfacePFs")
		{
			if ((interfacePFs != NULL) && (interfacePFs_Deletable == true))
				delete interfacePFs;
			interfacePFs = new SL_Interface_Fracture_PF();
			interfacePFs_Deletable = true;
			interfacePFs->Read_SL_Interface_Fracture_PF(in, 1);
		}
		else if (buf == "tPtSlns")
		{
			if ((tPtSlns != NULL) && (tPtSlns_Deletable == true))
				delete tPtSlns;
			Create_TemporaryStorage();
			tPtSlns->Read(in);
		}
		else if (buf == "ptTimeSequenceSlns")
		{
			if (ptTimeSequenceSlns_Deletable == true)
				delete ptTimeSequenceSlns;
			ptTimeSequenceSlns = new SL_interfacePPtData_Time_Seq();
			ptTimeSequenceSlns_Deletable = true;
			ptTimeSequenceSlns->Read(in);
		}
		else if (buf == "sigmaCFactor")
		{
			READ_NDOUBLE(in, buf, sigmaCFactor);
		}
		else if (buf == "deltaCFactor")
		{
			READ_NDOUBLE(in, buf, deltaCFactor);
		}
		else if (buf == "dcont")
		{
			READ_NDOUBLE(in, buf, dcont);
		}
		// other things are computed as a part of computation (same for tPtSlns but still may contain most recent displacements ...)
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	if (u_read || v_read || s_read)
	{
		VEC *ulPtr = NULL, *urPtr = NULL, *vlPtr = NULL, *vrPtr = NULL, *slPtr = NULL, *srPtr = NULL;
		if (u_read)
		{
			ulPtr = &ul;
			urPtr = &ur;
		}
		if (v_read)
		{
			vlPtr = &vl;
			vrPtr = &vr;
		}
		if (s_read)
		{
			slPtr = &sl;
			srPtr = &sr;
		}
		Initialize_setFromOutside_Traces_vel_stress_etc(slPtr, srPtr, vlPtr, vrPtr, ulPtr, urPtr);
	}
}

void SLInterfaceCalculator::Create_TemporaryStorage()
{
	tPtSlns = new SL_interface_Temp_PPtData();
	tPtSlns_Deletable = true;
	tPtSlns->sigmaC = sigmaCFactor * interfacePFs->gen_strength;
	tPtSlns->deltaC = deltaCFactor * interfacePFs->deltaC;
	tPtSlns->tau = deltaCFactor * interfacePFs->tau;
	if (interfacePFs->TRS_Tensile_Mode())
	{
		tPtSlns->dn = tPtSlns->deltaC;
		tPtSlns->dt = tPtSlns->dn / interfacePFs->beta_delU;
	}
	else
	{
		tPtSlns->dn = tPtSlns->tau * tPtSlns->sigmaC * ts_bulkProps->YAve_n;
		tPtSlns->dt = tPtSlns->tau * tPtSlns->sigmaC * ts_bulkProps->YAve_t;
	}
	tPtSlns->dsep = tPtSlns->dn * interfacePFs->rel_dsep;
	dcont = tPtSlns->dn * interfacePFs->rel_dcont;
}

void SLInterfaceCalculator::Initialize_setFromOutside_ElsticFractureProperties(SL_Elastic_InterfaceProperties*	ts_bulkPropsIn, SL_Interface_Fracture_PF* interfacePFsIn)
{
	if (ts_bulkPropsIn != NULL)
	{
		if ((ts_bulkProps != NULL) && (ts_bulkProps_Deletable))
			delete ts_bulkProps;
		ts_bulkProps = ts_bulkPropsIn;
		ts_bulkProps_Deletable = false;
	}

	if (interfacePFsIn != NULL)
	{
		if ((interfacePFs != NULL) && (interfacePFs_Deletable))
			delete interfacePFs;
		interfacePFs = interfacePFsIn;
		interfacePFs_Deletable = false;
	}
}

void SLInterfaceCalculator::Initialize_setFromOutside_Current_Step_EmptyPermanentPt_Storage(SL_interfacePPtData* pPtSlnsIn) //, double xIn)
{
//	x = xIn;
	if (pPtSlnsIn != NULL)
	{
		if ((pPtSlns != NULL) && (pPtSlns_Deletable))
			delete pPtSlns;
		pPtSlns = pPtSlnsIn;
		pPtSlns_Deletable = false;
	}
}

void SLInterfaceCalculator::Initialize_setFromOutside_Earlier_Steps_NonEmptyPermanentPts_Storage(SL_interfacePPtData_Time_Seq& ptTimeSequenceSlnsIn)
{
	ptTimeSequenceSlns = &ptTimeSequenceSlnsIn;
	ptTimeSequenceSlns_Deletable = false;
}

void SLInterfaceCalculator::Initialize_setFromOutside_Incoming_Characteristics_etc(const VEC* wlSide_rGoing_WO_in_situPtrIn, const VEC* wrSide_lGoing_WO_in_situPtrIn, const VEC* traceCurrentDisplacementLeftIn, const VEC* traceCurrentDisplacementRightIn, const VEC* traceCurrentVelLeftIn, const VEC* traceCurrentVelRightIn)
{
	if (wlSide_rGoing_WO_in_situPtrIn != NULL)
		CopyVec(*wlSide_rGoing_WO_in_situPtrIn, wlSide_rGoing_WO_in_situ);
	if (wrSide_lGoing_WO_in_situPtrIn != NULL)
		CopyVec(*wrSide_lGoing_WO_in_situPtrIn, wrSide_lGoing_WO_in_situ);
	bool provided_side_disp_traces = (traceCurrentDisplacementLeftIn != NULL);
	if (!provided_side_disp_traces)
		return;
	bool b_directlyCalculateBondedInFinalSln = DirectlyCalculateBondedInFinalSln();
	// no side traces of displacement are needed for non-fracture / non-BC cases
	if (b_directlyCalculateBondedInFinalSln)
		return;
	if ((ptTimeSequenceSlns != NULL) && (!ptTimeSequenceSlns->IsEmpty()))
	{
		THROW("traces of side displacements should only be provided when there is no earlier solutions. Providing traces of displacement is ONLY needed when there is no earlier solution (common for 1 point computations).\n");
	}
	if (tPtSlns == NULL)
		Create_TemporaryStorage();

	// now can put the values in temp-storage
	CopyVec(*traceCurrentDisplacementLeftIn, tPtSlns->sl_side_temp_ptData[SDL].u_downstream_latestValue);
	if (traceCurrentDisplacementRightIn != NULL)
		CopyVec(*traceCurrentDisplacementRightIn, tPtSlns->sl_side_temp_ptData[SDR].u_downstream_latestValue);
	bool compute_ss = false;
	if (traceCurrentVelLeftIn != NULL)
	{
		CopyVec(*traceCurrentVelLeftIn, tPtSlns->sl_side_temp_ptData[SDL].v_downstream_latestValue);
		compute_ss = true;
	}
	if (traceCurrentVelRightIn != NULL)
	{
		CopyVec(*traceCurrentVelRightIn, tPtSlns->sl_side_temp_ptData[SDR].v_downstream_latestValue);
		compute_ss = true;
	}
	if (compute_ss)
		ts_bulkProps->Compute_stresses_from_vs_and_characteristicss(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, traceCurrentVelLeftIn, traceCurrentVelRightIn,
			tPtSlns->sl_side_temp_ptData[SDL].sigma_downstream_latestValue, tPtSlns->sl_side_temp_ptData[SDR].sigma_downstream_latestValue);
}

void SLInterfaceCalculator::Initialize_setFromOutside_Traces_vel_stress_etc(const VEC* traceCurrentStress_WO_in_situ_LeftIn, const VEC* traceCurrentStress_WO_in_situ_RightIn, const VEC* traceCurrentVelLeftIn, const VEC* traceCurrentVelRightIn, VEC* traceCurrentDisplacementLeftIn, VEC* traceCurrentDisplacementRightIn)
{
	ts_bulkProps->Computecharacteristics_from_Traces_vel_stress_etc(traceCurrentStress_WO_in_situ_LeftIn, traceCurrentStress_WO_in_situ_RightIn, traceCurrentVelLeftIn, traceCurrentVelRightIn, wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ);
	bool provided_side_disp_traces = (traceCurrentDisplacementLeftIn != NULL);
	if (!provided_side_disp_traces)
		return;
	bool b_directlyCalculateBondedInFinalSln = DirectlyCalculateBondedInFinalSln();
	// no side traces of displacement are needed for non-fracture / non-BC cases
	if (b_directlyCalculateBondedInFinalSln)
		return;
	if ((ptTimeSequenceSlns != NULL) && (!ptTimeSequenceSlns->IsEmpty()))
	{
		THROW("traces of side displacements should only be provided when there is no earlier solutions. Providing traces of displacement is ONLY needed when there is no earlier solution (common for 1 point computations).\n");
	}
	if (tPtSlns == NULL)
		Create_TemporaryStorage();
	// now can put the values in temp-storage
	CopyVec(*traceCurrentDisplacementLeftIn, tPtSlns->sl_side_temp_ptData[SDL].u_downstream_latestValue);
	if (traceCurrentDisplacementRightIn != NULL)
		CopyVec(*traceCurrentDisplacementRightIn, tPtSlns->sl_side_temp_ptData[SDR].u_downstream_latestValue);
}

void SLInterfaceCalculator::Initialize_setFromOutside_Fracture_InhomogeneousFactors(double sigmaCFactorIn, double deltaCFactorIn)
{
	sigmaCFactor = sigmaCFactorIn;
	deltaCFactor = deltaCFactorIn;
}

AdaptivityS SLInterfaceCalculator::Main_Compute_OnePoint(bool& accept_point, IOF_type iof_finalVals, IOF_type iof_scalarVals, double space_or_time, ostream* outScalars, ostream* outFinalSln, ostream* outAdaptivity, ostream* outIterationConv)
{
	if (pPtSlns == NULL)
	{
		pPtSlns = new SL_interfacePPtData();
		pPtSlns_Deletable = true;
	}
	// setting defau
	between_steps_adaptivity_s.a_delt = g_slf_conf->uniform_del_t;
	between_steps_adaptivity_s.Update_AdaptFlag(g_slf_conf->a_flag_default);

	b_directlyCalculateBondedInFinalSln = DirectlyCalculateBondedInFinalSln();
	if (!b_directlyCalculateBondedInFinalSln)
	{
		if (tPtSlns == NULL)
			Create_TemporaryStorage();
		bool iterations_enough = false;
		while (!iterations_enough)
		{
			if (interfacePFs->IsSimpleOrtizModel())
			{
				if (ortiz1DTSRptr != NULL)
					delete ortiz1DTSRptr;
				ortiz1DTSRptr = new TSR1d_stepn_to_np1();
			}
			// form "most recent" traces of v, s from earlier values (iteration 0) and recently computed target values (iterations > 0)
			InitializeOrUpdate_Most_Current_Values();
			DBCHK(dbout << "\nitern" << itern << '\n';);

			bool checkConvergence = (ortiz1DTSRptr == NULL);
#if HAVE_SOURCE
			checkConvergence = true;
			// updating characteristics
			if (oneIntAlltimesPtr != NULL)
				oneIntAlltimesPtr->Compute_DownStream_Characteristics_wo_in_situ(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, pPtSlns->interface_time, current_timeIndex);
#endif			
			// compute star values and update the stats
			Compute_1Cycle_vsStar_Calculation();


			if (checkConvergence)
				iterations_enough = Check_withinstep_iteration_convergence();
			else
				iterations_enough = true;
			DBCHK(dbout << "iterations_enough" << iterations_enough << '\n';);
		}
	}
	else
		Compute_1Cycle_vsStar_Calculation();

	// at this point iterations are completed and can compute the errors
	Copy_temp_to_final_vals_Check_adaptivityFlag(accept_point);

	//checking if final run conditions are met ...
	if ((terminateConditionMet) || (interfacePFs->damageOffOnMix == sl_interfacial_damage_mixed) && (pPtSlns->interface_damage_final >= g_slf_conf->terminate_run_target_max_damage))
		between_steps_adaptivity_s.Update_AdaptFlag(a_terminate_run_correctly);

	if (outFinalSln != NULL)
		Output_SLInterfaceCalculator(accept_point, iof_finalVals, iof_scalarVals, space_or_time, outScalars, outFinalSln, outAdaptivity, outIterationConv);

	return between_steps_adaptivity_s;
}

void SLInterfaceCalculator::Output_SLInterfaceCalculator(bool accept_point, IOF_type iof_finalVals, IOF_type iof_scalarVals, double space_or_time, ostream* outScalars, ostream* outFinalSln, ostream* outAdaptivity, ostream* outIterationConv)
{
	if (accept_point)
	{
		double x = 0.0, t = pPtSlns->interface_time;
		bool has_ring_opened1D_al = false;
		if (oneIntAlltimesPtr != NULL)
		{
			x = oneIntAlltimesPtr->interface_x;
			has_ring_opened1D_al = oneIntAlltimesPtr->has_ring_opened1D_al;
		}

		if ((tPtSlns != NULL) && !b_directlyCalculateBondedInFinalSln && (outScalars != NULL))
			tPtSlns->Output_ScalarValues(*outScalars, iof_scalarVals, space_or_time);
		if (outFinalSln != NULL)
			pPtSlns->Output_FinalSolution(*outFinalSln, iof_finalVals, space_or_time, x, t, has_ring_opened1D_al);
	}
	if ((outAdaptivity != NULL) && (g_slf_conf->between_steps_adaptivity))
	{
		*outAdaptivity << between_steps_adaptivity_s << '\t';
		between_steps_timestep_check.PrintShort(*outAdaptivity, pPtSlns->interface_time, itern);
		*outAdaptivity << '\n';
	}
	if ((outIterationConv != NULL) && !b_directlyCalculateBondedInFinalSln && (interfacePFs->Need_BetweenIteration_ConvergenceCheck()))
	{
		within_step_convergence_check.PrintShort(*outIterationConv, pPtSlns->interface_time, itern);
		*outIterationConv << '\n';
	}
}

void SLInterfaceCalculator::InitializeOrUpdate_Most_Current_Values()
{
	bool iter0 = (itern == 0);
	bool hasEarlierVals = (ptTimeSequenceSlns != NULL);
	SL_interfacePPtData* step_n_PtSln = NULL;
	int szEarlierVals = 0;
	if (hasEarlierVals)
	{
		szEarlierVals = ptTimeSequenceSlns->Earlier_size();
		if (szEarlierVals <= 0)
			hasEarlierVals = false;
		else
			step_n_PtSln = ptTimeSequenceSlns->GetBackwardPosition(0);
	}

#if RING_PROBLEM
	double x = 0;
	if (oneIntAlltimesPtr != NULL)
		x = oneIntAlltimesPtr->interface_x;
	double factor_vr_source_sigma_theta = 0.5 * (ts_bulkProps->bulk_leftPtr->E_iso + ts_bulkProps->bulk_rightPtr->E_iso) / g_domain->ring_R;
	double ringp_n = g_SL_desc_data.GetRing_p(x, step_n_PtSln->interface_time), ringp_np1 = g_SL_desc_data.GetRing_p(x, pPtSlns->interface_time);
	double factor_sigmatheta_source_vr = -0.5 * (1.0 / ts_bulkProps->bulk_leftPtr->rho + 1.0 / ts_bulkProps->bulk_rightPtr->rho) / g_domain->ring_R;
#endif

	// taking care of values that are not obtained from ODEs
	if (!iter0)
	{
		CopyVec(tPtSlns->sl_side_temp_ptData[SDL].sigma_Star, tPtSlns->sl_side_temp_ptData[SDL].sigma_downstream_latestValue);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDR].sigma_Star, tPtSlns->sl_side_temp_ptData[SDR].sigma_downstream_latestValue);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDL].v_Star, tPtSlns->sl_side_temp_ptData[SDL].v_downstream_latestValue);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDR].v_Star, tPtSlns->sl_side_temp_ptData[SDR].v_downstream_latestValue);
	}
	else if (hasEarlierVals)
	{
//		if (ortiz1DTSRptr == NULL)
		{
			CopyVec(step_n_PtSln->sl_side_ptData[SDL].sigma_downstream_final, tPtSlns->sl_side_temp_ptData[SDL].sigma_downstream_latestValue);
			CopyVec(step_n_PtSln->sl_side_ptData[SDR].sigma_downstream_final, tPtSlns->sl_side_temp_ptData[SDR].sigma_downstream_latestValue);
			CopyVec(step_n_PtSln->sl_side_ptData[SDL].v_downstream_final, tPtSlns->sl_side_temp_ptData[SDL].v_downstream_latestValue);
			CopyVec(step_n_PtSln->sl_side_ptData[SDR].v_downstream_final, tPtSlns->sl_side_temp_ptData[SDR].v_downstream_latestValue);
		}
	}
#if 0	// intentionally commented out, so if the user has provided nonzero values for these, they are not overwritten
	else
	{
		setValue(tPtSlns->sl_side_temp_ptData[SDL].sigma_downstream_latestValue, 0.0);
		setValue(tPtSlns->sl_side_temp_ptData[SDR].sigma_downstream_latestValue, 0.0);
		setValue(tPtSlns->sl_side_temp_ptData[SDL].v_downstream_latestValue, 0.0);
		setValue(tPtSlns->sl_side_temp_ptData[SDR].v_downstream_latestValue, 0.0);
	}
#endif
	//// now taking care of ODE values
	// dealing with fDot = Sf (source of f) -> f_n+1 = f_n + factor_n * Sf_n + factor_n+1 * Sf_n+1
	if (tPtSlns != NULL)
	{
		if (interfacePFs->damageOffOnMix == sl_interfacial_damage_off)
			tPtSlns->interface_damageLatestValue = 0.0;
		else if (interfacePFs->damageOffOnMix == sl_interfacial_damage_on)
			tPtSlns->interface_damageLatestValue = 1.0;
	}
	if (hasEarlierVals)
	{
		double delT = pPtSlns->interface_time - step_n_PtSln->interface_time;
		current_delT = delT;
		// factor of source
		double factor_Step_n, factor_Step_np1 = 0.0;
		bool has_factor_Step_np1 = !iter0;
		if (has_factor_Step_np1)
		{
			factor_Step_n = g_slf_conf->one_step_update_alpha_it_gt0 * delT;
			factor_Step_np1 = delT - factor_Step_n;
		}
		else
			factor_Step_n = g_slf_conf->one_step_update_alpha_it0 * delT;

		if (ortiz1DTSRptr != NULL)
		{
			ortiz1DTSRptr->step_n.sides[SDL].u = step_n_PtSln->sl_side_ptData[SDL].u_downstream_final[0];
			ortiz1DTSRptr->step_n.sides[SDR].u = step_n_PtSln->sl_side_ptData[SDR].u_downstream_final[0];

			ortiz1DTSRptr->step_n.sides[SDL].v = step_n_PtSln->sl_side_ptData[SDL].v_downstream_final[0];
			ortiz1DTSRptr->step_n.sides[SDR].v = step_n_PtSln->sl_side_ptData[SDR].v_downstream_final[0];
			if ((oneIntAlltimesPtr != NULL) && (oneIntAlltimesPtr->has_ring_opened1D_al))
			{
				ortiz1DTSRptr->step_n.sides[SDR].v += g_domain->ring_opened1D_al;
				ortiz1DTSRptr->step_n.sides[SDR].u += step_n_PtSln->interface_time * g_domain->ring_opened1D_al;
			}
			ortiz1DTSRptr->step_n_max_delu = step_n_PtSln->maxEffDelU;
			ortiz1DTSRptr->step_n_tsrStage = step_n_PtSln->tsrStage;
			ortiz1DTSRptr->delT = delT;
		}
		else
		{

			// now computing values:
			// A. displacement from the two sides
			VEC* us_np1[2];
			us_np1[SDL] = &tPtSlns->sl_side_temp_ptData[SDL].u_downstream_latestValue;
			us_np1[SDR] = &tPtSlns->sl_side_temp_ptData[SDR].u_downstream_latestValue;

			for (int eside = 0; eside < NUM_SIDES; ++eside)
			{
				VEC* uSrc_n = &step_n_PtSln->sl_side_ptData[eside].v_downstream_final;
				VEC* u_n = &step_n_PtSln->sl_side_ptData[eside].u_downstream_final;
				VEC* u_np1 = us_np1[eside];
				for (int j = 0; j < DiM; ++j)
					(*u_np1)[j] = (*u_n)[j] + factor_Step_n * (*uSrc_n)[j];
				if (has_factor_Step_np1)
				{
					VEC* uSrc_np1 = &tPtSlns->sl_side_temp_ptData[eside].v_Star;
					for (int j = 0; j < DiM; ++j)
						(*u_np1)[j] += factor_Step_np1 * (*uSrc_np1)[j];
				}
			}
			// if del_u 0 causes too much penetration, we want to fix this
			double delu0 = (*(us_np1[SDR]))[0] - (*(us_np1[SDL]))[0];
			if ((oneIntAlltimesPtr != NULL) && (oneIntAlltimesPtr->has_ring_opened1D_al))
				delu0 += step_n_PtSln->interface_time * g_domain->ring_opened1D_al;
			double diff = dcont - delu0;
			if (diff > 0.0)
			{
				if ((g_slf_conf->between_steps_adaptivity) && (interfacePFs->contactOffOnMix == sl_contact_mixed))
				{
					double tol = g_slf_conf->between_steps_del_sep2cont_c;
					if (tol > 0.0)
					{
						double rat = diff / (tPtSlns->dsep - dcont);
						if (rat > tol)
						{
							between_steps_adaptivity_s.a_delt = current_delT * g_slf_conf->refinement_delt_factor;
							if (between_steps_adaptivity_s.a_delt > g_slf_conf->min_del_t)
								between_steps_adaptivity_s.Update_AdaptFlag(a_refine);
							else
								between_steps_adaptivity_s.Update_AdaptFlag(a_terminate_run_prematurely);
						}
					}
				}
				diff *= 0.5;
				(*(us_np1[SDR]))[0] += diff;
				(*(us_np1[SDL]))[0] -= diff;
			}

			// B. Damage
			if (interfacePFs->damageOffOnMix == sl_interfacial_damage_mixed)
			{
				double fSrc_n = step_n_PtSln->interface_damage_source_final;
				double f_n = step_n_PtSln->interface_damage_final;
				double* f_np1 = &tPtSlns->interface_damageLatestValue;
				*f_np1 = f_n + factor_Step_n * fSrc_n;
				if (has_factor_Step_np1)
				{
					double fSrc_np1 = tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_Dsrc];
					*f_np1 += factor_Step_np1 * fSrc_np1;
				}
				if (1.0 - *f_np1 <= interfacePFs->zeroTol4DamageEq1)
					*f_np1 = 1.0;
			}
			pPtSlns->maxEffDelU = step_n_PtSln->maxEffDelU;
			tPtSlns->tsrStageLatestValue = step_n_PtSln->tsrStage;
			pPtSlns->tsrStage = step_n_PtSln->tsrStage;
		}
#if RING_PROBLEM
		{
			// evolution of vr
			double fSrc_n = step_n_PtSln->v_r_source_final;
			double f_n = step_n_PtSln->v_r_final;
			double* f_np1 = &tPtSlns->v_r_latestValue;
			*f_np1 = f_n + factor_Step_n * fSrc_n;
			if (has_factor_Step_np1)
			{
				double* fSrc_np1 = &tPtSlns->v_r_source_latestValue;
				double sigmaTheta = tPtSlns->sl_side_temp_ptData[SDL].sigma_Star[0];
				*fSrc_np1 = ringp_np1 + factor_sigmatheta_source_vr * sigmaTheta;
				*f_np1 += factor_Step_np1 * *fSrc_np1;
			}
			double vrLatest = *f_np1;
			tPtSlns->sigma_theta_source_latestValue = factor_vr_source_sigma_theta * vrLatest;
		}
#endif
	}
	else
	{
		double delT = pPtSlns->interface_time;
		current_delT = delT;
		if (interfacePFs->damageOffOnMix == sl_interfacial_damage_on)
			tPtSlns->interface_damageLatestValue = 1.0;

#if RING_PROBLEM
		tPtSlns->v_r_source_latestValue = ringp_np1;
#endif
	}
}

bool SLInterfaceCalculator::Check_withinstep_iteration_convergence()
{
	// incrementing iterations
	++itern;
	// if the run does not require within iteration convergence check, return true.
	// this is true for linear simulations (D = 0)
	if (!interfacePFs->Need_BetweenIteration_ConvergenceCheck())
		return true;
	// checking if the error is small enough
	SL_Interface_Temp_PtData_OneSide* side_L = &tPtSlns->sl_side_temp_ptData[SDL], *side_R = &tPtSlns->sl_side_temp_ptData[SDR];
	bool acceptConvergence = true;
	within_step_convergence_check.SetZero();
	// A. traction
	if (g_slf_conf->del_s_type_within_step != ect_notActive)
	{
		double tol = g_slf_conf->within_step_iter_del_s_tol;
		if (g_slf_conf->del_s_type_within_step == ect_nonDimensional)
			tol *= tPtSlns->sigmaC;
		VEC diff;
		SubtractVec(side_L->sigma_Star, side_L->sigma_downstream_latestValue, diff);
		double normV = Norm2(diff);
		within_step_convergence_check.norm_del_s = normV;
		if (normV > tol)
			acceptConvergence = false;

#if DB_ON
		if (b_db_p)
		{
			dbout << "sigma_Star\t" << side_L->sigma_Star[0] << '\n';
			dbout << "sigma_downstream_latestValue\t" << side_L->sigma_downstream_latestValue[0] << '\n';
			dbout << "normV\t" << normV << '\t';
			dbout << "tol\t" << tol << '\t';
			dbout << "acceptConvergenceS\t" << acceptConvergence << '\n';
		}
#endif
	}
	// B. velocities
	if ((acceptConvergence) && (g_slf_conf->del_v_type_within_step != ect_notActive))
	{
		double tol = g_slf_conf->within_step_iter_del_v_tol;
		if (g_slf_conf->del_v_type_within_step == ect_nonDimensional)
			tol *= tPtSlns->sigmaC * ts_bulkProps->YAve_n;
		VEC diff;
		SubtractVec(side_L->v_Star, side_L->v_downstream_latestValue, diff);
		double normV = Norm2(diff);
		within_step_convergence_check.norm_del_v = normV;

#if DB_ON
		if (b_db_p)
		{
			dbout << "Lv_Star\t" << side_L->v_Star[0] << '\n';
			dbout << "Lv_downstream_latestValue\t" << side_L->v_downstream_latestValue[0] << '\n';
			dbout << "normV\t" << normV << '\t';
			dbout << "tol\t" << tol << '\t';
			dbout << "acceptConvergenceVL\t" << (normV > tol) << '\n';
		}
#endif
		if (normV > tol)
			acceptConvergence = false;
		else
		{
			SubtractVec(side_R->v_Star, side_R->v_downstream_latestValue, diff);
			double normV = Norm2(diff);
			within_step_convergence_check.norm_del_v = MIN(within_step_convergence_check.norm_del_v, normV);
			if (normV > tol)
				acceptConvergence = false;


			within_step_convergence_check.norm_del_v = normV;

#if DB_ON
			if (b_db_p)
			{
				dbout << "Rv_Star\t" << side_R->v_Star[0] << '\n';
				dbout << "Rv_downstream_latestValue\t" << side_R->v_downstream_latestValue[0] << '\n';
				dbout << "normV\t" << normV << '\t';
				dbout << "tol\t" << tol << '\t';
				dbout << "acceptConvergenceVR\t" << (normV > tol) << '\n';
			}
#endif
		}
	}
	// if accepted, no further iterations is required, true is returned
	if (acceptConvergence == true)
		return true;
	// still below maximum allowable iterations, so will return false for still a non-convergent point -> further iteration is needed
	if (itern < g_slf_conf->within_step_iter_max_num_iter)
	{
		DBCHK(dbout << "itern<g_slf_conf->within_step_iter_max_num_iter-return-false\n";);
		return false;
	}
	DBCHK(dbout << "return-true\n";);

	// now we have reached the maximum iterations and no convergence is achieved. Check options.
	// no refinement is needed after reaching max number of iterations, but still the point is considered converged
	if (!g_slf_conf->refine_after_max_num_iter)
		return true;
	// refinement is needed for the point
	// target time step is decided:
	double target_delt = current_delT * g_slf_conf->refine_after_max_num_iter_del_t_factor;
	between_steps_adaptivity_s.a_delt = target_delt;
	if (!g_slf_conf->between_steps_adaptivity)
	{
		between_steps_adaptivity_s.Update_AdaptFlag(a_terminate_run_prematurely);
		return true;
	}
	// check if this time step is too small
	if (target_delt > g_slf_conf->min_del_t)
		between_steps_adaptivity_s.Update_AdaptFlag(a_refine);
	else
		between_steps_adaptivity_s.Update_AdaptFlag(a_terminate_run_prematurely); // terminate run because the error is too small
	return true;
}

void SLInterfaceCalculator::Copy_temp_to_final_vals_Check_adaptivityFlag(bool& accept_point)
{
	if (incidentSide != ilt_noSided)
		UpdateStarValues_FromIncident();
	else if ((oneIntAlltimesPtr != NULL) && (oneIntAlltimesPtr->has_ring_opened1D_al))
		Update_Velocity_Left_Ring_Opened1D();

	bool hasEarlierVals = (ptTimeSequenceSlns != NULL);
	SL_interfacePPtData* lastSlnPt = NULL;
	int szEarlierVals = 0;
	if (hasEarlierVals)
	{
		szEarlierVals = ptTimeSequenceSlns->Earlier_size();
		if (szEarlierVals <= 0)
			hasEarlierVals = false;
		else
			lastSlnPt = ptTimeSequenceSlns->GetBackwardPosition(0);
		double delT = pPtSlns->interface_time - lastSlnPt->interface_time;
		current_delT = delT;
	}
	else
		current_delT = pPtSlns->interface_time;
	SL_interfacePPtData* currentSlnPt = pPtSlns;

	// A. first copy temp to final values
	// copying v, s, u
	if (!b_directlyCalculateBondedInFinalSln)
	{
		CopyVec(tPtSlns->sl_side_temp_ptData[SDL].v_Star, pPtSlns->sl_side_ptData[SDL].v_downstream_final);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDR].v_Star, pPtSlns->sl_side_ptData[SDR].v_downstream_final);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDL].sigma_Star, pPtSlns->sl_side_ptData[SDL].sigma_downstream_final);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDR].sigma_Star, pPtSlns->sl_side_ptData[SDR].sigma_downstream_final);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDL].u_downstream_latestValue, pPtSlns->sl_side_ptData[SDL].u_downstream_final);
		CopyVec(tPtSlns->sl_side_temp_ptData[SDR].u_downstream_latestValue, pPtSlns->sl_side_ptData[SDR].u_downstream_final);

		pPtSlns->interface_damage_source_final = 0.0;
		pPtSlns->interface_damage_final = tPtSlns->interface_damageLatestValue;
		if (tPtSlns->interfacePropPtr != NULL)
		{
			pPtSlns->interface_damage_source_final = tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_Dsrc];
			pPtSlns->maxEffDelU = MAX(pPtSlns->maxEffDelU, tPtSlns->interfacePropPtr->del_u_nt_parts.getEffective());
			pPtSlns->tsrStage = tPtSlns->tsrStageLatestValue;
		}
#if RING_PROBLEM
		pPtSlns->v_r_final = tPtSlns->v_r_latestValue;
		pPtSlns->sigma_theta_source_final = tPtSlns->sigma_theta_source_latestValue;
		pPtSlns->v_r_source_final = tPtSlns->v_r_source_latestValue;
#endif
	}
	else
	{
		/// only need to take care of u and ring problem stuff - velocity and traction are already put in the right place
		if (lastSlnPt != NULL)
		{
			double delT = current_delT;
			double factor_Step_n = g_slf_conf->one_step_update_alpha_it_gt0 * delT;
			double factor_Step_np1 = delT - factor_Step_n;

			// now computing values:
			// A. displacement from the two sides
			for (int eside = 0; eside < NUM_SIDES; ++eside)
			{
				VEC* uSrc_n = &lastSlnPt->sl_side_ptData[eside].v_downstream_final;
				VEC* u_n = &lastSlnPt->sl_side_ptData[eside].u_downstream_final;

				VEC* uSrc_np1 = &currentSlnPt->sl_side_ptData[eside].v_downstream_final;
				VEC* u_np1 = &currentSlnPt->sl_side_ptData[eside].u_downstream_final;
				for (int j = 0; j < DiM; ++j)
					(*u_np1)[j] = (*u_n)[j] + factor_Step_n * (*uSrc_n)[j] + factor_Step_np1 * (*uSrc_np1)[j];
			}
		}
		else
		{
			// now computing values:
			// A. displacement from the two sides
			double delT = pPtSlns->interface_time;
			for (int eside = 0; eside < NUM_SIDES; ++eside)
			{
				VEC* uSrc_np1 = &currentSlnPt->sl_side_ptData[eside].v_downstream_final;
				VEC* u_np1 = &currentSlnPt->sl_side_ptData[eside].u_downstream_final;
				for (int j = 0; j < DiM; ++j)
					(*u_np1)[j] = delT * (*uSrc_np1)[j];
			}
		}
#if RING_PROBLEM
		THROW("Ring problems don't get here as there is no left or right boundary condition as solution should not be computed directly (nonzero source). This may change in later versions when the implicit solution is obtained.\n");
#endif
	}

	// now the values are all transferred and if the run does not need adaptivity we can return true
	accept_point = true;
	if (!g_slf_conf->between_steps_adaptivity)
	{
		between_steps_adaptivity_s.Update_AdaptFlag(a_none);
		return;
	}
	AdaptivityF a_flag = between_steps_adaptivity_s.get_a_flag();
	// it's possible for iterative scheme to fail prematurely, i.e. reaching max iteration count, and requesting run termination or refine based on g_slf_conf parameters
	if ((a_flag == a_refine) || (a_flag == a_terminate_run_prematurely))
	{
		accept_point = false;
		return;
	}

	if (lastSlnPt == NULL)
	{
		between_steps_adaptivity_s.Update_AdaptFlag(a_none);
		return;
	}

	// run needs adaptivity, so check errors
	// this value is used to check all the errors are small enough
//	bool acceptAdaptivity = true;
	double tol, normV, ratio, error_2_tol_max = 0.0;
	if (g_slf_conf->del_s_type != ect_notActive)
	{
		tol = g_slf_conf->between_steps_del_s_tol;
		if (g_slf_conf->del_s_type == ect_nonDimensional)
			tol *= tPtSlns->sigmaC;
		VEC diff;
		SubtractVec(currentSlnPt->sl_side_ptData[SDL].sigma_downstream_final, lastSlnPt->sl_side_ptData[SDL].sigma_downstream_final, diff);
		normV = Norm2(diff);
		between_steps_timestep_check.norm_del_s = normV;
		ratio = normV / tol;
		error_2_tol_max = MAX(error_2_tol_max, ratio);
//		if (ratio >= 1.0)
//			acceptAdaptivity = false;
	}
	// B. velocities
	if (g_slf_conf->del_v_type != ect_notActive)
	{
		tol = g_slf_conf->between_steps_del_v_tol;
		if (g_slf_conf->del_v_type == ect_nonDimensional)
			tol *= tPtSlns->sigmaC * ts_bulkProps->YAve_n;
		VEC diff;
		SubtractVec(currentSlnPt->sl_side_ptData[SDL].v_downstream_final, lastSlnPt->sl_side_ptData[SDL].v_downstream_final, diff);
		normV = Norm2(diff);
		between_steps_timestep_check.norm_del_v = normV;
		ratio = normV / tol;
		error_2_tol_max = MAX(error_2_tol_max, ratio);
//		if (ratio >= 1.0)
//			acceptAdaptivity = false;

		SubtractVec(currentSlnPt->sl_side_ptData[SDR].v_downstream_final, lastSlnPt->sl_side_ptData[SDR].v_downstream_final, diff);
		normV = Norm2(diff);
		between_steps_timestep_check.norm_del_v = MAX(between_steps_timestep_check.norm_del_v, normV);
		ratio = normV / tol;
		error_2_tol_max = MAX(error_2_tol_max, ratio);
//		if (ratio >= 1.0)
//			acceptAdaptivity = false;
	}
	if (!b_directlyCalculateBondedInFinalSln)
	{
		tol = g_slf_conf->between_steps_del_damage_tol;
		if (tol > 0.0)
		{
			normV = fabs(currentSlnPt->interface_damage_final - lastSlnPt->interface_damage_final);
			between_steps_timestep_check.del_damage = normV;
			ratio = normV / tol;
			error_2_tol_max = MAX(error_2_tol_max, ratio);
			//		if (ratio >= 1.0)
			//			acceptAdaptivity = false;
		}
		tol = g_slf_conf->between_steps_damage_source_tau_tol;
		if (tol > 0.0)
		{
//			normV = fabs(currentSlnPt->interface_damage_source_final * tPtSlns->tau);
			normV = fabs(currentSlnPt->interface_damage_source_final * current_delT);
			between_steps_timestep_check.damage_source_tau = normV;
			ratio = normV / tol;
			error_2_tol_max = MAX(error_2_tol_max, ratio);
			//		if (ratio >= 1.0)
			//			acceptAdaptivity = false;
		}
		if (tPtSlns->contactOffOnMix == sl_contact_mixed)
		{
			tol = g_slf_conf->between_steps_del_sep2cont_c;
			if (tol > 0.0)
			{
				double del_uBefore = lastSlnPt->sl_side_ptData[SDR].u_downstream_final[0] - lastSlnPt->sl_side_ptData[SDL].u_downstream_final[0];
				double del_uNow = tPtSlns->interfacePropPtr->del_u_nt_parts.vec_n;
				double delBeforetoNow = del_uBefore - del_uNow;
				if (delBeforetoNow > 0.0) // going to contact
				{
					double approx_delc = delBeforetoNow / (tPtSlns->dsep - dcont);
					if (approx_delc > tol)
					{
						between_steps_timestep_check.del_sep2cont_c = approx_delc;
						ratio = approx_delc / tol;
						error_2_tol_max = MAX(error_2_tol_max, ratio);
					}
				}
			}
		}
	}
	g_slf_conf->Update_AdaptivityFlag_From_current_delT_max_err2tol(current_delT, error_2_tol_max, between_steps_adaptivity_s, accept_point);
	return;
}

void SLInterfaceCalculator::Compute_1Cycle_vsStar_Calculation()
{
	//1. Computing characteristics with in-situ stresse
	if (!b_directlyCalculateBondedInFinalSln)
		ts_bulkProps->Compute_Downstream_characteristic_From_in_situ_Stress(interfacePFs->has_in_situ, interfacePFs->in_situ_stress,
			wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, wlSide_rGoing, wrSide_lGoing);
	
	if (incidentSide != ilt_noSided)
		UpdateCharacteristics_FromIncident();

	if (ortiz1DTSRptr != NULL)
		Compute_1ConvergedImplcit_vsStar_TSR_OrtizSimple_Calculation();
	else
	{
		//2. Compute the star values
		Compute_Star_sv_from_RiemannItoIII_toStar();

		if (b_directlyCalculateBondedInFinalSln)
			return;
		//3. Computing powers
		// these are area fraction averaged power dissipations
		if (interfacePFs->HasAnyContactFriction())
		{
			double ave_powerIDissipation_n, ave_powerIDissipation_s;
			SL_Interface_Temp_PtData_IntContFrac* interfacePropPtr = tPtSlns->interfacePropPtr;
			double aIII = interfacePropPtr->interface_scalar_Vals[s1i_absSep];
			ave_powerIDissipation_n = aIII * interfacePropPtr->interface_scalar_Vals[s1i_dpn_separation];
			ave_powerIDissipation_s = aIII * interfacePropPtr->interface_scalar_Vals[s1i_dps_separation];

#if DiM2a3_F	
			double aII = interfacePropPtr->interface_scalar_Vals[s1i_absSlip];
			ave_powerIDissipation_s += aII * interfacePropPtr->interface_scalar_Vals[s1i_dps_slip];
#endif
			interfacePropPtr->interface_scalar_Vals[s1i_dpn_mean] = ave_powerIDissipation_n;
			interfacePropPtr->interface_scalar_Vals[s1i_dps_mean] = ave_powerIDissipation_s;

			// now computing effective ones from star values
			double star_powerIDissipation_n, star_powerIDissipation_s = 0.0;
			VEC *sStar, *vStars[NUM_SIDES];
			sStar = &tPtSlns->sl_side_temp_ptData[SDL].sigma_Star;
			vStars[SDL] = &tPtSlns->sl_side_temp_ptData[SDL].v_Star;
			vStars[SDR] = &tPtSlns->sl_side_temp_ptData[SDR].v_Star;
			star_powerIDissipation_n = (*sStar)[0] * ((*(vStars[SDR]))[0] - (*(vStars[SDL]))[0]);
			interfacePropPtr->interface_scalar_Vals[s1i_dpn_star] = star_powerIDissipation_n;
#if DiM2a3_F	
			for (int i = 1; i < DiM; ++i)
				star_powerIDissipation_s += (*sStar)[i] * ((*(vStars[SDR]))[i] - (*(vStars[SDL]))[i]);
			interfacePropPtr->interface_scalar_Vals[s1i_dps_star] = star_powerIDissipation_s;
			interfacePropPtr->interface_scalar_Vals[s1i_dp_star] = star_powerIDissipation_n + star_powerIDissipation_s;
#else
			interfacePropPtr->interface_scalar_Vals[s1i_dps_star] = 0;
			interfacePropPtr->interface_scalar_Vals[s1i_dp_star] = star_powerIDissipation_n;
#endif
		}
	}
	// 4. Taking care of in-situ stress and copying values to permanent storage
	if (interfacePFs->has_in_situ)
	{
		VEC* sStarL = &tPtSlns->sl_side_temp_ptData[SDL].sigma_Star;
		VEC* inSitu = &interfacePFs->in_situ_stress;
		for (int i = 0; i < DiM; ++i)
			(*sStarL)[i] -= (*inSitu)[i];
		CopyVec(*sStarL, tPtSlns->sl_side_temp_ptData[SDR].sigma_Star);
	}
}

void SLInterfaceCalculator::UpdateCharacteristics_FromIncident()
{
#if !USE_ISO_ASSUMPTION
	THROW("The function needs to be updated\n");
#endif
	// wo in-situ is used as these interfaces are fully bonded for incident
	if (incidentSide == ilt_left)
	{

		for (unsigned int i = 0; i < DiM; ++i)
			wlSide_rGoing_WO_in_situ[i] += 2.0 * incientIncomingStress[i];
		return;
	}
	if (incidentSide == ilt_right)
	{

		for (unsigned int i = 0; i < DiM; ++i)
			wrSide_lGoing_WO_in_situ[i] += 2.0 * incientIncomingStress[i];
		return;
	}
	cout << "incidentSide\t" << incidentSide << '\n';
	THROW("Invalid incient option\n");
}

void SLInterfaceCalculator::UpdateStarValues_FromIncident()
{
#if !USE_ISO_ASSUMPTION
	THROW("The function needs to be updated\n");
#endif
	if (incidentSide == ilt_left) 
	{
		VEC* sStarL = &pPtSlns->sl_side_ptData[SDL].sigma_downstream_final;

		for (unsigned int i = 0; i < DiM; ++i)
			(*sStarL)[i] -= incientIncomingStress[i];

		VEC* vStarL = &pPtSlns->sl_side_ptData[SDL].v_downstream_final;
		for (unsigned int i = 0; i < DiM; ++i)
			(*vStarL)[i] += incientIncomingStress[i] * ts_bulkProps->bulk_leftPtr->inv_c_rhos[i];
	}
	else if (incidentSide == ilt_right)
	{
		VEC* sStarR = &pPtSlns->sl_side_ptData[SDR].sigma_downstream_final;

		for (unsigned int i = 0; i < DiM; ++i)
			(*sStarR)[i] -= incientIncomingStress[i];

		VEC* vStarR = &pPtSlns->sl_side_ptData[SDR].v_downstream_final;
		for (unsigned int i = 0; i < DiM; ++i)
			(*vStarR)[i] -= incientIncomingStress[i] * ts_bulkProps->bulk_rightPtr->inv_c_rhos[i];
	}
}

void SLInterfaceCalculator::Update_Velocity_Left_Ring_Opened1D()
{
#if RING_PROBLEM
	return;
#endif
	VEC *vel;
	if (b_directlyCalculateBondedInFinalSln)
		vel = &pPtSlns->sl_side_ptData[SDR].v_downstream_final;
	else
	{
		vel = &tPtSlns->sl_side_temp_ptData[SDR].v_Star;
		VEC* u = &tPtSlns->sl_side_temp_ptData[SDR].u_downstream_latestValue;
		(*u)[0] -= g_domain->ring_opened1D_al * pPtSlns->interface_time;
	}
	(*vel)[0] -= g_domain->ring_opened1D_al;
}

void SLInterfaceCalculator::Compute_1ConvergedImplcit_vsStar_TSR_OrtizSimple_Calculation()
{
	ortiz1DTSRptr->ws[SDL] = wlSide_rGoing[0];
	ortiz1DTSRptr->ws[SDR] = wrSide_lGoing[0];
	SL_Bulk_Properties *bulk_left_outPtr, *bulk_right_outPtr;
	ts_bulkProps->Get_Bulks2Side(bulk_left_outPtr, bulk_right_outPtr);
	ortiz1DTSRptr->Zs[SDL] = bulk_left_outPtr->c_rhos[0];
	ortiz1DTSRptr->Zs[SDR] = bulk_right_outPtr->c_rhos[0];

//	ortiz1DTSRptr->Zeff = 1.0 / ts_bulkProps->YAve_n;
	ortiz1DTSRptr->sigmaC = tPtSlns->sigmaC;
	ortiz1DTSRptr->deltaC = tPtSlns->deltaC;
	double sI, sStar, uLStar, vLStar, uRStar, vRStar, maxDelu;
	TSR_loadingStages solutionStage;
	TSR_OrtizStages stageShort;
#if DB_ON
	if (b_db_p)
		ortiz1DTSRptr->Write_MainInputs(dbout);
//	if (g_time > 7.0)	{		fstream out("TestFiles/TestOrtiz1D_debug.txt", ios::out);		ortiz1DTSRptr->Write_MainInputs(out);		out.close();	}
#endif
	sI = ortiz1DTSRptr->Compute_step_np1_Solution(sStar, uLStar, vLStar, uRStar, vRStar, maxDelu, solutionStage, stageShort);
	
	tPtSlns->sl_side_temp_ptData[SDL].sigma_mode_Star[rmode_stick][0] = sI;
	tPtSlns->sl_side_temp_ptData[SDR].sigma_mode_Star[rmode_stick][0] = sI;

	double delu = uRStar - uLStar;
	tPtSlns->sl_side_temp_ptData[SDL].sigma_Star[0] = sStar;				tPtSlns->sl_side_temp_ptData[SDR].sigma_Star[0] = sStar;
	tPtSlns->sl_side_temp_ptData[SDL].v_Star[0] = vLStar;					tPtSlns->sl_side_temp_ptData[SDR].v_Star[0] = vRStar;
	tPtSlns->sl_side_temp_ptData[SDL].u_downstream_latestValue[0] = uLStar;	tPtSlns->sl_side_temp_ptData[SDR].u_downstream_latestValue[0] = uRStar;
	pPtSlns->maxEffDelU = maxDelu;

	bool isSep = ((sStar > 0) || (delu > ortiz1DTSRptr->toldc));
	tPtSlns->damageOffOnMix = sl_interfacial_damage_on;
	tPtSlns->contactOffOnMix = sl_contact_off;
	tPtSlns->slipOffOnMix = sl_slip_off;

	// the reason this is not stored in "interfacePropPtr"is that it can be used to make the pointer null or not
	tPtSlns->interface_damageLatestValue = 1.0;
	tPtSlns->tsrStageLatestValue = ortiz1DTSRptr->step_np1_tsrStage;

	if (interfacePFs->HasAnyContactFriction())
	{
		if (tPtSlns->interfacePropPtr == NULL)
			tPtSlns->interfacePropPtr = new SL_Interface_Temp_PtData_IntContFrac();

		tPtSlns->interfacePropPtr->del_u_nt_parts.vec_n = delu;
		tPtSlns->interfacePropPtr->del_u_nt_parts.vec_n_pos = MAX(delu, 0.0);

		tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_D] = 1.0;
		tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_effectiveStrs] = sI;
		double delv = vRStar - vLStar;
		tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_magnitude_delv] = delv;
		double dissPower = delv * sStar;
		tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_dpn_mean] = dissPower;
		tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_dpn_star] = dissPower;
		tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_dp_star] = dissPower;


		if (isSep)
		{
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_absStick] = 0.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_contactRel] = 0.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_stickRel] = 0.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_absSep] = 1.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_dpn_separation] = dissPower;
		}
		else
		{
			tPtSlns->contactOffOnMix = sl_contact_on;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_absStick] = 1.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_contactRel] = 1.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_stickRel] = 1.0;
			tPtSlns->interfacePropPtr->interface_scalar_Vals[s1i_absSep] = 0.0;
		}
	}

#if DB_ON
	if (b_db_p)
	{
		dbout << ortiz1DTSRptr->toldc << "toldc\n";
		dbout << "solutionStage\t" << solutionStage << '\n';
		dbout << "stageShort\t" << stageShort << '\n';
		dbout << "un\t";
		dbout << ortiz1DTSRptr->step_n.sides[SDL].u << '\t' << ortiz1DTSRptr->step_n.sides[SDR].u << '\n';
		dbout << "vn\t";
		dbout << ortiz1DTSRptr->step_n.sides[SDL].v << '\t' << ortiz1DTSRptr->step_n.sides[SDR].v << '\n';
		dbout << "w\t";
		dbout << ortiz1DTSRptr->ws[SDL] << '\t' << ortiz1DTSRptr->ws[SDR] << '\n';
		dbout << "maxdelun\t";
		dbout << ortiz1DTSRptr->step_n_max_delu << '\n';
		dbout << "tsrSn\t";
		dbout << ortiz1DTSRptr->step_n_tsrStage << '\n';
		dbout << "sigmaC\t";
		dbout << ortiz1DTSRptr->sigmaC << '\n';
		dbout << "delT\t";
		dbout << ortiz1DTSRptr->delT << '\n';
		dbout << "sI\t";
		dbout << sI << '\n';
		dbout << "sStar\t";
		dbout << sStar << '\n';
		dbout << "maxDelu\t";
		dbout << maxDelu << '\t';
		dbout << "pPtSlns->maxEffDelU\t" << pPtSlns->maxEffDelU  << '\t';
		if (interfacePFs->HasAnyContactFriction())
			dbout << "tPtSlns\t" << tPtSlns->interfacePropPtr->del_u_nt_parts.getEffective();
		dbout << '\n';
		dbout << "uStar(L,R,del)\t";
		dbout << uLStar << '\t' << uRStar << '\t' << uRStar - uLStar << '\n';
		dbout << "vStar(L,R,del)\t";
		dbout << vLStar << '\t' << vRStar << '\t' << vRStar - vLStar << '\n';
	}
#endif

}

void SLInterfaceCalculator::Compute_Star_sv_from_RiemannItoIII_toStar()
{
	if (b_directlyCalculateBondedInFinalSln)
	{
		Compute_Riemann_I_Bonded();
		return;
	}
	tPtSlns->damageOffOnMix = interfacePFs->damageOffOnMix;
	tPtSlns->contactOffOnMix = interfacePFs->contactOffOnMix;
	tPtSlns->slipOffOnMix = interfacePFs->slipOffOnMix;

	// step 1: setting pointer to Riemann values and absolute areas
	VEC *sStars[NUM_SIDES], *vStars[NUM_SIDES];
	sStars[SDL] = &tPtSlns->sl_side_temp_ptData[SDL].sigma_Star;
	sStars[SDR] = &tPtSlns->sl_side_temp_ptData[SDR].sigma_Star;
	vStars[SDL] = &tPtSlns->sl_side_temp_ptData[SDL].v_Star;
	vStars[SDR] = &tPtSlns->sl_side_temp_ptData[SDR].v_Star;

	SolidInterfaceShortType mode;
	VEC* sRs[3];
	double tmp;
	// indexed as [indexAorB][ri]
	VEC* vRs[NUM_SIDES][3];

	// relative areas for bonded (stick), slip, and separation 
	double aI, aII, aIII;

	sRs[s1i_absStick] = &tPtSlns->sl_side_temp_ptData[SDL].sigma_mode_Star[rmode_stick];
	sRs[s1i_absSep] = &tPtSlns->sl_side_temp_ptData[SDL].sigma_mode_Star[rmode_sep];

	vRs[SDL][s1i_absStick] = &tPtSlns->sl_side_temp_ptData[SDL].v_mode_Star[rmode_stick];
	vRs[SDL][s1i_absSep] = &tPtSlns->sl_side_temp_ptData[SDL].v_mode_Star[rmode_sep];

	vRs[SDR][s1i_absStick] = &tPtSlns->sl_side_temp_ptData[SDR].v_mode_Star[rmode_stick];
	vRs[SDR][s1i_absSep] = &tPtSlns->sl_side_temp_ptData[SDR].v_mode_Star[rmode_sep];
#if DiM2a3_F
	sRs[s1i_absSlip] = &tPtSlns->sl_side_temp_ptData[SDL].sigma_mode_Star[rmode_slip];
	vRs[SDL][s1i_absSlip] = &tPtSlns->sl_side_temp_ptData[SDL].v_mode_Star[rmode_slip];
	vRs[SDR][s1i_absSlip] = &tPtSlns->sl_side_temp_ptData[SDR].v_mode_Star[rmode_slip];
#endif
	// step 1.I Computing I (bonded + stick) solution
	Compute_Riemann_I_Bonded();

	// if D == 0
	// no fracture: Nothing further to be done
	if (!interfacePFs->HasAnyContactFriction()) // damageOffOnMix == sl_interfacial_damage_off)
	{
		mode = s1i_absStick;
		CopyVec(*(sRs[mode]), *(sStars[SDL]));
		CopyVec(*(sRs[mode]), *(sStars[SDR]));
		CopyVec(*(vRs[SDL][mode]), *(vStars[SDL]));
		CopyVec(*(vRs[SDR][mode]), *(vStars[SDR]));
		return;
	}

	// if pointer of scalar values is null, it's made, del_u, sigma_I nt values computed, and del_v nt processed is set to false
	double delu0Change = 0.0;
	if ((oneIntAlltimesPtr != NULL) && (oneIntAlltimesPtr->has_ring_opened1D_al))
		delu0Change = g_domain->ring_opened1D_al * pPtSlns->interface_time;

	tPtSlns->MakeReady_For_ContactFractureRuns(interfacePFs->beta_delU, interfacePFs->beta_traction, *(sRs[s1i_absStick]), delu0Change);
	SL_Interface_Temp_PtData_IntContFrac* interfacePropPtr = tPtSlns->interfacePropPtr;
	interfacePropPtr->interface_scalar_Vals[s1i_fricCoef] = interfacePFs->Get_Const_friction_Coef();

	// from here 0 < D < 1 or D = 1
	if (1.0 - tPtSlns->interface_damageLatestValue <= interfacePFs->zeroTol4DamageEq1)
	{
		tPtSlns->interface_damageLatestValue = 1.0;
		tPtSlns->damageOffOnMix = sl_interfacial_damage_on;
	}
	double D = tPtSlns->interface_damageLatestValue;
	double omD = 1.0 - D;
	interfacePropPtr->interface_scalar_Vals[s1i_D] = D;

	// step 1.III computing separation solution
	double relCont;
	tPtSlns->contactOffOnMix = Compute_Riemann_III_Separation(tPtSlns->damageOffOnMix, relCont);
	interfacePropPtr->interface_scalar_Vals[s1i_contactRel] = relCont;

	if (tPtSlns->contactOffOnMix == sl_contact_off) // so aRelContact is zero and all the interface is between separation and bonded (stick) modes
	{
		// relCont == 0
		tPtSlns->slipOffOnMix = sl_slip_off;
		interfacePropPtr->interface_scalar_Vals[s1i_stickRel] = 0.0;

		if (tPtSlns->damageOffOnMix == sl_interfacial_damage_on) // fully damaged and entirely in separation mode
		{
			mode = s1i_absSep;
			CopyVec(*(sRs[mode]), *(sStars[SDL]));
			CopyVec(*(sRs[mode]), *(sStars[SDR]));
			CopyVec(*(vRs[SDL][mode]), *(vStars[SDL]));
			CopyVec(*(vRs[SDR][mode]), *(vStars[SDR]));

			interfacePropPtr->interface_scalar_Vals[s1i_absStick] = 0.0;
			interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = 0.0;
			interfacePropPtr->interface_scalar_Vals[s1i_absSep] = 1.0;
			return;
		}
		// now 0 < D < 1, aRelContact (relCont) = 0, aI = 1 - D, aII = 0, aIII = D -> interpolation between stick and separation modes
		interfacePropPtr->interface_scalar_Vals[s1i_absStick] = omD;
		interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = 0.0;
		interfacePropPtr->interface_scalar_Vals[s1i_absSep] = D;

		for (int i = 0; i < DiM; ++i)
		{
			tmp = omD * (*(sRs[s1i_absStick]))[i] + D * (*(sRs[s1i_absSep]))[i];
			(*(sStars[SDL]))[i] = tmp;
			(*(sStars[SDR]))[i] = tmp;
			tmp = omD * (*(vRs[SDL][s1i_absStick]))[i] + D * (*(vRs[SDL][s1i_absSep]))[i];
			(*(vStars[SDL]))[i] = tmp;
			tmp = omD * (*(vRs[SDR][s1i_absStick]))[i] + D * (*(vRs[SDR][s1i_absSep]))[i];
			(*(vStars[SDR]))[i] = tmp;
		}
		return;
	}
	// from here there is some contact (relCont > 0) [Already know D >= 0]
	// Short: D != 0, aRelativeContact (relCont) != 0

	// relative area of contact
	// check on whether we have slip or stick:
	// step 1.II computing slip
	onOffT slip_bOnOffMixed = allOffV;
	double relStick = 1.0;
	if (tPtSlns->slipOffOnMix != sl_slip_off)
		slip_bOnOffMixed = Compute_Riemann_II_Slip(relStick, tPtSlns->slipOffOnMix);
	interfacePropPtr->interface_scalar_Vals[s1i_stickRel] = relStick;

	// there are two cases, slip on or off
	if (slip_bOnOffMixed == allOffV)
	{
		// lambda = 1 (all of contact is stick) -> 
		//	aI = 1 - D + D * relCont3
		//	aII = 0
		//  aIII = D * (1 - relCont)
		// two cases 
		// all stick
		if (tPtSlns->contactOffOnMix == sl_contact_on) // all in contact and in contact all in stick 
		{
			// relCont = 1 -> aI = 1, aII = 0, aIII = 0
			interfacePropPtr->interface_scalar_Vals[s1i_absStick] = 1.0;
			interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = 0.0;
			interfacePropPtr->interface_scalar_Vals[s1i_absSep] = 0.0;

			mode = s1i_absStick;
			CopyVec(*(sRs[mode]), *(sStars[SDL]));
			CopyVec(*(sRs[mode]), *(sStars[SDR]));
			CopyVec(*(vRs[SDL][mode]), *(vStars[SDL]));
			CopyVec(*(vRs[SDR][mode]), *(vStars[SDR]));
		}
		else // 0 < relCont < 1, general case above -> stick - separation
		{
			aI = omD + D * relCont;
			aII = 0;
			aIII = D * (1 - relCont);
			interfacePropPtr->interface_scalar_Vals[s1i_absStick] = aI;
			interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = aII;
			interfacePropPtr->interface_scalar_Vals[s1i_absSep] = aIII;

			for (int i = 0; i < DiM; ++i)
			{
				tmp = aI * (*(sRs[s1i_absStick]))[i] + aIII * (*(sRs[s1i_absSep]))[i];
				(*(sStars[SDL]))[i] = tmp;
				(*(sStars[SDR]))[i] = tmp;
				tmp = aI * (*(vRs[SDL][s1i_absStick]))[i] + aIII * (*(vRs[SDL][s1i_absSep]))[i];
				(*(vStars[SDL]))[i] = tmp;
				tmp = aI * (*(vRs[SDR][s1i_absStick]))[i] + aIII * (*(vRs[SDR][s1i_absSep]))[i];
				(*(vStars[SDR]))[i] = tmp;
			}
		}
		return;
	}
	else if (slip_bOnOffMixed == allOnV) // the case that within contact all is slip
	{
		// aI = 1 - D;
		// aII = D * relCont;
		// aIII = D * (1 - relCont);

		// two cases, fully damaged or not
		if (tPtSlns->damageOffOnMix == sl_interfacial_damage_on) // D = 1  (slip - separation)
		{
			// aI = 0, aII = relCont, aIII = 1 - relCont
			if (tPtSlns->contactOffOnMix == sl_contact_on) // all conact -> relCont = 1
			{
				// aI = 0, aII = 1, aIII = 0;
				mode = s1i_absSlip;
				CopyVec(*(sRs[mode]), *(sStars[SDL]));
				CopyVec(*(sRs[mode]), *(sStars[SDR]));
				CopyVec(*(vRs[SDL][mode]), *(vStars[SDL]));
				CopyVec(*(vRs[SDR][mode]), *(vStars[SDR]));

				interfacePropPtr->interface_scalar_Vals[s1i_absStick] = 0.0;
				interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = 1.0;
				interfacePropPtr->interface_scalar_Vals[s1i_absSep] = 0.0;
				return;
			}
			if (tPtSlns->contactOffOnMix == sl_contact_mixed) // general case, interpolation between slip and separation
			{
				aI = 0, aII = relCont, aIII = 1 - relCont;
				interfacePropPtr->interface_scalar_Vals[s1i_absStick] = aI;
				interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = aII;
				interfacePropPtr->interface_scalar_Vals[s1i_absSep] = aIII;

				for (int i = 0; i < DiM; ++i)
				{
					tmp = aII * (*(sRs[s1i_absSlip]))[i] + aIII * (*(sRs[s1i_absSep]))[i];
					(*(sStars[SDL]))[i] = tmp;
					(*(sStars[SDR]))[i] = tmp;
					tmp = aII * (*(vRs[SDL][s1i_absSlip]))[i] + aIII * (*(vRs[SDL][s1i_absSep]))[i];
					(*(vStars[SDL]))[i] = tmp;
					tmp = aII * (*(vRs[SDR][s1i_absSlip]))[i] + aIII * (*(vRs[SDR][s1i_absSep]))[i];
					(*(vStars[SDR]))[i] = tmp;
				}
				return;
			}
			if (tPtSlns->contactOffOnMix == sl_contact_off) // This case should have been handled before
			{
				THROW("Contact is off! We should not have gotten to this point in the function.\n")
			}
		}
		else // 0 < D < 1, lamba = 0	-> aI = 1 - D, aII =  D * relCont, aIII = D * (1 - relCont)
		{
			interfacePropPtr->interface_scalar_Vals[s1i_absStick] = omD;
			if (tPtSlns->contactOffOnMix == sl_contact_on) // all conact -> relCont = 1
			{
				aI = 1 - D, aII = D, aIII = 0;
				interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = D;
				interfacePropPtr->interface_scalar_Vals[s1i_absSep] = 0;

				for (int i = 0; i < DiM; ++i)
				{
					tmp = aI * (*(sRs[s1i_absStick]))[i] + aII * (*(sRs[s1i_absSlip]))[i];
					(*(sStars[SDL]))[i] = tmp;
					(*(sStars[SDR]))[i] = tmp;
					tmp = aI * (*(vRs[SDL][s1i_absStick]))[i] + aII * (*(vRs[SDL][s1i_absSlip]))[i];
					(*(vStars[SDL]))[i] = tmp;
					tmp = aI * (*(vRs[SDR][s1i_absStick]))[i] + aII * (*(vRs[SDR][s1i_absSlip]))[i];
					(*(vStars[SDR]))[i] = tmp;
				}
			}
			else if (tPtSlns->contactOffOnMix == sl_contact_mixed) // general case
			{
				// 0 < D < 1, lambda = 0, relCont != 0 or 1
				aI = omD, aII = D * relCont, aIII = D * (1 - relCont);
				interfacePropPtr->interface_scalar_Vals[s1i_absSlip] = aII;
				interfacePropPtr->interface_scalar_Vals[s1i_absSep] = aIII;

				for (int i = 0; i < DiM; ++i)
				{
					tmp = aI * (*(sRs[s1i_absStick]))[i] + aII * (*(sRs[s1i_absSlip]))[i] + aIII * (*(sRs[s1i_absSep]))[i];
					(*(sStars[SDL]))[i] = tmp;
					(*(sStars[SDR]))[i] = tmp;
					tmp = aI * (*(vRs[SDL][s1i_absStick]))[i] + aII * (*(vRs[SDL][s1i_absSlip]))[i] + aIII * (*(vRs[SDL][s1i_absSep]))[i];
					(*(vStars[SDL]))[i] = tmp;
					tmp = aI * (*(vRs[SDR][s1i_absStick]))[i] + aII * (*(vRs[SDR][s1i_absSlip]))[i] + aIII * (*(vRs[SDR][s1i_absSep]))[i];
					(*(vStars[SDR]))[i] = tmp;
				}
			}
			else if (tPtSlns->contactOffOnMix == sl_contact_off) // This case should have been handled before
			{
				THROW("Contact is off! We should not have gotten to this point in the function.\n")
			}
			return;
		}
	}
	else if (slip_bOnOffMixed == partialV)
	{
		if (interfacePFs->slipOffOnMix != sl_slip_mixed)
			THROW("(slip_bOnOffMixed != sl_slip_mixed)\n");
		THROW("We don't interpolate between stick and slip solutions, so slip_bOnOffMixed == sl_slip_mixed is not implemented now.\n");
	}
}

void SLInterfaceCalculator::SetStressComputeVelocityStarFromStressStar(const VEC& sigmaStar, RiemannMode_StorageT rMode)
{
	for (int sidei = 0; sidei < NUM_SIDES; ++sidei)
		CopyVec(sigmaStar, tPtSlns->sl_side_temp_ptData[sidei].sigma_mode_Star[rMode]);
	ts_bulkProps->Compute_vStars_from_sigmaStar_ws(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, sigmaStar,
		tPtSlns->sl_side_temp_ptData[SDL].v_mode_Star[rMode],
		tPtSlns->sl_side_temp_ptData[SDR].v_mode_Star[rMode]);
}

bool SLInterfaceCalculator::DirectlyCalculateBondedInFinalSln() const
{
	return ((ts_bulkProps->interfaceLoc != ilt_twoSided) || (!interfacePFs->HasAnyContactFriction()));
}

// const VEC& uL, const VEC& uR, 
// beta_delU = interfacePFs->beta_delU;
// VEC delU;
// for (int i = 0; i < DiM; ++i)
//	delU[i] = uR[i] - uL[i];
// del_u_nt_parts.Compute_derivedValues_From_vec(delU, beta_delU);

// double beta_traction = interfacePFs->beta_traction;
// sigma_I_nt_parts.Compute_derivedValues_From_vec(sigmaI, beta_traction);

void Compute_Riemann_III_Separation_Damage_Rate_Aux(SLFF_InterfacialDamageModeType damageOffOnMix, SLFF_ContactType contactOffOnMix, Vec_nt& sigma_I_nt_parts, Vec_nt& del_u_nt_parts, const VEC& vL, const VEC& vR, double maxEffDelU, TSR_loadingStages tsrStagePrevStepValue, TSR_loadingStages& tsrStageLatestValue, double damage,
	double sigmaC, double deltaC, double tau, const SL_Interface_Fracture_PF* interfacePFs, double dn, double dt, 
	Vec_nt& del_v_nt_parts, double& effS, double& DDot, double& DTarget, double& DSrc_Allix, double & DSrc_RBH, VEC& sigmaIII, bool& terminateConditionMet)
{
	DDot = 0.0;
	DTarget = damage;
	DSrc_Allix = 0.0;
	DSrc_RBH = 0.0;
	double strengthNormalizer;
	setValue(sigmaIII, 0.0);
	effS = Compute_Interface_EffectiveStress(interfacePFs, sigmaC, sigma_I_nt_parts, strengthNormalizer);
	if (contactOffOnMix == sl_contact_on)
		return;
	double kCompressive = 10.0 * sigmaC / deltaC;
//	if (interfacePFs->damageTractionModel == sl_interfacial_damage_traction_zero)
	if (interfacePFs->damageTractionModel == sl_interfacial_damage_traction_TSR)
	{
		ComputeTSR_Tractions(interfacePFs->tsrModel, del_u_nt_parts, sigma_I_nt_parts, interfacePFs->beta_delU, deltaC, sigmaC, kCompressive, sigmaIII, maxEffDelU, tsrStagePrevStepValue, tsrStageLatestValue, terminateConditionMet);
		DBCHK(dbout << "sIII\t" << sigmaIII[0] << '\n';);
		DBCHK(dbout << "tsrStageLatestValueFinal\t" << tsrStageLatestValue << '\n';);
	}

	// now computing damage rate
	if ((damageOffOnMix != sl_interfacial_damage_mixed) || (damage >= interfacePFs->mxV4SrcZero))
		return;

	// computing damage rate:
	//! 1. Allix's part
	DamageDTargetRelT d_dTargetRel;
	double pos_damTarget_minus_dam = 0.0;
	double damageBounded = MIN(damage, 1.0);

	// computing target damage 
	double c = interfacePFs->c0;
	if (!interfacePFs->cConst)
	{
		double omk = 1.0 - interfacePFs->kc;
		if (interfacePFs->pc == 1)
			c = (omk * (1.0 - damageBounded) + interfacePFs->kc) * interfacePFs->c0;
		else
		{
			double omd = 1.0 - damageBounded;
			c = (omk * pow(omd, interfacePFs->pc) + interfacePFs->kc) * interfacePFs->c0;
		}
	}
	if (!interfacePFs->nConst)
	{
		double omk = 1.0 - interfacePFs->kn;
		if (interfacePFs->pn == 1)
			c = (omk * (1.0 - damageBounded) + interfacePFs->kn) * interfacePFs->c0;
		else
		{
			double omd = 1.0 - damageBounded;
			c = (omk * pow(omd, interfacePFs->pn) + interfacePFs->kn) * interfacePFs->c0;
		}
	}
	double denomInverse = 1.0 / (1.0 - c);
	double xn = effS / strengthNormalizer;
	double xp = (xn - c) * denomInverse;
	int xp_Status = 1;
	if (xp < 0)
	{
		DTarget = 0.0;
		d_dTargetRel = damage_target_zero;
	}
	else
	{
		DTarget = xp;
		if (interfacePFs->fx_boundedBy1)
			DTarget = MIN(DTarget, 1.0);
		pos_damTarget_minus_dam = DTarget - damageBounded;
		if (pos_damTarget_minus_dam > 0.0)
			d_dTargetRel = damage_target_gtt_damage;
		else
		{
			d_dTargetRel = damage_target_lesst_damage;
			pos_damTarget_minus_dam = 0.0;
		}
	}

	DSrc_Allix = 0.0;
	if ((interfacePFs->hfunT == AllixHExp) && (d_dTargetRel == damage_target_gtt_damage))
	{
		double expt = exp(-interfacePFs->Hpara0 * pos_damTarget_minus_dam);
		DSrc_Allix = 1.0 / tau * (1.0 - expt);
	}

	if (interfacePFs->rbhvDTT == RBH_VelDDot_none)
	{
		DSrc_RBH = 0.0;
		return;
	}
	VEC delV;
	for (int i = 0; i < DiM; ++i)
		delV[i] = vR[i] - vL[i];
	double beta_delU = interfacePFs->beta_delU;
	del_v_nt_parts.Compute_derivedValues_From_vec(delV, beta_delU);
	// scales:

	// displacement scales for Allix's model
	double dvn = dn * interfacePFs->dvnRatio;
	double dvnInv = 1.0 / dvn;

	if (interfacePFs->rbhvDTT == RBH_VelDDot_vnp_only)
		DSrc_RBH = dvnInv * del_v_nt_parts.vec_n_pos;
#if DiM2a3_F
	else
	{
		if (interfacePFs->rbhvDTT == RBH_VelDDot_vel_effective)
			DSrc_RBH = dvnInv * del_v_nt_parts.getEffective();
		else
		{
			double dvt = dt * interfacePFs->dvtRatio;

			if (interfacePFs->rbhvDTT == RBH_VelDDot_vtabs_only)
				DSrc_RBH = (1.0 / dvt) * del_v_nt_parts.vec_t_abs;
			else if (interfacePFs->rbhvDTT == RBH_VelDDot_vnp_vtabs)
				DSrc_RBH = dvnInv * del_v_nt_parts.vec_n_pos + (1.0 / dvt) * del_v_nt_parts.vec_t_abs;
		}
	}
#endif
}

void ComputeTSR_Tractions(SLFF_TSRType tsrModel, Vec_nt& del_u_nt_parts, Vec_nt& sigma_I_nt_parts, double beta_delU, double deltaC, double sigmaC, double kCompressive, VEC& S, double max_delU, TSR_loadingStages tsrStagePrevStepValue, TSR_loadingStages& tsrStageLatestValue, bool& terminateConditionMet)
{
	if (tsrModel == tsr_Zero)
	{
		setValue(S, 0.0);
		return;
	}
#if DiM3
	THROW("Update the function\n");
#endif
	if (tsrModel == tsr_Xu_Needleman)
	{
		double dup0 = del_u_nt_parts.vec_n;
		double dup1 = 0;
#if DiM2a3_F
#if DiM2	
		dup1 = del_u_nt_parts.vec_t[0];
#else
		dup1 = del_u_nt_parts.vec_t_abs;
#endif
#endif
		double normDelU = dup0 * dup0 + dup1 * dup1;
		normDelU = sqrt(normDelU);
		double factor;
		bool beforePeakLoad = ((max_delU > 0.0) && (normDelU < max_delU));
		if (beforePeakLoad)
		{
			factor = normDelU / max_delU;
			dup0 /= factor;
#if DiM2a3_F
			dup1 /= factor;
#endif
		}
		dup0 /= deltaC;
		dup1 /= deltaC;
		double expp = 1.0 - dup0 - dup1 * dup1;
		double expv = exp(expp);
		S[0] = sigmaC * dup0 * expv;
#if DiM2a3_F
		S[1] = 2.0 * sigmaC * dup1 * (1.0 + dup0) * expv;
#endif	
		if ((max_delU > 1.0) && (Norm2(S) / sigmaC <= 0.001))
			terminateConditionMet = true;

		if (beforePeakLoad)
		{
			S[0] *= factor;
#if DiM2a3_F
			S[1] *= factor;
#endif
		}
		return;
	}
	if (tsrModel == tsr_Ortiz)
	{
		DBCHK(dbout << "sigma_I_nt_parts(n,+)\t" << sigma_I_nt_parts.vec_n << "\t" << sigma_I_nt_parts.vec_n_pos << '\n';);
		DBCHK(dbout << "del_u_nt_parts(n,+)\t" << del_u_nt_parts.vec_n << "\t" << del_u_nt_parts.vec_n_pos << '\n';);
		DBCHK(dbout << "max_delU\t" << max_delU << '\n';);
		DBCHK(dbout << "sigmaC\t" << sigmaC << '\n';);
		DBCHK(dbout << "deltaC\t" << deltaC << '\n';);
		DBCHK(dbout << "kCompressive\t" << kCompressive << '\n';);
		DBCHK(dbout << "tsrStagePrevStepValue\t" << tsrStagePrevStepValue << '\n';);
		DBCHK(dbout << "tsrStageLatestValue\t" << tsrStageLatestValue << '\n';);
		
		// fully damaged
		if (max_delU > deltaC)
		{
			setValue(S, 0.0);
			DBCHK(dbout << "delU>delC-S0\n";);

			if (del_u_nt_parts.vec_n < 0.0) // want to get out of penetration mode
			{
				if (sigma_I_nt_parts.vec_n < 0.0) // if sigmI > 0, since s* (S) == 0, we'll get out of penetration, so need to only take care of this for sigmaI < 0
					S[0] = MIN(sigma_I_nt_parts.vec_n, del_u_nt_parts.vec_n * kCompressive);
				tsrStageLatestValue = tsr1_dam1_delun_neg;
				DBCHK(dbout << "delunNeg\n";);
			}
			else
			{
				terminateConditionMet = true;
				tsrStageLatestValue = tsr1_dam1_delun_pos;
			}
			return;
		}
	//	Pandolfi_Ortiz_1999_Finite element simulation of ring expansion and fragmentation The capturing of length and time scales through cohesive models of fracture.pdf
		if (max_delU < 1e-10) // no opening
		{
			DBCHK(dbout << "max_delU0\n";);

			double effS;
#if DiM1
			effS = sigma_I_nt_parts.vec_n_pos;
#else
			double tmp = sigma_I_nt_parts.vec_t_abs / beta_delU;
			effS = sqrt(sigma_I_nt_parts.vec_n_pos * sigma_I_nt_parts.vec_n_pos + tmp * tmp);
#endif
			double ratio = effS / sigmaC;

			DBCHK(dbout << "effS\t" << effS << "\tratio\t" << ratio << '\n';);

			if (ratio > 1.0)
			{
				double inv_ratio = 1.0 / ratio;
				for (int i = 0; i < DiM; ++i)
					S[i] = inv_ratio * sigma_I_nt_parts.vec[i];
				if (sigma_I_nt_parts.vec_n < 0)
					S[0] = sigma_I_nt_parts.vec_n;
				tsrStageLatestValue = tsrl_dam0_sigI_over_sigC;
				DBCHK(dbout << "ratio>1\tinv_ratio\t" << inv_ratio << '\n';);
			}
			else
			{
				CopyVec(sigma_I_nt_parts.vec, S);
				if (sigma_I_nt_parts.vec_n > 0)
					tsrStageLatestValue = tsrl_dam0_sigI_less_sigC;
				else
					tsrStageLatestValue = tsrl_dam0_sigI_neg;
			}
//			if (del_u_nt_parts.vec_n < 0)
//				S[0] = MIN(S[0], del_u_nt_parts.vec_n * kCompressive);
			return;
		}
		// now the cohesive interface is inserted and opening is experienced
		double eff_delu = del_u_nt_parts.getEffective();
		DBCHK(dbout << "eff_delu\t" << eff_delu << '\n';);

		// this is at some fraction (< 1) of max_delU to have this region below which it's definitely unloading
		// above which opening state determines the branch
		double tol_delc = 0.04 * deltaC;
		double delta1 = MAX(max_delU - tol_delc, 0.0);
		DBCHK(dbout << "tol_delc\t" << tol_delc << '\n';);
		DBCHK(dbout << "delta1\t" << delta1 << '\n';);

		int isLoading = -1; // undecided | 1 means it's on loading branch, 0 on unloading
		if (eff_delu > max_delU)
		{
			isLoading = 1;
			DBCHK(dbout << "delU>maxDelU(L=1)\n";);
		}
		else
		{
			if ((eff_delu < delta1) && (eff_delu >= tol_delc))
			{
				isLoading = 0;
				DBCHK(dbout << "delU<delu1,>tolDelC(L=0)\n";);
			}
		}
		double scalar_stress_t;
		if (isLoading != 0) // so it's loading (eff_delu > max_delU) OR it's in undecided region
		{
			double rel_eff_delu = eff_delu / deltaC;
			DBCHK(dbout << "(L!=0)\trel_eff_delu\t" << rel_eff_delu << "\n";);

			// using the linear model
			if (rel_eff_delu < 1.0)
			{
				scalar_stress_t = sigmaC * (1.0 - rel_eff_delu);
				DBCHK(dbout << "(effDelU<1)scalar_stress_t\t" << scalar_stress_t << "\n";);
			}
			else
			{
				tsrStageLatestValue = tsr1_dam1_delun_pos;
				setValue(S, 0.0);
				DBCHK(dbout << "(effDelU>1)S0\n";);
				if (del_u_nt_parts.vec_n < 0)
				{
					S[0] = MIN(S[0], del_u_nt_parts.vec_n * kCompressive);
					tsrStageLatestValue = tsr1_dam1_delun_neg;
					DBCHK(dbout << "del_u_nt_parts.vec_n < 0,throw\n";);
					THROW("For 1D we shouldn't get here. For 2D, this may result in unphysical response. Check if the statement after if makes sense. I may be reasonable.\n");
				}
				return;
			}
			if (isLoading == -1) // still not determined
			{
				// computing velocity jump:
				bool isOpeningMode = ((sigma_I_nt_parts.vec_n - scalar_stress_t) > 0.0);
				DBCHK(dbout << "sI\t" << sigma_I_nt_parts.vec_n << "\tscalar_stress_t\t" << scalar_stress_t << '\n';);
				DBCHK(dbout << "(L==-1)\tisOpeningMode" << isOpeningMode << '\n';);
				if (!isOpeningMode)
				{
					isLoading = 0;
					if (max_delU < eff_delu)
					{
						DBCHK(dbout << "max_delU->updatd\tmax_delU\t" << max_delU << "\teff_delu\t" << eff_delu << '\n';);
						max_delU = eff_delu;
					}
					DBCHK(dbout << "(L->0)isOpeningMode\n";);
				}
			}
		}	
		// not-loading
		double scalar_stress_t_unloading_solution = -1;
		if (isLoading != 1)
		{
			double ratio_eff_delu = eff_delu / max_delU;
			scalar_stress_t_unloading_solution = ratio_eff_delu * (sigmaC * (1.0 - max_delU / deltaC));
			DBCHK(dbout << "(L!=1B)\tratio_eff_delu\t" << ratio_eff_delu << "\tsclrS\t" << scalar_stress_t_unloading_solution << '\n';);
		}
		if (isLoading == 0)
		{
			scalar_stress_t = scalar_stress_t_unloading_solution;
			DBCHK(dbout << "CHKscalar_stress_t\t" << scalar_stress_t << '\n';);
		}
		else if (isLoading == -1)
		{
			DBCHK(dbout << "FL-1\ttsrStagePrevStepValue\t" << tsrStagePrevStepValue << '\n';);

			// both solutions are valid, find the one that is closer to the previous solution
			double scalar_stress_t_loading_solution = scalar_stress_t;
#if !DiM1
			THROW("This impementation is only valid for 1D\n");
#endif
			if ((tsrStagePrevStepValue == tsrl_damMid_unloadingBranchPosDelv) || (tsrStagePrevStepValue == tsrl_damMid_unloadingBranchNegDelv) ||
				(tsrStagePrevStepValue == tsrl_damMid_delun_neg_sigI_neg) || (tsrStagePrevStepValue == tsrl_damMid_delun_neg_sigI_pos)
				)
			{
				isLoading = 0;
				scalar_stress_t = scalar_stress_t_unloading_solution;
				DBCHK(dbout << "FL-1->L0,scalar_stress_t\t" << scalar_stress_t << '\n';);
			}
			else if ((tsrStagePrevStepValue == tsrl_damMid_loadingBranch) ||
					(tsrStagePrevStepValue == tsrl_dam0_sigI_less_sigC) || (tsrStagePrevStepValue == tsrl_dam0_sigI_over_sigC))
			{
				isLoading = 1;
				scalar_stress_t = scalar_stress_t_loading_solution;
				DBCHK(dbout << "FL-1->L1,scalar_stress_t\t" << scalar_stress_t << '\n';);
			}
			else
			{
				cout << "tsrStagePrevStepValue\t" << (int)tsrStagePrevStepValue << '\n';
				THROW("Invalid tsrStagePrevStepValue\n");
			}
		}
#if DiM1
		S[0] = scalar_stress_t;
#else
		S[0] = scalar_stress_t * del_u_nt_parts.vec_n_pos / eff_delu;
		double factor = scalar_stress_t / eff_delu * beta_delU * beta_delU;
		for (int i = 1; DiM; ++i)
			S[i] = factor * del_u_nt_parts.vec[i];
#endif
		DBCHK(dbout << "S[0]\n" << S[0] << '\n';);
		if (del_u_nt_parts.vec_n < 0) // want to prevent further penetration
		{
#if 0 // checked: never get here
			static double tol = 1e-3 * sigmaC;
			if (fabs(S[0]) > tol)
			{
				cout << S[0] << '\n';
				THROW("S[0] should not be nonzero\n");
			}
#endif
			DBCHK(dbout << "delunNeg\tval\t" << del_u_nt_parts.vec_n << '\n';);

			if (sigma_I_nt_parts.vec_n < 0.0) // delta_n is zero, sigmaI < 0, S[0] = 0, so this results in continued penetration. Want to prevent this
			{
				double delunk = del_u_nt_parts.vec_n * kCompressive;
				S[0] = MIN(sigma_I_nt_parts.vec_n, delunk);
				tsrStageLatestValue = tsrl_damMid_delun_neg_sigI_neg;
				DBCHK(dbout << "sINeg\tval\t" << sigma_I_nt_parts.vec_n << "\tdelunk\t" << delunk << "\tS\t" << S[0] << '\n';);
			}
			else
			{
				tsrStageLatestValue = tsrl_damMid_delun_neg_sigI_pos;
				DBCHK(dbout << "elseCs\ttsrStageLatestValue\t" << tsrStageLatestValue << '\n';);
			}
		}
		else if (isLoading == 1)
		{
			tsrStageLatestValue = tsrl_damMid_loadingBranch;
			DBCHK(dbout << "L1:tsrl_damMid_loadingBranch\n";);
		}
		else if (isLoading == 0)
		{
			bool isOpeningMode = ((sigma_I_nt_parts.vec_n - S[0]) > 0.0);
			DBCHK(dbout << "L0:isOpeningModeB" << isOpeningMode << "\n";);
			if (isOpeningMode)
				tsrStageLatestValue = tsrl_damMid_unloadingBranchPosDelv;
			else
				tsrStageLatestValue = tsrl_damMid_unloadingBranchNegDelv;
		}
		else
		{
			cout << "isLoading\t" << isLoading << '\n';
			THROW("isLoading is probably -1, should have been 0 or 1 by this time\n");
		}
		return;
	}
	if (tsrModel == tsr_linear)
	{
		THROW("Linear model not implemented\n");
	}
}


double Compute_Interface_EffectiveStress(const SL_Interface_Fracture_PF* interfacePFs, double sigmaC, Vec_nt& sigma_I_nt_parts, double& strengthNormalizer)
{
	SLFF_Eeffective_StsType effType = interfacePFs->effType;
	if (effType == effSs_Rankin)
		return sigma_I_nt_parts.vec_n_pos;
	strengthNormalizer = interfacePFs->gen_strength;
	if (effType == effSs_sqrt)
		return sigma_I_nt_parts.getEffective();
#if DiM2a3_F
	if (effType == effSs_Tresca)
		return sigma_I_nt_parts.vec_t_abs;
	if ((effType == effSs_mc) || (effType == effSs_mct))
	{
		double c = interfacePFs->c * sigmaC;

		double kscalar = interfacePFs->k;
		double sn = interfacePFs->sn;
		if (sn < 0)
		{
			double phi = atan(kscalar);
			//			double cphi = cos(phi), sphi = sin(phi);
			//			double snMC = 2.0 * cphi / (1.0 + sphi) * c;
			double snMC = c / kscalar;
			sn *= -snMC;
		}
		else
			sn *= sigmaC;

		double factor = sn / c;
		double effS = sigma_I_nt_parts.vec_t_abs;
		effS += kscalar * sigma_I_nt_parts.vec_n;
		if (effS < 0.0)
			effS = 0.0;

		if ((effType == effSs_mct) && (effS < sigma_I_nt_parts.vec_n_pos))
			effS = sigma_I_nt_parts.vec_n_pos;
		return  effS;
	}
#endif
	cout << "effType\t" << effType << '\n';
	THROW("invalid effective stress type\n");
}


//const VEC& uL, const VEC& uR, const VEC& vL, const VEC& vR, double maxEffDelU, double damage,
// 	const VEC& uL, const VEC& uR, const VEC& vL, const VEC& vR, double maxEffDelU, double damage,
//Vec_nt& del_u_nt_parts, Vec_nt& del_v_nt_parts, double& DDot, double& DTarget, double& DSrc_Allix, double & DSrc_RBH, VEC& sigmaIII)

onOffT Compute_Riemann_II_Slip_Aux(SLFF_SlipType slipTypeIn, Vec_nt& sigma_I_nt_parts, Vec_nt& del_u_nt_parts,
	const SL_Interface_Fracture_PF* interfacePFs, const SL_Elastic_InterfaceProperties*	ts_bulkProps,
	VEC& sigmaII, double& k, double& thetaTauFriction, double& delvIIAbs, double& thetaDelv, double& powerIDissipationSlip)
{
	// this is dissipation power from slip
	powerIDissipationSlip = 0.0;
	thetaTauFriction = 0.0;
	delvIIAbs = 0.0;
	thetaDelv = 0.0;
	setValue(sigmaII, 0.0);
	k = interfacePFs->damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax;
	if (slipTypeIn == sl_slip_off)
		return allOffV;

	//! 1. Computing normal and shear parts from mode I traction solution
	double abs_sIt_val = 0.0;
#if DiM2a3_F
	abs_sIt_val = sigma_I_nt_parts.vec_t_abs;
#endif

	if (abs_sIt_val < 1e-10)
		return allOffV;
	double sIn_val = sigma_I_nt_parts.vec_n;

	double damaged_sn, damaged_c;
	//! 2.a (slip traction solution) this is the case that it's easier to get slip solution
	if (ts_bulkProps->isIsoInterface)
	{
		// normal direction matches Riemann I
		sigmaII[0] = sigma_I_nt_parts.vec[0];

		double delV_Angle = 0.0;
		Compute_Directional_FrictionParameters(del_u_nt_parts, delV_Angle, interfacePFs, k, damaged_sn, damaged_c);
		if (sIn_val > damaged_sn)
			return allOffV;
		double sIIt_limit = damaged_c - k * sIn_val;
		if (sIIt_limit < 0)
			return allOffV;
		double diff = abs_sIt_val - sIIt_limit;
		if (diff < 0)
			return allOffV;
		delvIIAbs = 2.0 * ts_bulkProps->YAve_t * diff;
		powerIDissipationSlip = delvIIAbs * sIIt_limit;
#if DiM2
		if (sigma_I_nt_parts.vec_t[0] >= 0.0)
		{
			sigmaII[1] = sIIt_limit;
			thetaTauFriction = 0.0;
		}
		else
		{
			sigmaII[1] = -sIIt_limit;
			thetaTauFriction = PI;
		}
#endif
#if DiM3
		Vc_dm1 dir_sIIt;
		double norm_sIt = Normalize_m1(sigma_I_nt_parts.vec_t, dir_sIIt);
		thetaTauFriction = atan2(dir_sIIt[1], dir_sIIt[0]);
		for (int i = 0; i < DiMm1; ++i)
			sigmaII[i + 1] = sIIt_limit * dir_sIIt[i];
#endif
		thetaDelv = thetaTauFriction;
		return allOnV;
	}
#if !USE_ISO_ASSUMPTION
	// 2.b,c anisotropic cases
#if DiM2
		// eTheta is the direction of tauII (1, -1) in 2D
	double tauI = sigma_I_nt_parts.vec_t[0];
	bool pos_eTheta = (tauI >= 0);
	double tauITheta;

	double alphaThetaTheta = ts_bulkProps->alphaMat[0][0];
	double betaTheta;
	if (pos_eTheta)
	{
		thetaDelv = 0;
		tauITheta = tauI;
		betaTheta = ts_bulkProps->betaVec[0];
	}
	else
	{
		thetaDelv = PI;
		tauITheta = -tauI;
		betaTheta = -ts_bulkProps->betaVec[0];
	}
	thetaTauFriction = thetaDelv;
	Compute_Directional_FrictionParameters(del_u_nt_parts, thetaDelv, interfacePFs, k, damaged_sn, damaged_c);
	double denom = alphaThetaTheta - k * betaTheta;
	double num = tauITheta - damaged_c + k * sIn_val;
	delvIIAbs = num / denom;
	if (delvIIAbs <= 0)
		return allOffV;
	double tauIIhetaVal = tauITheta - alphaThetaTheta * delvIIAbs;
	if (tauIIhetaVal < 0)
		return allOffV;
	double sIIVal = sIn_val + betaTheta * delvIIAbs;
	if (sIIVal >= damaged_sn)
			return allOffV;
	powerIDissipationSlip = delvIIAbs * tauIIhetaVal;
	sigmaII[0] = sIIVal;
	if (pos_eTheta)
		sigmaII[1] = tauIIhetaVal;
	else
		sigmaII[1] = -tauIIhetaVal;
	return allOnV;
#endif
/// 3D
#if DiM3
	double theta0 = atan2(sigma_I_nt_parts.vec_t[1], sigma_I_nt_parts.vec_t[0]);
	Slip3D s3d;
	for (int i = 0; i < 2; ++i)
		s3d.tauItVal[i] = sigma_I_nt_parts.vec_t[i];
	s3d.sInval = sIn_val;
		
	int pos = 0;
	double del = PI / s3d.numPhiSegments;
	int numPts = s3d.numPhiSegments + 1;
	s3d.slp3ds.resize(numPts);
	s3d.slp3ds[pos].theta = theta0;
	SetSlip3DParameters(interfacePFs, ts_bulkProps, s3d, del_u_nt_parts, s3d.slp3ds[pos++]);

	double chnge = 0.0;
	int num = s3d.numPhiSegments / 2;
	for (int i = 0; i < num; ++i)
	{
		chnge += del;
		s3d.slp3ds[pos++].theta = theta0 + del;
		s3d.slp3ds[pos++].theta = theta0 - del;
	}
	s3d.indexMaxPower = -1;
	s3d.theta4MaxPower = 0;
	s3d.maximumPower = 0;
	Slip3D_1Dir * slip3d1dir;
	for (int i = 0; i < numPts; ++i)
	{
		slip3d1dir = &s3d.slp3ds[i];
		if (SetSlip3DParameters(interfacePFs, ts_bulkProps, s3d, del_u_nt_parts, *slip3d1dir))
		{
			if (slip3d1dir->power > s3d.maximumPower)
			{
				s3d.maximumPower = slip3d1dir->power;
				s3d.indexMaxPower = i;
				s3d.theta4MaxPower = slip3d1dir->theta;
			}
		}
	}
	if (s3d.indexMaxPower < 0)
		return allOffV;

	// this is the case that we have a valid slip solution
	Slip3D_1Dir* s3d1Dir = &s3d.slp3ds[s3d.indexMaxPower];
	powerIDissipationSlip = s3d.maximumPower;
	thetaDelv = s3d1Dir->theta;
	Compute_Directional_FrictionParameters(del_u_nt_parts, thetaDelv, interfacePFs, k, damaged_sn, damaged_c);
	sigmaII[0] = s3d1Dir->sIInVal;
	double v[2];
	for (int i = 0; i < 2; ++i)
		v[i] = s3d1Dir->n[i] * s3d1Dir->delvII;
	delvIIAbs = s3d1Dir->delvII;

	for (int i = 0; i < 2; ++i)
	{
		double tmp = s3d.tauItVal[i];
		for (int j = 0; j < 2; ++j)
			tmp -= ts_bulkProps->alphaMat[i][j] * v[j];
		sigmaII[i + 1] = tmp;
	}
	thetaTauFriction = atan2(sigmaII[2], sigmaII[1]);
	return allOnV;
#endif
#else
	THROW("For anisotropic material USE_ISO_ASSUMPTION_PRE should be 0\n");
#endif
}

Slip3D::Slip3D()
{
	numPhiSegments = 180;
	indexMaxPower = -1;
	theta4MaxPower = 0;
	maximumPower = 0.0;
	tauItVal[0] = 0.0;
	tauItVal[1] = 0.0;
	sInval = 0;
	delPhi = 1.0 / numPhiSegments;
}

ostream & operator<<(ostream & out, const Slip3D_1Dir & dat)
{
	out << "validSlip\t" << dat.validSlip << '\t';
	out << "power\t" << dat.power << '\t';
	out << "n\t" << dat.n[0] << '\t' << dat.n[1] << '\t';
	out << "theta\t" << dat.theta << '\t';
	out << "betaTheta\t" << dat.betaTheta << '\t';
	out << "alphaThetaTheta\t" << dat.alphaThetaTheta << '\t';
	out << "damaged_sn\t" << dat.damaged_sn << '\t';
	out << "damaged_c\t" << dat.damaged_c << '\t';
	out << "k\t" << dat.k << '\t';
	out << "delvII\t" << dat.delvII << '\t';
	out << "sIInVal\t" << dat.sIInVal << '\t';
	out << "tauIhetaVal\t" << dat.tauIhetaVal << '\t';
	out << "tauIIhetaVal\t" << dat.tauIIhetaVal << '\n';
	return out;
}

ostream & operator<<(ostream & out, const Slip3D & dat)
{
	out << "indexMaxPower\t" << dat.indexMaxPower << '\t';
	out << "maximumPower\t" << dat.maximumPower << '\t';
	out << "sInval\t" << dat.sInval << '\t';
	out << "numPhiSegments\t" << dat.numPhiSegments << '\t';
	out << "delPhi\t" << dat.delPhi << '\t';
	out << "tauItVal\t" << dat.tauItVal[0] << '\t' << dat.tauItVal[1] << '\n';
	int sz = dat.slp3ds.size();
	out << "slp3ds\n";
	for (int i = 0; i < sz; ++i)
		out << dat.slp3ds[i] << '\t';
	out << "\n";
	return out;
}

void Compute_Directional_FrictionParameters(Vec_nt& del_u_nt_parts, double delV_Angle, const SL_Interface_Fracture_PF* interfacePFs, double& k, double& damaged_sn, double& damaged_c)
{
	if (interfacePFs->slipOffOnMix == sl_slip_on)
	{
		k = 0.0;
		damaged_sn = 0.0;
		damaged_c = 0.0;
		return;
	}

	// these values may take nonzero values in future, but for now I set them to zero
	// < 0 means it's computed from c (c / k)
	damaged_sn = -1.0;
	damaged_c = 0.0;
	double k_upper, k_lower;
	Compute_Directional_k(delV_Angle, interfacePFs, k_upper, k_lower);

	// This is not fully correct, the tangential traction in the direction of slipis needed
	double directional_tangential_delU = 0.0;
#if DiM2a3_F
	if (interfacePFs->frictionModel != sl_friction_constant)
		directional_tangential_delU = del_u_nt_parts.vec_t_abs;
#endif

	bool dk_dDeluZero = true;

	GetDirectionalFrictionCoefficientFromRaw_k(interfacePFs, k_upper, k_lower, directional_tangential_delU, k);

	if (k > 0)
	{
		if (damaged_sn < 0.0)
			damaged_sn = damaged_c / k;
		else
			damaged_sn = MIN(damaged_sn, damaged_c / k);
	}
	else
	{
		damaged_c = 0.0;
		damaged_sn = 0.0;
	}
}

void Compute_Directional_k(double delV_Angle, const SL_Interface_Fracture_PF* interfacePFs, double& k_upper, double& k_lower)
{
	if (interfacePFs->isIsoFrictionParameters)
	{
		k_upper = interfacePFs->damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax;
		k_lower = k_upper;
		if (interfacePFs->frictionModel != sl_friction_constant)
			k_lower = interfacePFs->damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax;
	}
	else
	{
#if DiM2
		if (fabs(delV_Angle) < 1e-2)
		{
			k_upper = interfacePFs->damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax;
			k_lower = k_upper;
			if (interfacePFs->frictionModel != sl_friction_constant)
				k_lower = interfacePFs->damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax;
		}
		else
		{
			k_upper = interfacePFs->damaged_kUpper_2D_neg_y_dir_or_3D_angularMin;
			k_lower = k_upper;
			if (interfacePFs->frictionModel != sl_friction_constant)
				k_lower = interfacePFs->damaged_kLower_2D_neg_y_dir_or_3D_angularMin;
		}
		return;
#endif
#if DiM3
		double k_upperOverAngle = interfacePFs->damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax, k_lowerOverAngle = interfacePFs->damaged_kUpper_2D_neg_y_dir_or_3D_angularMin;
		double A = 0.5 * (k_upperOverAngle + k_lowerOverAngle);
		double B = A - k_lowerOverAngle;
		double theta0 = interfacePFs->theta4MaxValRad;
		k_upper = A + B * cos(delV_Angle - theta0);
		k_lower = k_upper;
		if (interfacePFs->frictionModel != sl_friction_constant)
		{
			k_lowerOverAngle = interfacePFs->damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax, k_lowerOverAngle = interfacePFs->damaged_kLower_2D_neg_y_dir_or_3D_angularMin;
			A = 0.5 * (k_lowerOverAngle + k_lowerOverAngle);
			B = A - k_lowerOverAngle;
			k_lower = A + B * cos(delV_Angle - theta0);
		}
#endif
	}
}

void GetDirectionalFrictionCoefficientFromRaw_k(const SL_Interface_Fracture_PF* interfacePFs, double k_upper, double k_lower, double directional_tangential_delU, double& k)
{
	if (interfacePFs->frictionModel == sl_friction_constant)
	{
		k = k_upper;
		return;
	}
	else if (interfacePFs->frictionModel == sl_friction_slip_weakeningLinear)
	{
		double e = interfacePFs->damaged_para0_epsSW * interfacePFs->damaged_slip_dt;
		if (interfacePFs->damaged_para0_epsSW < 0 || interfacePFs->damaged_para0_epsSW > 0.5)
			THROW("regularization of transition in expressed as a ratio of Dc, it needs to be positive (lenghts cannot be negative). It also needs to be less than 0.5 otherwise we regularization regions overlap in the linear section!");

		if (interfacePFs->damaged_slip_dt <= 0)
			THROW("Dc has to be positive!");

		double W_slip = (k_lower - k_upper) / interfacePFs->damaged_slip_dt;
		double slipDist = directional_tangential_delU;

		double dk_d_directional_delu;
		if (slipDist < e)
			C1_interpolateC(k_upper, 0, 0, k_upper + W_slip * e, W_slip, e, directional_tangential_delU, k, dk_d_directional_delu);
		else if ((slipDist >= e) && (slipDist <= interfacePFs->damaged_slip_dt - e))
		{
			// linear section
			k = k_upper + W_slip * slipDist;
			dk_d_directional_delu = W_slip;
		}
		else if ((slipDist > interfacePFs->damaged_slip_dt - e) && (slipDist < interfacePFs->damaged_slip_dt + e))
		{
			// 2nd regularization section
			C1_interpolateC(k_upper + W_slip * (interfacePFs->damaged_slip_dt - e), W_slip, interfacePFs->damaged_slip_dt - e, k_lower, 0, interfacePFs->damaged_slip_dt + e, slipDist, k, dk_d_directional_delu);
		}
		else if (slipDist >= interfacePFs->damaged_slip_dt + e)
		{
			// constant section
			k = k_lower;
			dk_d_directional_delu = 0.0;
		}
		else
		{
			THROW("Invalid slip distance");
		}
	}
	else if (interfacePFs->frictionModel == sl_friction_slip_weakeningPower)
	{
		// law is of the form f_sw = k_0 + A*exp(-n * d) // where k_0, A, and n are constants, and d is the slip along the direction of the face
		double f_0 = k_lower; // cohIProp->friction_model_first_parameter();
		double A = k_upper - k_lower;
		double n = interfacePFs->damaged_para1_n;
		double xn, dxn_dx;
		xn = n * directional_tangential_delU / interfacePFs->damaged_slip_dt;
		dxn_dx = n / interfacePFs->damaged_slip_dt;

		double tmp = A * exp(-xn);
		k = f_0 + tmp;
	}
	else
	{
		cout << "frictionModel\n" << interfacePFs->frictionModel << '\n';
		THROW("Invalid friction model type");
	}
}

// ian mcnamara added
void C1_interpolateC(const double p0, const double m0, const double x0, const double p1, const double m1, const double x1, const double x, double & value, double & slope)
{
	assert(x1 > x0);
	assert(x >= x0);
	assert(x <= x1);

	double t = (x - x0) / (x1 - x0);
	double DtDx = 1 / (x1 - x0);
	double p0_f = 2 * pow(t, 3) - 3 * pow(t, 2) + 1;
	double m0_f = pow(t, 3) - 2 * pow(t, 2) + t;
	double p1_f = -2 * pow(t, 3) + 3 * pow(t, 2);
	double m1_f = pow(t, 3) - pow(t, 2);
	value = p0 * p0_f + m0 * (x1 - x0)*m0_f + p1 * p1_f + m1 * (x1 - x0)*m1_f;

	double Dp0_fDt = 6 * pow(t, 2) - 6 * t;
	double Dm0_fDt = 3 * pow(t, 2) - 4 * t + 1;
	double Dp1_fDt = -6 * pow(t, 2) + 6 * t;
	double Dm1_fDt = 3 * pow(t, 2) - 2 * t;
	slope = p0 * Dp0_fDt*DtDx + m0 * (x1 - x0)*Dm0_fDt*DtDx + p1 * Dp1_fDt*DtDx + m1 * (x1 - x0)*Dm1_fDt*DtDx;


}

bool SetSlip3DParameters(const SL_Interface_Fracture_PF* interfacePFs, const SL_Elastic_InterfaceProperties* ts_bulkProps, Slip3D& s3d, Vec_nt& del_u_nt_parts, Slip3D_1Dir& s3d1Dir)
{
#if !USE_ISO_ASSUMPTION
	double theta = s3d1Dir.theta;
	s3d1Dir.n[0] = cos(theta), s3d1Dir.n[1] = sin(theta);

	Compute_Directional_FrictionParameters(del_u_nt_parts, theta, interfacePFs, s3d1Dir.k, s3d1Dir.damaged_sn, s3d1Dir.damaged_c);
	s3d1Dir.tauIhetaVal = 0.0;
	s3d1Dir.betaTheta = 0.0;
	for (int i = 0; i < 2; ++i)
	{
		s3d1Dir.tauIhetaVal += s3d1Dir.n[i] * s3d.tauItVal[i];
		s3d1Dir.betaTheta += s3d1Dir.n[i] * ts_bulkProps->betaVec[i];
	}
	s3d1Dir.alphaThetaTheta = 0.0;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			s3d1Dir.alphaThetaTheta += s3d1Dir.n[i] * s3d1Dir.n[j] * ts_bulkProps->alphaMat[i][j];

	double kVal = s3d1Dir.k, sInVal = s3d.sInval;
	double denom = s3d1Dir.alphaThetaTheta - kVal * s3d1Dir.betaTheta;
	double num = s3d1Dir.tauIhetaVal - s3d1Dir.damaged_sn + kVal * sInVal;
	s3d1Dir.delvII = num / denom;
	if (s3d1Dir.delvII <= 0)
	{
		s3d1Dir.validSlip = false;
		return s3d1Dir.validSlip;
	}
	s3d1Dir.tauIIhetaVal = s3d1Dir.tauIhetaVal - s3d1Dir.alphaThetaTheta * s3d1Dir.delvII;
	if (s3d1Dir.tauIIhetaVal < 0)
	{
		s3d1Dir.validSlip = false;
		return s3d1Dir.validSlip;
	}
	s3d1Dir.sIInVal = sInVal + s3d1Dir.betaTheta * s3d1Dir.delvII;
	if (s3d1Dir.sIInVal >= s3d1Dir.damaged_sn)
	{
		s3d1Dir.validSlip = false;
		return s3d1Dir.validSlip;
	}
	s3d1Dir.power = s3d1Dir.delvII * s3d1Dir.tauIIhetaVal;
	s3d1Dir.validSlip = true;
	return s3d1Dir.validSlip;
#else
	THROW("Function should only be called for anisotropic case\n");
#endif
}

void SLInterfaceCalculator::Set_tPtSlns(SL_interface_Temp_PPtData* tPtSlnsIn, bool tPtSlnsDeletableIn)
{
	tPtSlns = tPtSlnsIn;
	tPtSlns_Deletable = tPtSlnsDeletableIn;
}

void SLInterfaceCalculator::Set_pPtSlns(SL_interfacePPtData* pPtSlnsIn, bool pPtSlnsDeletableIn)
{
	pPtSlns = pPtSlnsIn;
	pPtSlns_Deletable = pPtSlnsDeletableIn;
}

Periodic1IntrFrag_TimeStageStorage::Periodic1IntrFrag_TimeStageStorage()
{
	b_set = 0;
	pitss = pit_sigma0;
	timeIndex = 0;
	perConf = NULL;
}

void Periodic1IntrFrag_TimeStageStorage::CopyData(const Periodic1IntrFrag_TimeStageStorage& other)
{
	b_set = other.b_set;
	pitss = other.pitss;
	perConf = other.perConf;
	timeIndex = other.timeIndex;
	vals = other.vals;
	spatial_stress_vals = other.spatial_stress_vals;
	spatial_vel_vals = other.spatial_vel_vals;
}

void Periodic1IntrFrag_TimeStageStorage::Initialize(PITSS pitssIn, const string& nameBase, Periodic1IntrFrag* perConfIn)
{
	pitss = pitssIn;
	perConf = perConfIn;
	string name = nameBase + "_criterion_", ser;
	unsigned int ipitss = (int)pitss;
	toString(ipitss, ser);
	name += ser;
	name += "_";
	name += getName(pitss);
	name += ".txt";
	fstream in(name.c_str(), ios::in);
	if (!in.is_open())
	{
		out.open(name.c_str(), ios::app);
		out << "success\ta\tl\tlog10(a)\tlog10(l)\ttimeIndex";
		for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
			out << "\t" << getName((OneSegmentPFT)i);
		for (unsigned int i = 0; i < perConf->stressHeader.size(); ++i)
			out << "\t" << perConf->stressHeader[i];
		for (unsigned int i = 0; i < perConf->velHeader.size(); ++i)
			out << "\t" << perConf->velHeader[i];
		out << '\n';
	}
	else
	{
		in.close();
		out.open(name.c_str(), ios::app);
	}
}

void Periodic1IntrFrag_TimeStageStorage::Print()
{
	out << b_set << '\t';
	out << perConf->a << '\t';
	out << perConf->l << '\t';
	out << perConf->log10a << '\t';
	out << perConf->log10l << '\t';
	if (b_set != 0)
	{
		out << timeIndex;
		for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
			out << "\t" << vals[i];
		for (unsigned int i = 0; i < perConf->velHeader.size(); ++i)
			out << "\t" << spatial_stress_vals[i];
		for (unsigned int i = 0; i < perConf->stressHeader.size(); ++i)
			out << "\t" << spatial_vel_vals[i];
		out << '\n';
	}
	else
	{
		string nan = "nan";
		out << -1;
		for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
			out << "\t" << nan;
		for (unsigned int i = 0; i < perConf->stressHeader.size(); ++i)
			out << "\t" << nan;
		for (unsigned int i = 0; i < perConf->velHeader.size(); ++i)
			out << "\t" << nan;
		out << '\n';
	}
	out.close();
}

bool GetZeroCrosserValue(double tF, double tOFC, double vF, double vtOFCbar, double& tF_minus_tOFC, double& crosserValue, bool OFC_v_decreasing)
{
	if (tOFC < 0) // this time has not reached
	{
		tF_minus_tOFC = std::numeric_limits<double>::quiet_NaN();
		crosserValue = vtOFCbar - vF;
		if (!OFC_v_decreasing)
			crosserValue *= -1;
		return false;
	}
	tF_minus_tOFC = tF - tOFC;
	crosserValue = tF_minus_tOFC;
	return true;
}

Periodic1IntrFrag::Periodic1IntrFrag()
{
	setEmpty_Periodic1IntrFrag();
}

void Periodic1IntrFrag::setEmpty_Periodic1IntrFrag()
{
	max_bsigma = -0.1 * DBL_MAX;
	max_bsigma_pw = max_bsigma;
	aIndex = 0;
	lIndex = 0;
	aIndex_secondary = 0;
	lIndex_secondary = 0;
	a = 10.0;
	l = 2.0;
	relTol = 1e-4;
	energy_diss_per_length_At_t_dilute_set = false;
	numSpatialSubsegments_BulkInterfacePoints_Print_4PP = 10;
	useRepeatedSimpsonRuleForHigherOrders_4PP = true;
	delvsZeroObserved = false;
	nameBase = "periodic1IntrFrag";
	step4_segment_vsigma_output = 1;
	if (g_1_interface_low_disk_space != 0)
		step4_segment_vsigma_output = 10000;
	terminateState = PITSS_none;
	isExtrinsic = true;
	energyScale = 0.5;
}

void Periodic1IntrFrag::Set_TSR_Model(SLFF_TSRType	tsrModelIn)
{
	isExtrinsic = IsExtrinsic(tsrModelIn);
	energyScale = GetEnergyConstantFactor(tsrModelIn);
}

void Periodic1IntrFrag::Initialize_Periodic1IntrFrag(bool prematureExit)
{
//	setEmpty_Periodic1IntrFrag();
	fstream inla("l10a.txt", ios::in);
	if (inla.is_open())
	{
		double tmp;
		inla >> tmp; // is log10 of a
		a = pow(10.0, tmp);
		int tmpi;
		inla >> tmpi;
		if (!inla.eof())
			aIndex = tmpi;
		inla >> tmpi;
		if (!inla.eof())
			aIndex_secondary = tmpi;
	}
	fstream inll("l10l.txt", ios::in);
	if (inll.is_open())
	{
		double tmp;
		inll >> tmp; // is log10 of l
		l = pow(10.0, tmp);
		int tmpi;
		inll >> tmpi;
		if (!inll.eof())
			lIndex = tmpi;
		inll >> tmpi;
		if (!inll.eof())
			lIndex_secondary = tmpi;
	}

	la = l * a;
	aInv = 1.0 / a;
	if (!isExtrinsic)
		t0 = 0.0;
	else
		t0 = 1.0 / a;

	string aIndex_s, lIndex_s;
	toString(aIndex, aIndex_s);
	toString(lIndex, lIndex_s);
	nameOne_a_SharedAll_l = nameBase + "_aI_" + aIndex_s;
	string ser;
	toString(aIndex_secondary, ser);
	nameOne_a_SharedAll_l += "_";
	nameOne_a_SharedAll_l += ser;

	nameOne_a_One_l = nameOne_a_SharedAll_l + "_lI_" + lIndex_s;
	toString(lIndex_secondary, ser);
	nameOne_a_One_l += "_";
	nameOne_a_One_l += ser;

	// (x): equation number x from Zhu06
	log10a = rlog10(a);
	log10l = rlog10(l);
	zhu6_epsilonDotScale = 1.0 / energyScale;//(13)
	a_p = a / zhu6_epsilonDotScale; // (11) normalized loadingrate from Zhu06
	zhu6_sBarScale = energyScale; // (21)
	l_p = l / zhu6_sBarScale;

	l_Zhu6a = zhu6_sBarScale * 4.5 / (1.0 + 6.00 * pow(a_p, 2.0 / 3.0)); //(28)
	l_Zhu6b = zhu6_sBarScale * 4.5 / (1.0 + 0.77 * pow(a_p, 0.25) + 5.4 * pow(a_p, 0.75)); // (29a)
	l_Grady = zhu6_sBarScale * pow(24.0 / a_p / a_p, 1.0 / 3.0); // (1')
	l_Glenn = zhu6_sBarScale * 4.0 / a_p * sinh(1.0 / 3 * asinh(1.5 * a_p));// (2')

#if 0 // older calculations
	double l_Grady_b = pow(24.0 * energyScale / a / a, 1.0 / 3.0);
	double alpha = 3.0 / a / a;
	double beta = alpha * energyScale / 2.0; // 3.0 G / 2.0 / a / a;
	double phi = asinh(beta * pow(3.0 / alpha, 1.5));
	double l_Glenn_b = 4.0 * sinh(phi / 3.0) * sqrt(alpha / 3.0);
#endif

	// reference wherein t0 = Z deltaC/2/sigmaC (Zhu06, (12)), and s0 is different ...
	tFailure = -1.0;
	if (isExtrinsic) // Ortiz model - Zhu06, (16)
	{
		// equation is Zhu06, (16), a_p(exp(t_p) - t_p - 1) = 1.0
		if (a_p > 0.499)
			t_dilute_approx_p = sqrt(2.0 / a_p);
		else
			t_dilute_approx_p = log(1.0 / a_p);

		double tol_t = t_dilute_approx_p * 1e-9;
		unsigned int cntr = 0;
		double tmp, fl, fpl, delt;
		t_dilute_p = t_dilute_approx_p;
		while (cntr++ < 1000)
		{
			tmp = exp(t_dilute_p) - 1.0;
			fl = a_p * (tmp - t_dilute_p) - 1.0;
			fpl = a_p * tmp;
			delt = -computeRatio(fl, fpl);
			t_dilute_p += delt;
			//		cout << cntr << '\t' << t_dilute_p << '\t' << fl << '\t' << delt << '\n';
			if ((delt < tol_t) && (fl < 1e-9))
				break;
		}
		if (cntr == 1000)
		{
			THROW("Cannot find the approximate solution. Can comment out throw and get to the next line.\n");
			t_dilute_p = t_dilute_approx_p;
		}
		delu4Sigma0 = 1.0;
		t_SigmaMax_dilute = aInv;
		t_SigmaMax_real = aInv;

		l_dilute_p = t_dilute_p;
		l_dilute_approx_p = t_dilute_approx_p;

		// getting back to real fragment size, using sBar (Zhu06 eqn (21)), zhu6_sBarScale = GFactor = (0.5 for Ortiz, e for Xu-Needleman)
		l_dilute_approx = zhu6_sBarScale * l_dilute_approx_p;
		l_dilute = zhu6_sBarScale * l_dilute_p;

		t_dilute_zeroStressCheck = l_dilute + aInv - t0;

		t_dilute = energyScale * t_dilute_p;
		t_dilute_approx = energyScale * t_dilute_approx_p;
		if (diluteFractureModel)
			tFailure = l_dilute + aInv;

		/// now incrementing t_dilutes to match start from time zero
		t_dilute += aInv;
		t_dilute_approx += aInv;

		t_dilute_p += aInv / energyScale;
		t_dilute_approx_p += aInv / energyScale;

		// f(delta) + 0.5 deltaDot = at
		// vsolid = v(x = L/2) - v(x = 0) = al/2 - deltaDot/2 = 0
		// -> l = 2t - 2 f(delta) / a, where t, delta are taken for a reference point, and t is total time from time zero
		// A) initial stage of unloading:
		// t = 1/a, f(delta) = 1 -> l = 2/a - 2/a = 0 ! Not appropriate! This is similar to what Drugan does for XuNeedleman at max load but it does not apply here
		// B) final stage  (delta = 1, sigma = 0)
		// -> l = 2 t_final - 2 * 0 / a = 2 * t final
		l_SigmaRef_DruganFinal_Dilute = 2.0 * t_dilute;
		l_SigmaRef_DruganMax_Dilute = std::numeric_limits<double>::quiet_NaN();
		if (prematureExit)
			return;
	}
	else // intrinsic Xu-Needleman, Drugan, eqn (36) and also (31)
	{
		tsr_xn.a = a;
		tsr_xn.Z = 1.0;
		tsr_xn.sigmaC = 1.0;
		tsr_xn.deltaC = 1.0;
		tsr_xn.sigmaCFactor4Zero = 0.001; // 0.01
		tsr_xn.delTFactor = 0.0001; // 0.00001;
		string name = nameOne_a_One_l + "_TSR_XuNeedleman.txt";
		fstream out(name.c_str(), ios::out);
		unsigned int step = 1;
		if (a < 1)
			step = 100;
		tsr_xn.Compute(&out, 100);

		t_dilute = tsr_xn.tSigmaZero;
		t_dilute_approx = t_dilute;
		
		l_dilute = t_dilute;
		l_dilute_approx = t_dilute_approx;

		// other normalization
		t_dilute_p = t_dilute / energyScale;
		t_dilute_approx_p = t_dilute_approx / energyScale;

		l_dilute_p = l_dilute / zhu6_sBarScale;
		l_dilute_approx_p = l_dilute_approx / zhu6_sBarScale;

		/////////////////////////////////////////////////
		t_dilute_zeroStressCheck = t_dilute; // this is the time that for dilute model stress is "0" e.g. around 0.01 stress
		delu4Sigma0 = tsr_xn.delnSigmaZero;
		
		////////////////////////////////////////////////
		// Drugan model
		t_SigmaMax_dilute = tsr_xn.tSigmaMax;
		// getting segment size from Drugan idea:
		// vSolid = al/2 - deltaDot/2 = 0 -> deltDot = al
		// -> deltaDot = al
		// from characteristics
		// sigma(delta) + Z/2 deltaDot = at,
		// at tR, delta = 1 -> sigma = 1; Moreover we have Z = 1.0 ->
		// 1 + 0.5 deltaDot = at -> 
		// 1 + 0.5 al = a tR -> l = 2(a tR - 1)/a = 2 tR - 2/a // Drugan equation (27)
		l_SigmaRef_DruganMax_Dilute = 2.0 * t_SigmaMax_dilute - 2.0 / a;
		l_SigmaRef_DruganFinal_Dilute = 2.0 * t_dilute;

		if (diluteFractureModel)
			tFailure = t_dilute;

		if (prematureExit)
			return;

#if 0
		// Drugan (31), (36) a_p(e(-t_p) + t_p - 1) = 1  NOTE: t' is measured fro mtime zero NOT time of sigmaC
		// Eventually time to sigmaC is subtracted
		if (a_p > 0.499)
			t_dilute_approx_p = sqrt(2.0 / a_p);
		else
			t_dilute_approx_p = 1.0 / a_p;

		double tol_t = t_dilute_approx_p * 1e-9;
		unsigned int cntr = 0;
		double tmp, fl, fpl, delt;
		t_dilute_p = t_dilute_approx_p;
		while (cntr++ < 1000)
		{
			tmp = exp(-t_dilute_p) - 1.0;
			fl = a_p * (tmp + t_dilute_p) - 1.0;
			fpl = -a_p * tmp;
			delt = -computeRatio(fl, fpl);
			t_dilute_p += delt;
			//		cout << cntr << '\t' << t_dilute_p << '\t' << fl << '\t' << delt << '\n';
			if ((delt < tol_t) && (fl < 1e-9))
				break;
		}
		if (cntr == 1000)
		{
			THROW("Cannot find the approximate solution. Can comment out throw and get to the next line.\n");
			t_dilute_p = t_dilute_approx_p;
		}
#endif
	}
	t_SigmaZero_minus_SigmaMax_dilute = t_dilute - t_SigmaMax_dilute;

	diluteFractureModel = (l >= 0.9999 * l_dilute);
	//t0 = aInv;
	double suggested_delt_a = aInv;
//	double suggested_delt_l = 1.0 / l;	//	double suggested_delt_dilute_l = 1.0 / l_dilute;
//	double del_t = MIN(MIN(suggested_delt_a, suggested_delt_l), suggested_delt_dilute_l);
	double del_t = MIN(suggested_delt_a, 1.0);
	del_t = MIN(del_t, l);
	suggeste_delt = relTol * del_t;
	double fact = 10.0;
	if (!isExtrinsic)
		fact = 20.0;
	if (diluteFractureModel)
		suggested_final_time4ZeroSigmaC = t_dilute_zeroStressCheck;
	else
		suggested_final_time4ZeroSigmaC = fact * t_dilute_zeroStressCheck * l_dilute / l; // 10.0 * l_dilute * l_dilute / l;

	energy_diss_per_length = energyScale / l;
	energy_diss_per_length_At_t_dilute = energyScale / l;

	//	sigmaMaxPrecalculated = 1.0 + 0.5 * la;
//	if (la > )
	SetNewtonCotes_Points_AndWeights(numSpatialSubsegments_BulkInterfacePoints_Print_4PP, useRepeatedSimpsonRuleForHigherOrders_4PP, weights, xs_wrCenter);
	sz_xs = xs_wrCenter.size();
	velVals.resize(sz_xs);
	stressVals.resize(sz_xs);

	half_l = 0.5 * l;
	for (unsigned int i = 0; i < sz_xs; ++i)
	{
		xs_wrCenter[i] *= half_l;
//		weights[i] *= half_l; // not scaled by size as we are calculating averages anyway (if half_l is used here, a division by l/2 is needed later)
	}
	currentStepVals.resize(OneSegmentPFT_SIZE);
	fill(currentStepVals.begin(), currentStepVals.end(), 0.0);
	currentStepVals[pft_ivsolid] = 0.5 * la;
	if (t0 >= 1e-9)
	{
		currentStepVals[pft_time] = t0;
		currentStepVals[pft_isigma] = 1.0;
		currentStepVals[pft_ichar] = 1.0;
		currentStepVals[pft_irelusolid] = 1.0; //0.5 l / 0.5 l
		currentStepVals[pft_bepsilon] = 1.0;
		currentStepVals[pft_bsigma] = 1.0;
		currentStepVals[pft_bsigma_maxpw] = 1.0;
		currentStepVals[pft_bU_PL] = 0.5;
		double Ut0 = 0.5 * sigma_t0 * sigma_t0;
		currentStepVals[pft_bEneSource_PL] = Ut0;
	}
	currentStepVals[pft_ilog10Ened] = rlog10(currentStepVals[pft_iEned]);
	currentStepVals[pft_ilog10Ened_PL] = rlog10(currentStepVals[pft_iEned_PL]);
	stat4Vals.resize(OneSegmentPFT_SIZE);
	string name = nameOne_a_One_l + "_outAllVals.txt";
	outAllVals.open(name.c_str(), ios::out);
	name = nameOne_a_One_l + "_out_Stress.txt";
	outStress.open(name.c_str(), ios::out);
	name = nameOne_a_One_l + "_out_Vel.txt";
	outVel.open(name.c_str(), ios::out);
	for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
	{
		string name = getName((OneSegmentPFT)i);
		stat4Vals[i].setName(name, name);
	}
	velHeader.clear();
	stressHeader.clear();
	outStress << "timeIndex\ttime\tnumCycles";
	outVel << "timeIndex\ttime\tnumCycles";
	for (unsigned int i = 0; i < sz_xs; ++i)
	{
		string name_x = "x_i_", tmp;
		toString(i, tmp);
		name_x += tmp;
		name_x += "_xr2c_";
		double tmpd = xs_wrCenter[i] / half_l;
		toString(tmpd, tmp);
		name_x += tmp;
		outStress << '\t' << name_x;
		outVel << '\t' << name_x;
		velHeader.push_back("v_" + name_x);
		stressHeader.push_back("s_" + name_x);
	}
	outStress << '\n';
	outVel << '\n';
	outStress << "0\t" << t0 << "\t0.0";
	outVel << "0\t" << t0 << "\t0.0";
	for (unsigned int i = 0; i < sz_xs; ++i)
	{
		outStress << '\t' << sigma_t0;
		outVel << '\t' << 0.0;
	}
	outStress << '\n';
	outVel << '\n';

	outAllVals << "timeIndex";
	for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
		outAllVals << '\t' << getName((OneSegmentPFT)i);
	outAllVals << '\n';
	timeIndexNew = 0,
	outAllVals << timeIndexNew++;
	for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
		outAllVals << '\t' << currentStepVals[i];
	outAllVals << '\n';

	for (unsigned int i = 0; i < PITSS_SIZE; ++i)
		stageSlns[i].Initialize((PITSS)i , nameOne_a_SharedAll_l, this);
}

PITSS  Periodic1IntrFrag::UpdateStats(double timeNew, double delu, double delv, double sigma)
{
	double timeNewAbsolute = timeNew + t0;
	double timeRelSigmaMax = timeNewAbsolute;
	if (isExtrinsic)
		timeRelSigmaMax -= aInv;
	vector<PITSS> pits;
	bool printSegmentvsigma = (timeIndexNew % step4_segment_vsigma_output == 0);
	if (g_1_interface_low_disk_space == 1)
		printSegmentvsigma = false;
	double halfla = 0.5 * la;
	double prevTime = currentStepVals[pft_time];
	double delT = timeNewAbsolute - prevTime;
	currentStepVals[pft_time] = timeNewAbsolute;
	double numCycles = timeRelSigmaMax / l;
	currentStepVals[pft_numCyclesAfterCrackOpening] = numCycles;
	double prevDelu = currentStepVals[pft_idelu];
	currentStepVals[pft_idelu] = delu;
	double prevDelv = currentStepVals[pft_idelv];
	currentStepVals[pft_idelv] = delv;
	double prevSigma = currentStepVals[pft_isigma];
	currentStepVals[pft_isigma] = sigma;
	currentStepVals[pft_ichar] = sigma + 0.5 * delv;
	currentStepVals[pft_iEne] += 0.5 * (sigma + prevSigma) * (delu - prevDelu);
	currentStepVals[pft_iEner] = 0.5 * (sigma * delu);
	currentStepVals[pft_iEned] = currentStepVals[pft_iEne] - currentStepVals[pft_iEner];
	currentStepVals[pft_ilog10Ened] = rlog10(currentStepVals[pft_iEned]);
	double delvs = halfla - 0.5 * delv;
	currentStepVals[pft_ivsolid] = delvs;
	double reldus = (halfla * timeNewAbsolute - 0.5 * delu) / half_l;
	currentStepVals[pft_irelusolid] = reldus;
	bool current_delvsZero = false;
	if (delvsZeroObserved == false)
	{
		if (delv < 1e-13 * la)
		{
			current_delvsZero = true;
			delvsZeroObserved = true;
		}
	}
	currentStepVals[pft_bepsilon] = a * timeNewAbsolute;
	double bsigma = 0.0, bsigma_maxpw = sigma, K = 0.0, U = 0.0;

	if (printSegmentvsigma)
	{
		outStress << timeIndexNew << '\t' << timeNewAbsolute << '\t' << numCycles;
		outVel << timeIndexNew << '\t' << timeNewAbsolute << '\t' << numCycles;
	}
	for (unsigned int i = 0; i < sz_xs; ++i)
	{
		double x = xs_wrCenter[i];
		double weight = weights[i];
		bool lastPoint = (fabs(x - half_l) < 0.0001 * half_l);
		// x is w.r.t. the center of the segment
		// A is the point on the left  side of point (x, timeNew) -> left side, right-going characteristics
		// B is the point on the right side of point (x, timeNew) -> right side, left-going characteristics
		double sigma_pw, vTheta_pw;
		if (!lastPoint)
		{
			double sigmaA, vA, sigmaB, vB, vL, vR;
			double tA = timeNew - (x + half_l); // /c but wave speed = 1
			if (tA < 0)
			{
				sigmaA = sigma_t0;
				vA = a * (x - timeNew); // timeNew * c but c = 1.0
			}
			else
			{
				bool found = g_seq_short->GetPt(tA, vL, vR, sigmaA);
				if (!found)
				{
					THROW("Could not find the point!\n");
				}
				vA = vR - halfla; // -0.5 la is to take care of -0.5l shift to the left relative to the center of the segment
			}

			double tB = timeNew - (half_l - x); // /c but wave speed = 1
			if (tB < 0)
			{
				sigmaB = sigma_t0;
				vB = a * (x + timeNew); // timeNew * c but c = 1.0
			}
			else
			{
				bool found = g_seq_short->GetPt(tB, vL, vR, sigmaB);
				if (!found)
				{
					double tm;
					cout << "FirstPt\t" << g_seq_short->GetPt(0, tm, vL, vR, sigma);
					cout << "\time\t" << tm << "\tvL\t" << vL << "\tvR\t" << vR << "\tsigma\t" << sigma << '\n';
					cout << "LastPt\t" << g_seq_short->GetPt(-1, tm, vL, vR, sigma);
					cout << "\time\t" << tm << "\tvL\t" << vL << "\tvR\t" << vR << "\tsigma\t" << sigma << '\n';
					cout << "tB\t" << tB << '\n';
					if (per_if != NULL)
						cout << "a\t" << per_if->a << "\tl\t" << per_if->l << '\n';
					THROW("Could not find the point!\n");
				}
				vB = vL + halfla; // +0.5 la is to take care of +0.5l shift to the right relative to the center of the segment
			}
			double wL_rightGoing = sigmaA - vA;
			double wR_leftGoing = sigmaB + vB;
			sigma_pw = 0.5 * (wR_leftGoing + wL_rightGoing);
			vTheta_pw = 0.5 * (wR_leftGoing - wL_rightGoing) - a * x; // -ax is to get to vTheta from v'Theta, see Zhou_2006_Molinari_Ramesh_Analysis of the brittle fragmentation of an expanding ring.pdf eqn 9
		}
		else
		{
			sigma_pw = sigma;
			vTheta_pw = -0.5 * delv;
		}
		stressVals[i] = sigma_pw;
		velVals[i] = vTheta_pw;

		bsigma += weight * sigma_pw;
		K += weight * vTheta_pw * vTheta_pw;
		U += weight * sigma_pw * sigma_pw;
		if (sigma_pw > bsigma_maxpw)
		{
			bsigma_maxpw = sigma_pw;
		}
		if (printSegmentvsigma)
		{
			outStress << '\t' << sigma_pw;
			outVel << '\t' << vTheta_pw;
		}
	}
	if (bsigma > max_bsigma)
	{
		max_bsigma = bsigma;
		pits.push_back(pit_timeSigmaAveMax);
	}
	if (bsigma_maxpw > max_bsigma_pw)
	{
		max_bsigma_pw = bsigma_maxpw;
		pits.push_back(pit_timeSigmaPWMax);
	}
	if (printSegmentvsigma)
	{
		outStress << '\n';
		outVel << '\n';
	}
	K *= 0.5;
	U *= 0.5;
	double bsigmaPrev = currentStepVals[pft_bsigma];
	currentStepVals[pft_bsigma] = bsigma;
	currentStepVals[pft_bsigma_maxpw] = bsigma_maxpw;
	currentStepVals[pft_bK_PL] = K;
	currentStepVals[pft_bU_PL] = U;
	currentStepVals[pft_bEneSource_PL] += 0.5 * a * delT * (bsigma + bsigmaPrev);
	currentStepVals[pft_iEned_PL] = currentStepVals[pft_iEned] / l;
	currentStepVals[pft_ilog10Ened_PL] = rlog10(currentStepVals[pft_iEned_PL]);
	// initial K = 0.0
	double U0 = 0;
	currentStepVals[pft_bEneN_PL] = U0 + currentStepVals[pft_bEneSource_PL]  - (K + U + currentStepVals[pft_iEned_PL]);
	// currentStepVals[pft_iEned]  = EneD = total energy dissipation
	// currentStepVals[pft_bEneSource] = is total energy from stress source (from vR term)
	// K(t) + U(t) - K(0) - U(0) = -EneD + EneSource - EneN

	if ((g_1_interface_low_disk_space == 0) || ((g_1_interface_low_disk_space == 2) && (timeIndexNew % 1000 == 0)))
	{
		outAllVals << timeIndexNew;
			for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
				outAllVals << '\t' << currentStepVals[i];
			outAllVals << '\n';
	}
	++timeIndexNew;
	for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
		stat4Vals[i].update(currentStepVals[i], timeNewAbsolute);

	PITSS retVal = PITSS_none;
	if (delu4Sigma0 - delu < 1e-6)
	{
		retVal = pit_sigma0;
		pits.push_back(pit_sigma0);
	}
	PITSS tmp = pit_sigmaMax;
	if ((!isExtrinsic) && (1.0 - delu < 1e-6) && (stageSlns[tmp].b_set == 0))
	{
		stageSlns[tmp].b_set = 1;
		pits.push_back(tmp);
		t_SigmaMax_real = timeNew;
	}
	vector<PITSS> tmps;
	tmp = PITSS_none;
	if (currentStepVals[pft_ivsolid] < 0.0)
	{
		tmp = pit_vSolidNegative;
		tmps.push_back(tmp);
	}
	if (1e-3 < 1.0 - reldus)
	{
		tmp = pit_reluSolid1;
		tmps.push_back(tmp);
	}
	if (reldus <= 0.0)
	{
		tmp = pit_reluSolid0;
		tmps.push_back(tmp);
	}
	if ((reldus + 1.0) <= 1e-3)
	{
//		retVal = pit_reluSolidm1;
		tmp = pit_reluSolidm1;
		tmps.push_back(tmp);
	}
	if (timeRelSigmaMax >= l_dilute)
	{
		tmp = pit_timeDilute;
		tmps.push_back(tmp);
	}
	if (timeRelSigmaMax >= l)
	{
		tmp = pit_timeInteraction;
		tmps.push_back(tmp);
	}
	for (unsigned int jj = 0; jj < tmps.size(); ++jj)
	{
		tmp = tmps[jj];
		if ((tmp != PITSS_none) && (stageSlns[tmp].b_set == 0))
			pits.push_back(tmp);
	}
	for (unsigned int i = 0; i < pits.size(); ++i)
		Updata_Periodic1IntrFrag_TimeStageStorage(pits[i]);
	if (terminateState != pit_sigma0)
		terminateState = retVal;

	set_energy_diss_per_length_4_t_dilute(timeNew, delu);
	return retVal;
}

void Periodic1IntrFrag::Finalize_Periodic1IntrFrag()
{
	string name = nameOne_a_One_l + "_outStats.txt";
	fstream out(name.c_str(), ios::out);
	for (unsigned int i = 0; i < OneSegmentPFT_SIZE; ++i)
	{
#if 0
		string name = nameOne_a_One_l + "_outStats_", ser;
		toString(i, ser);
		name += ser;
		name += "_";
		name += getName((OneSegmentPFT)i);
		name += ".txt";
		fstream out(name.c_str(), ios::out);
#endif
		out << stat4Vals[i];
	}
	if (stageSlns[pit_timeDilute].b_set == 0)
	{
		stageSlns[pit_timeDilute].CopyData(stageSlns[terminateState]);
		stageSlns[pit_timeDilute].b_set = 2;
	}
	if (stageSlns[pit_timeInteraction].b_set == 0)
	{
		stageSlns[pit_timeInteraction].CopyData(stageSlns[terminateState]);
		stageSlns[pit_timeInteraction].b_set = 2;
	}
	for (unsigned int i = 0; i < PITSS_SIZE; ++i)
		stageSlns[i].Print();
}

void Periodic1IntrFrag::Updata_Periodic1IntrFrag_TimeStageStorage(PITSS pitss)
{
	Periodic1IntrFrag_TimeStageStorage* stagePtr = &stageSlns[pitss];
	stagePtr->b_set = 1;
	stagePtr->vals = currentStepVals;
	stagePtr->spatial_stress_vals = stressVals;
	stagePtr->spatial_vel_vals = velVals;
	stagePtr->timeIndex = timeIndexNew;
}

void Periodic1IntrFrag::GetPrimaryFragmentSizes(vector<double>& primaryFragmentSizes)
{
	primaryFragmentSizes.resize(FragmentSizeT_SIZE);
	primaryFragmentSizes[fst_dilute] = l_dilute;
	primaryFragmentSizes[fst_dilute_approx] = l_dilute_approx;
	primaryFragmentSizes[fst_zhu6a] = l_Zhu6a;
	primaryFragmentSizes[fst_zhu6b] = l_Zhu6b;
	primaryFragmentSizes[fst_Glenn] = l_Glenn;
	primaryFragmentSizes[fst_Grady] = l_Grady;
}

void Periodic1IntrFrag::set_energy_diss_per_length_4_t_dilute(double time, double deltau)
{
	double timeAbsolute = time + t0;
	if ((tFailure < 0.0) && (delu4Sigma0 - deltau <= 1e-5))
		tFailure = timeAbsolute;
	if (energy_diss_per_length_At_t_dilute_set)
		return;
	if (timeAbsolute >= t_dilute)
	{
		energy_diss_per_length_At_t_dilute_set = true;
		energy_diss_per_length_At_t_dilute = currentStepVals[pft_iEned] / l;
		return;
	}
	if (diluteFractureModel)
	{
		energy_diss_per_length_At_t_dilute_set = true;
		energy_diss_per_length_At_t_dilute = energy_diss_per_length;
		return;
	}
}

void Periodic1IntrFrag::Output_Periodic1IntrFrag_Header(ostream& out)
{
	out << "terminateFlag" << "\t";
	out << "iterminateFlag" << "\t";
	
	out << "failureState" << "\t";
	out << "ifailureState" << "\t";

	out << "aIndex" << '\t';
	out << "lIndex" << '\t';
	out << "lIndex_secondary" << '\t';

	out << "a" << '\t';
	out << "l" << '\t';
	out << "t_SigmaMax_real" << '\t';
	out << "log10(a)" << '\t';
	out << "log10(l)" << '\t';
	out << "log10(t_SigmaMax_real)" << '\t';

	// dp base on value of 2 difference which is the normalization for both Zhu and Drugan
	// ap and bp is more reasonable as they use G. These used in Zhu06 and match b2 values for Ortiz model
	// but for Xu-Needleman a' = e * a, l' = l / e | whereas ab2 = 0.5 * a, lb2 = l / 0.5 = 2 l

	// 2: factors of two used for tau, a, ... regardless
	// G: G is used (e.g. 0.5 for Ortiz, e for Xu-Needleman)
	out << "a2p" << '\t';
	out << "l2p" << '\t';
	out << "t_SigmaZero_real_2p" << '\t';
	out << "log10(a2p)" << '\t';
	out << "log10(l2p)" << '\t';
	out << "log10(t_SigmaZero_real_2p)" << '\t';

	out << "aGp" << '\t';
	out << "lGp" << '\t';
	out << "t_SigmaZero_real_Gp" << '\t';
	out << "log10(aGp)" << '\t';
	out << "log10(lGp)" << '\t';
	out << "log10(t_SigmaZero_real_Gp)" << '\t';

	out << "t_SigmaZero_minus_SigmaMax_real" << '\t';
	out << "t_SigmaMax_dilute" << '\t';
	out << "t_SigmaZero_dilute" << '\t';
	out << "t_SigmaZero_minus_SigmaMax_dilute" << '\t';

	out << "t_dilute_approx" << '\t';
	out << "t_dilute" << '\t';
	out << "l_dilute_approx" << '\t';
	out << "l_dilute" << '\t';
	out << "l_SigmaRef_DruganMax_Dilute" << '\t';
	out << "l_SigmaRef_DruganFinal_Dilute" << '\t';
	out << "l_Zhu6a" << '\t';
	out << "l_Zhu6b" << '\t';
	out << "l_Grady" << '\t';
	out << "l_Glenn" << '\t';

	out << "t_SigmaMax_real_Gp" << '\t';
	out << "t_SigmaZero_minus_SigmaMax_real_Gp" << '\t';
	out << "t_SigmaMax_dilute_Gp" << '\t';
	out << "t_dilute_Gp" << '\t';
	out << "t_SigmaZero_minus_SigmaMax_dilute_Gp" << '\t';

	out << "t_dilute_approx_Gp" << '\t';
	out << "t_dilute_Gp" << '\t';
	out << "l_dilute_approx_Gp" << '\t';
	out << "l_dilute_Gp" << '\t';
	out << "l_SigmaRef_DruganMax_Dilute_Gp" << '\t';
	out << "l_SigmaRef_DruganFinal_Dilute_Gp" << '\t';
	out << "l_Zhu6a_Gp" << '\t';
	out << "l_Zhu6b_Gp" << '\t';
	out << "l_Grady_Gp" << '\t';
	out << "l_Glenn_Gp" << '\t';

	out << "energy_diss_per_length_At_t_dilute" << '\t';
	out << "log10(energy_diss_per_length_At_t_dilute)" << '\t';
	out << "energy_diss_per_length" << '\t';
	out << "log10(energy_diss_per_length)" << '\t';

	out << "inputEnergyFromMeanStress_PL" << '\t';
	out << "energy_diss_real_inp_E_t_final" << '\t';
	out << "log10(energy_diss_real_inp_E_t_final)" << '\t';

	out << "inputEnergyFromMeanStress_MaxPossVal_PL" << '\t';
	out << "energy_diss_max_inp_E_t_final" << '\t';
	out << "log10(energy_diss_max_inp_E_t_final)" << '\t';

	out << "t_vs0" << '\t';
	out << "delt_vs0Final" << '\t';
	out << "crosser_vs0Final" << '\t';
	out << "t_relus_0" << '\t';
	out << "delt_relus_0" << '\t';
	out << "crosser_relus_0" << '\t';
	out << "t_relus_m1" << '\t';
	out << "delt_relus_m1" << '\t';
	out << "crosser_relus_m1" << '\t';
	out << "tSigmaCMax" << '\t';
	out << "delt_vs0Max" << '\t';
	out << "crosser_vs0Max" << '\t';

	out << "diluteFractureModel" << '\t';
	out << "isExtrinsic" << '\t';
	out << "energyScale" << '\t';
	out << "zhu6_epsilonDotScale" << '\t';
	out << "zhu6_sBarScale" << '\t';
	out << "delu4Sigma0" << '\t';
	out << "la" << '\t';
	out << "aInv" << '\t';
	out << "sigma_t0" << '\t';
	out << "t0" << '\t';

//	out << "tFailure" << '\t';

	out << "failure_Ncycle" << '\t';
	out << "failure_relus" << '\t';
	out << "failure_vs" << '\t';
	out << "failure_bsigma" << '\t';
	out << "failure_bsigma_pwMax" << '\t';
	out << "failure_K" << '\t';
	out << "failure_U" << '\t';
	out << "failure_phi" << '\t';
	out << "failure_K2phi" << '\t';
	out << "failure_EneN" << '\t';

	out << "bsigmaMax_val" << '\t';
	out << "bsigmaMax_time" << '\t';
	out << "bsigmaMax_pw_val" << '\t';
	out << "bsigmaMax_pw_time" << '\t';

	// v solid zero
	out << "vs0_time" << '\t';
	out << "vs0_timex2" << '\t';
	out << "vs0_Ncycle" << '\t';
	out << "vs0_relus" << '\t';
	out << "vs0_vs" << '\t';
	out << "vs0_energy_diss_per_length" << '\t';
	out << "vs0_log10(energy_diss_per_length)" << '\t';
	out << "vs0_energy_diss/phi0" << '\t';
	out << "vs0_delu" << '\t';
	out << "vs0_delv" << '\t';
	out << "vs0_sigma" << '\t';
	out << "vs0_bulk_sigma" << '\t';
	out << "vs0_bulk_sigma_maxpw" << '\t';
	out << "vs0_K" << '\t';
	out << "vs0_U" << '\t';

	out << "suggeste_delt" << '\t';
	out << "suggested_final_time4ZeroSigmaC" << '\t';
	out << "t_dilute_zeroStressCheck\t";

	// added 07/15/2024
	// zero_sigc_approx_: cohesive unloading stress is treated as zero (approximate model) to come up with an exact solution
	out << "zero_sigc_approx_bsigmaMax_val" << '\t';
	out << "zero_sigc_approx_bsigmaMax_time" << '\t';
	out << "zero_sigc_approx_bsigmaMax_pw_val" << '\t';
	out << "zero_sigc_approx_bsigmaMax_pw_time" << '\t';

	out << "bsigmaMax_val_exact_2_approx" << '\t';
	out << "bsigmaMax_time_exact_2_approx" << '\t';
	out << "bsigmaMax_pw_val_exact_2_approx" << '\t';
	out << "bsigmaMax_pw_time_exact_2_approx" << '\t';

	out << "log10_failure_Ncycle" << '\t';
	
	// ratio of actual input energy to maximum it could have been at t_final
	out << "real_inp_E_2_max_inp_E_t_final" << '\t';
	// energy dissipation to input energy is printed above. These are K, U, phi (= K + U) and numerical error to real input energy
	out << "K_real_inp_E_t_final" << '\t';
	out << "U_real_inp_E_t_final" << '\t';
	out << "phi_real_inp_E_t_final" << '\t';
	out << "phiN_real_inp_E_t_final" << '\t';

	// energy dissipation to input energy is printed above. These are K, U, phi (= K + U) and numerical error to max possible input energy
	out << "K_max_inp_E_t_final" << '\t';
	out << "U_max_inp_E_t_final" << '\t';
	out << "phi_max_inp_E_t_final" << '\t';
	out << "phiN_max_inp_E_t_final" << '\t';
}

void Periodic1IntrFrag::Output_Periodic1IntrFrag(ostream& out)
{
	static string nan = "nan";
	PITSS pit = pit_sigma0;
	Periodic1IntrFrag_TimeStageStorage *stageSln = &stageSlns[pit];
	if (stageSln->b_set == 0)
	{
		pit = pit_reluSolidm1;
		stageSln = &stageSlns[pit];
	}
	if (stageSln->b_set == 0)
		THROW("Neither state is set, perhaps need to increase time\n");
	Periodic1IntrFrag_TimeStageStorage *stageSlnF = stageSln;

	double t_SigmaZero_real = stageSln->vals[pft_time];
	t_SigmaZero_minus_SigmaMax_real = t_SigmaZero_real - t_SigmaMax_real;

	static double nand = std::numeric_limits<double>::quiet_NaN();
	double tF_minus_tOFC, crosserValue, vF = 0.0, vtOFCbar = 0.0, tOFC = -1;
	double t_vs0 = nand, t_relus_0 = nand, t_relus_m1 = nand;
	bool OFC_v_decreasing = true;
	pit = pit_vSolidNegative;
	stageSln = &stageSlns[pit];

	// vs = 0 criterion (Drugan - max time)
	delt_vs0Max = std::numeric_limits<double>::quiet_NaN(), crosser_vs0Max = std::numeric_limits<double>::quiet_NaN();
	double tSigmaCMax = aInv;
	if (!isExtrinsic)
	{
		Periodic1IntrFrag_TimeStageStorage *stageSlnMax = &stageSlns[pit_sigmaMax];
		if (stageSlnMax->b_set)
		{
			double tF = stageSlnMax->vals[pft_time];
			tSigmaCMax = tF;
			vF = stageSlnMax->vals[pft_ivsolid];
			if (stageSln->b_set == 1)
			{
				tOFC = stageSln->vals[pft_time];
				t_vs0 = tOFC;
			}
			GetZeroCrosserValue(tF, tOFC, vF, vtOFCbar, tF_minus_tOFC, crosserValue, OFC_v_decreasing);
			delt_vs0Max = tF_minus_tOFC, crosser_vs0Max = crosserValue;
		}
	}

	// vs = 0 criterion (Drugan - final time)
	vF = stageSlnF->vals[pft_ivsolid];
	if (stageSln->b_set == 1)
	{
		tOFC = stageSln->vals[pft_time];
		t_vs0 = tOFC;
	}
	GetZeroCrosserValue(t_SigmaZero_real, tOFC, vF, vtOFCbar, tF_minus_tOFC, crosserValue, OFC_v_decreasing);
	delt_vs0Final = tF_minus_tOFC, crosser_vs0Final = crosserValue;

	// relus = 0 criterion (mine, bar back to its original length)
	pit = pit_reluSolid0;
	stageSln = &stageSlns[pit];
	vF = stageSlnF->vals[pft_irelusolid];
	tOFC = -1.0;
	if (stageSln->b_set == 1)
	{ 
		tOFC = stageSln->vals[pft_time];
		t_relus_0 = tOFC;
	}
	GetZeroCrosserValue(t_SigmaZero_real, tOFC, vF, vtOFCbar, tF_minus_tOFC, crosserValue, OFC_v_decreasing);
	delt_relus_0 = tF_minus_tOFC, crosser_relus_0 = crosserValue;

	// relus = -1 criterion (mine, bar goes to nondimensional strain = -1 (eps/(sigma_C/E) = -1)
	pit = pit_reluSolidm1;
	vtOFCbar = -1.0;
	stageSln = &stageSlns[pit];
	vF = stageSlnF->vals[pft_irelusolid];	// already retrieved above, but keeping it to be safe
	tOFC = -1.0;
	if (stageSln->b_set == 1)
	{ 
		tOFC = stageSln->vals[pft_time];
		t_relus_m1 = tOFC;
	}
	GetZeroCrosserValue(t_SigmaZero_real, tOFC, vF, vtOFCbar, tF_minus_tOFC, crosserValue, OFC_v_decreasing);
	delt_relus_m1 = tF_minus_tOFC, crosser_relus_m1 = crosserValue;
	stageSln = stageSlnF;
	out << setprecision(18);
	out << terminateState << "\t";
	out << (int)terminateState << "\t";
	out << pit << '\t';
	out << (int)pit << '\t';

	out << aIndex << '\t';
	out << lIndex << '\t';
	out << lIndex_secondary << '\t';

	inputEnergyFromMeanStress_PL = stageSln->vals[pft_bEneSource_PL];
	// 0.5 E (a tF)^2
	inputEnergyFromMeanStress_MaxPossVal_PL = (a * t_SigmaZero_real);
	inputEnergyFromMeanStress_MaxPossVal_PL *= (0.5 * inputEnergyFromMeanStress_MaxPossVal_PL);

	energy_diss_real_inp_E_t_final = energy_diss_per_length / inputEnergyFromMeanStress_PL;
	energy_diss_max_inp_E_t_final = energy_diss_per_length / inputEnergyFromMeanStress_MaxPossVal_PL;

	// max value that I could have taken above, if sigma = at -> IMax = a^2t^2 E / 2 -- t is the final time

	out << a << '\t';
	out << l << '\t';
	out << t_SigmaZero_real << '\t';
	out << rlog10(a) << '\t';
	out << rlog10(l) << '\t';
	out << rlog10(t_SigmaZero_real) << '\t';

	double timeFactor = 1.0 / energyScale;

	// b2 base on value of 2 difference which is the normalization for both Zhu and Drugan
	// ap and bp is more reasonable as they use G. These used in Zhu06 and match b2 values for Ortiz model
	// but for Xu-Needleman a' = e * a, l' = l / e | whereas ab2 = 0.5 * a, lb2 = l / 0.5 = 2 l
	double ab2 = 0.5 * a;
	double lb2 = 2.0 * l;
	out << ab2 << '\t';
	out << lb2 << '\t';
	out << 2.0 * t_SigmaZero_real << '\t';
	out << rlog10(ab2) << '\t';
	out << rlog10(lb2) << '\t';
	out << rlog10(2.0 * t_SigmaZero_real) << '\t';

	out << a_p << '\t';
	out << l_p << '\t';
	out << timeFactor * t_SigmaZero_real << '\t';
	out << rlog10(a_p) << '\t';
	out << rlog10(l_p) << '\t';
	out << rlog10(timeFactor * t_SigmaZero_real) << '\t';

	out << t_SigmaZero_minus_SigmaMax_real << '\t';
	out << t_SigmaMax_dilute << '\t';
	out << t_dilute << '\t';
	out << t_SigmaZero_minus_SigmaMax_dilute << '\t';
	out << t_dilute_approx << '\t';
	out << t_dilute << '\t';
	out << l_dilute_approx << '\t';
	out << l_dilute << '\t';
	out << l_SigmaRef_DruganMax_Dilute << '\t';
	out << l_SigmaRef_DruganFinal_Dilute << '\t';
	out << l_Zhu6a << '\t';
	out << l_Zhu6b << '\t';
	out << l_Grady << '\t';
	out << l_Glenn << '\t';

	out << timeFactor * t_SigmaMax_real << '\t';
	out << timeFactor * t_SigmaZero_minus_SigmaMax_real << '\t';
	out << timeFactor * t_SigmaMax_dilute << '\t';
	out << timeFactor * t_dilute << '\t';
	out << timeFactor * t_SigmaZero_minus_SigmaMax_dilute << '\t';
	out << t_dilute_approx_p << '\t';
	out << t_dilute_p << '\t';
	out << l_dilute_approx_p << '\t';
	out << l_dilute_p << '\t';
	double factor = 1.0 / zhu6_sBarScale;
	out << factor * l_SigmaRef_DruganMax_Dilute << '\t';
	out << factor * l_SigmaRef_DruganFinal_Dilute << '\t';
	out << factor * l_Zhu6a << '\t';
	out << factor * l_Zhu6b << '\t';
	out << factor * l_Grady << '\t';
	out << factor * l_Glenn << '\t';

	out << energy_diss_per_length_At_t_dilute << '\t';
	out << rlog10(energy_diss_per_length_At_t_dilute) << '\t';
	out << energy_diss_per_length << '\t';
	out << rlog10(energy_diss_per_length) << '\t';

	out << inputEnergyFromMeanStress_PL << '\t';
	out << energy_diss_real_inp_E_t_final << '\t';
	out << rlog10(energy_diss_real_inp_E_t_final) << '\t';

	out << inputEnergyFromMeanStress_MaxPossVal_PL << '\t';
	out << energy_diss_max_inp_E_t_final << '\t';
	out << rlog10(energy_diss_max_inp_E_t_final) << '\t';

	out << t_vs0 << '\t';
	out << delt_vs0Final << '\t';
	out << crosser_vs0Final << '\t';
	out << t_relus_0 << '\t';
	out << delt_relus_0 << '\t';
	out << crosser_relus_0 << '\t';
	out << t_relus_m1 << '\t';
	out << delt_relus_m1 << '\t';
	out << crosser_relus_m1 << '\t';
	out << tSigmaCMax << '\t';
	out << delt_vs0Max << '\t';
	out << crosser_vs0Max << '\t';

	out << diluteFractureModel << '\t';
	out << isExtrinsic << '\t';
	out << energyScale << '\t';
	out << zhu6_epsilonDotScale << '\t';
	out << zhu6_sBarScale << '\t';
	out << delu4Sigma0 << '\t';
	out << la << '\t';
	out << aInv << '\t';
	out << sigma_t0 << '\t';
	out << t0 << '\t';

	double failure_Ncycle = stageSln->vals[pft_numCyclesAfterCrackOpening];
	out << failure_Ncycle << '\t';
	out << stageSln->vals[pft_irelusolid] << '\t';
	out << stageSln->vals[pft_ivsolid] << '\t';
	out << stageSln->vals[pft_bsigma] << '\t';
	out << stageSln->vals[pft_bsigma_maxpw] << '\t';
	double K = stageSln->vals[pft_bK_PL], U = stageSln->vals[pft_bU_PL];
	double phi = K + U, K2phi = K / phi;
	out << K << '\t';
	out << U << '\t';
	out << phi << '\t';
	out << K2phi << '\t';
	double phiN = stageSln->vals[pft_bEneN_PL];
	out << phiN << '\t';

	double bsigmaMax_val = std::numeric_limits<double>::quiet_NaN();
	double bsigmaMax_time = std::numeric_limits<double>::quiet_NaN();
	double bsigmaMax_pw_val = std::numeric_limits<double>::quiet_NaN();
	double bsigmaMax_pw_time = std::numeric_limits<double>::quiet_NaN();

	bool sigmaMaxAveSet = false;
	if (stageSlns[pit_timeSigmaAveMax].b_set != 0)
	{
		bsigmaMax_val = stageSlns[pit_timeSigmaAveMax].vals[pft_bsigma];
		bsigmaMax_time = stageSlns[pit_timeSigmaAveMax].vals[pft_time];
		out << bsigmaMax_val << '\t' << bsigmaMax_time << '\t';
		sigmaMaxAveSet = true;
	}
	else
		out << nan << '\t' << nan << '\t';

	bool sigmaMaxPWSet = false;
	if (stageSlns[pit_timeSigmaPWMax].b_set != 0)
	{
		bsigmaMax_pw_val = stageSlns[pit_timeSigmaPWMax].vals[pft_bsigma_maxpw];
		bsigmaMax_pw_time = stageSlns[pit_timeSigmaPWMax].vals[pft_time];
		out << bsigmaMax_pw_val << '\t' << bsigmaMax_pw_time << '\t';
		sigmaMaxPWSet = true;
	}
	else
		out << nan << '\t' << nan << '\t';

	pit = pit_vSolidNegative;
	stageSln = &stageSlns[pit];
	if (stageSln->b_set == 1)
	{
		stageSln = &stageSlns[pit];
		out << stageSln->vals[pft_time] << '\t';
		out << 2.0 * stageSln->vals[pft_time] << '\t';
		out << stageSln->vals[pft_numCyclesAfterCrackOpening] << '\t';
		out << stageSln->vals[pft_irelusolid] << '\t';
		out << stageSln->vals[pft_ivsolid] << '\t';
		out << stageSln->vals[pft_iEned_PL] << '\t';
		out << stageSln->vals[pft_ilog10Ened_PL] << '\t';
		out << stageSln->vals[pft_iEned] / energyScale << '\t';
		out << stageSln->vals[pft_idelu] << '\t';
		out << stageSln->vals[pft_idelv] << '\t';
		out << stageSln->vals[pft_isigma] << '\t';
		out << stageSln->vals[pft_bsigma] << '\t';
		out << stageSln->vals[pft_bsigma_maxpw] << '\t';
		out << stageSln->vals[pft_bK_PL] << '\t';
		out << stageSln->vals[pft_bU_PL] << '\t';
	}
	else
	{
		out << -1.0 << '\t';
		for (unsigned int i = 0; i < 13; ++i)
			out << nan << '\t';
	}
	out << suggeste_delt << '\t';
	out << suggested_final_time4ZeroSigmaC << '\t';
	out << t_dilute_zeroStressCheck << '\t';

	// 2024/07/15
	double time_a = 1 / a, time_l = l; // l / c but c = 1
	double zero_sigc_approx_bsigmaMax_pw_time = time_a + 0.5 * time_l;
	double zero_sigc_approx_bsigmaMax_pw_val = 1.0 + 0.5 * time_l / time_a;
	double zero_sigc_approx_bsigmaMax_time = time_a;
	double zero_sigc_approx_bsigmaMax_val = 1.0;
	double delt_4sigmax = 0.25 * (time_l - 2.0 * time_a);
	if (delt_4sigmax > 0.0)
	{
		zero_sigc_approx_bsigmaMax_time += delt_4sigmax;
		zero_sigc_approx_bsigmaMax_val += 2.0 * delt_4sigmax / time_a / time_l;
	}
	out << zero_sigc_approx_bsigmaMax_val << '\t';
	out << zero_sigc_approx_bsigmaMax_time << '\t';
	out << zero_sigc_approx_bsigmaMax_pw_val << '\t';
	out << zero_sigc_approx_bsigmaMax_pw_time << '\t';

	if (sigmaMaxAveSet)
	{
		double bsigmaMax_val_exact_2_approx = bsigmaMax_val / zero_sigc_approx_bsigmaMax_val;
		double bsigmaMax_time_exact_2_approx = bsigmaMax_time / zero_sigc_approx_bsigmaMax_time;
		out << bsigmaMax_val_exact_2_approx << '\t';
		out << bsigmaMax_time_exact_2_approx << '\t';
	}
	else
	{
		out << nan << '\t' << nan << '\t';
	}
	if (sigmaMaxPWSet)
	{
		double bsigmaMax_pw_val_exact_2_approx = bsigmaMax_pw_val / zero_sigc_approx_bsigmaMax_pw_val;
		double bsigmaMax_pw_time_exact_2_approx = bsigmaMax_pw_time / zero_sigc_approx_bsigmaMax_pw_time;
		out << bsigmaMax_pw_val_exact_2_approx << '\t';
		out << bsigmaMax_pw_time_exact_2_approx << '\t';
	}
	else
	{
		out << nan << '\t' << nan << '\t';
	}
	//log10_failure_Ncycle
	out << rlog10(failure_Ncycle) << '\t';

	// ratio of actual input energy to maximum it could have been at t_final
	double real_inp_E_2_max_inp_E_t_final = inputEnergyFromMeanStress_PL / inputEnergyFromMeanStress_MaxPossVal_PL;
	out << real_inp_E_2_max_inp_E_t_final << '\t';

	// energy dissipation to input energy is printed above. These are K, U, phi (= K + U) and numerical error to real input energy
	double K_real_inp_E_t_final = K / inputEnergyFromMeanStress_PL;
	out << K_real_inp_E_t_final << '\t';
	double U_real_inp_E_t_final = U / inputEnergyFromMeanStress_PL;
	out << U_real_inp_E_t_final << '\t';
	double phi_real_inp_E_t_final = phi / inputEnergyFromMeanStress_PL;
	out << phi_real_inp_E_t_final << '\t';
	double phiN_real_inp_E_t_final = phiN / inputEnergyFromMeanStress_PL;
	out << phiN_real_inp_E_t_final << '\t';

	// energy dissipation to input energy is printed above. These are K, U, phi (= K + U) and numerical error to max possible input energy
	double K_max_inp_E_t_final = K / inputEnergyFromMeanStress_MaxPossVal_PL;
	out << K_max_inp_E_t_final << '\t';
	double U_max_inp_E_t_final = U / inputEnergyFromMeanStress_MaxPossVal_PL;
	out << U_max_inp_E_t_final << '\t';
	double phi_max_inp_E_t_final = phi / inputEnergyFromMeanStress_MaxPossVal_PL;
	out << phi_max_inp_E_t_final << '\t';
	double phiN_max_inp_E_t_final = phiN / inputEnergyFromMeanStress_MaxPossVal_PL;
	out << phiN_max_inp_E_t_final << '\n';
}



