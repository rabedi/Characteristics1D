#include"SLFracTypes.h"
#include"commonMacros.h"

bool IsMC(SLFF_Eeffective_StsType effS)
{
	return ((effS == effSs_mc) || (effS == effSs_mct) || (effS == effSs_mctc));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_scalarRandomMultiplierType dat, string& name)
{
	if (dat == slff_sclrm_none)
	{
		name = "none";
		return;
	}
	if (dat == slff_sclrm_sx_dconst)
	{
		name = "sx_dconst";
		return;
	}
	if (dat == slff_sclrm_sx_Gconst)
	{
		name = "sx_Gconst";
		return;
	}
	if (dat == slff_sclrm_sx_tauconst)
	{
		name = "sx_tauconst";
		return;
	}
	if (dat == slff_sclrm_Gx_sconst)
	{
		name = "Gx_sconst";
		return;
	}
	if (dat == slff_sclrm_Gx_dconst)
	{
		name = "Gx_dconst";
		return;
	}
	if (dat == slff_sclrm_Gx_tauconst)
	{
		name = "Gx_tauconst";
		return;
	}
	if (dat == slff_sclrm_ICDamagex)
	{
		name = "ICDamagex";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_scalarRandomMultiplierType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_scalarRandomMultiplierType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_scalarRandomMultiplierType_SIZE; ++i)
	{
		getName(SLFF_scalarRandomMultiplierType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_scalarRandomMultiplierType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_scalarRandomMultiplierType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_scalarRandomMultiplierType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_BulkDamagePFType dat, string& name)
{
	if (dat == noLengthBulkDamage)
	{
		name = "noLengthBulkDamage";
		return;
	}
	if (dat == AllixTimeDelay)
	{
		name = "AllixTimeDelay";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_BulkDamagePFType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_BulkDamagePFType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_BulkDamagePFType_SIZE; ++i)
	{
		getName(SLFF_BulkDamagePFType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_BulkDamagePFType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_BulkDamagePFType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_BulkDamagePFType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFracPropT dat, string& name)
{
	if (dat == SLFF_rho)
	{
		name = "SLFF_rho";
		return;
	}
	if (dat == SLFF_E)
	{
		name = "SLFF_E";
		return;
	}
	if (dat == SLFF_nu)
	{
		name = "SLFF_nu";
		return;
	}
	if (dat == SLFF_Cn)
	{
		name = "SLFF_Cn";
		return;
	}
	if (dat == SLFF_Ct)
	{
		name = "SLFF_Ct";
		return;
	}
	if (dat == SLFF_Zn)
	{
		name = "SLFF_Zn";
		return;
	}
	if (dat == SLFF_Zt)
	{
		name = "SLFF_Zt";
		return;
	}
	if (dat == SLFF_cn)
	{
		name = "SLFF_cn";
		return;
	}
	if (dat == SLFF_ct)
	{
		name = "SLFF_ct";
		return;
	}
	if (dat == SLFF_cR)
	{
		name = "SLFF_cR";
		return;
	}
	if (dat == SLF_s)
	{
		name = "s";
		return;
	}
	if (dat == SLF_sn)
	{
		name = "sn";
		return;
	}
	if (dat == SLF_st)
	{
		name = "st";
		return;
	}
	if (dat == SLF_sc)
	{
		name = "sc";
		return;
	}
	if (dat == SLF_c)
	{
		name = "c";
		return;
	}
	if (dat == SLF_beta)
	{
		name = "beta";
		return;
	}
	if (dat == SLF_eta)
	{
		name = "eta";
		return;
	}
	if (dat == SLF_brittleness)
	{
		name = "brittleness";
		return;
	}
	if (dat == SLF_Gn)
	{
		name = "Gn";
		return;
	}
	if (dat == SLF_Gt)
	{
		name = "Gt";
		return;
	}
	if (dat == SLF_YScaleElastic)
	{
		name = "YScaleElastic";
		return;
	}
	if (dat == SLF_YScalePF)
	{
		name = "YScalePF";
		return;
	}
	if (dat == SLF_dn)
	{
		name = "dn";
		return;
	}
	if (dat == SLF_dt)
	{
		name = "dt";
		return;
	}
	if (dat == SLF_bn)
	{
		name = "bn";
		return;
	}
	if (dat == SLF_bt)
	{
		name = "bt";
		return;
	}
	if (dat == SLF_dvn)
	{
		name = "dvn";
		return;
	}
	if (dat == SLF_dvt)
	{
		name = "dvt";
		return;
	}
	if (dat == SLF_lOverbn)
	{
		name = "lOverbn";
		return;
	}
	if (dat == SLF_lOverbt)
	{
		name = "lOverbt";
		return;
	}
	if (dat == SLF_taun)
	{
		name = "taun";
		return;
	}
	if (dat == SLF_taut)
	{
		name = "taut";
		return;
	}
//	if (dat == SLF_tauRateModel)
//	{
//		name = "tauRateModel";
//		return;
//	}
	if (dat == SLF_ln)
	{
		name = "ln";
		return;
	}
	if (dat == SLF_lt)
	{
		name = "lt";
		return;
	}
	if (dat == SLF_vn)
	{
		name = "vn";
		return;
	}
	if (dat == SLF_vt)
	{
		name = "vt";
		return;
	}
	if (dat == SLF_dcont)
	{
		name = "dcont";
		return;
	}
	if (dat == SLF_dsep)
	{
		name = "dsep";
		return;
	}
	if (dat == SLF_k)
	{
		name = "k";
		return;
	}
	if (dat == SLF_damaged_kUpper)
	{
		name = "damaged_kUpper";
		return;
	}
	if (dat == SLF_damaged_kLower)
	{
		name = "damaged_kLower";
		return;
	}
	if (dat == SLF_damaged_kUpper_para1)
	{
		name = "damaged_kUpper_para1";
		return;
	}
	if (dat == SLF_damaged_kLower_para1)
	{
		name = "damaged_kLower_para1";
		return;
	}
	if (dat == SLF_damaged_kUpper_para2)
	{
		name = "damaged_kUpper_para2";
		return;
	}
	if (dat == SLF_damaged_kLower_para2)
	{
		name = "damaged_kLower_para2";
		return;
	}
	if (dat == SLF_damaged_kUpper_para3)
	{
		name = "damaged_kUpper_para3";
		return;
	}
	if (dat == SLF_damaged_kLower_para3)
	{
		name = "damaged_kLower_para3";
		return;
	}
	if (dat == SLF_damaged_para0_epsSW)
	{
		name = "damaged_para0_epsSW";
		return;
	}
	if (dat == SLF_damaged_para1_n)
	{
		name = "damaged_para1_n";
		return;
	}
	if (dat == SLF_damaged_para2)
	{
		name = "damaged_para2";
		return;
	}
	if (dat == SLF_damaged_c)
	{
		name = "damaged_c";
		return;
	}
	if (dat == SLF_damaged_sn)
	{
		name = "damaged_sn";
		return;
	}
	if (dat == SLF_phi)
	{
		name = "phi";
		return;
	}
	if (dat == SLF_phiDegree)
	{
		name = "phiDegree";
		return;
	}
	if (dat == SLF_dt_slip)
	{
		name = "dt_slip";
		return;
	}
	if (dat == SLF_beta_delU)
	{
		name = "beta_delU";
		return;
	}
	if (dat == SLF_factor_k)
	{
		name = "factor_k";
		return;
	}
	if (dat == SLF_td_kn)
	{
		name = "td_kn";
		return;
	}
	if (dat == SLF_td_kc)
	{
		name = "td_kc";
		return;
	}
	if (dat == SLF_Hpara0_Allix_a)
	{
		name = "Hpara0_Allix_a";
		return;
	}
	if (dat == SLF_pf_timescale)
	{
		name = "pf_timescale";
		return;
	}
	if (dat == SLF_pf_rhoTilde)
	{
		name = "pf_rhoTilde";
		return;
	}
	if (dat == SLF_pf_dTilde)
	{
		name = "pf_dTilde";
		return;
	}
	if (dat == SLF_pf_ETilde)
	{
		name = "pf_ETilde";
		return;
	}
	if (dat == SLF_InSitu_sn)
	{
		name = "InSitu_sn";
		return;
	}
	if (dat == SLF_InSitu_st0)
	{
		name = "InSitu_st0";
		return;
	}
	if (dat == SLF_InSitu_st1)
	{
		name = "InSitu_st1";
		return;
	}
	name = "none";
//	cout << "dat2\t" << (int)dat << '\n';
}

ostream& operator<<(ostream& out, SLFracPropT  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFracPropT& dat)
{
	string name2;
	for (int i = 0; i < SLFracPropT_SIZE; ++i)
	{
		getName(SLFracPropT(i), name2);
		if (name2 == name)
		{
			dat = SLFracPropT(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFracPropT& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFracFlagT dat, string& name)
{
	if (dat == SLFF_scalarRandomMultiplierV)
	{
		name = "scalarRandomMultiplier";
		return;
	}
	if (dat == SLFF_isBulkDamage)
	{
		name = "isBulkDamage";
		return;
	}
	if (dat == SLFF_BulkDamagePFV)
	{
		name = "BulkDamagePF";
		return;
	}
	if (dat == SLFF_Eeffective_StsV)
	{
		name = "Eeffective_Sts";
		return;
	}
	if (dat == SLFF_Eeffective_StnV)
	{
		name = "Eeffective_Stn";
		return;
	}
	if (dat == SLFF_EpsEl2SigmaV)
	{
		name = "EpsEl2Sigma";
		return;
	}
	if (dat == SLFF_EpsE2YbarV)
	{
		name = "EpsE2Ybar";
		return;
	}
	if (dat == SLFF_YbarEnergyV)
	{
		name = "YbarEnergy";
		return;
	}
	if (dat == SLFF_OmegaDV)
	{
		name = "OmegaD";
		return;
	}
	if (dat == SLFF_AlphaDV)
	{
		name = "AlphaD";
		return;
	}
	if (dat == SLFF_TSRV)
	{
		name = "TSR";
		return;
	}
	if (dat == SLFF_PF_velOptionV)
	{
		name = "PF_velOptionV";
		return;
	}
	if (dat == SLFF_PF_fieldOptionSrcLV)
	{
		name = "PF_fieldOptionSrcLV";
		return;
	}
	if (dat == SLFF_InterfaceTypeV)
	{
		name = "InterfaceType";
		return;
	}
	if (dat == SLFF_InterfaceHasInSituV)
	{
		name = "InterfaceHasInSitu";
		return;
	}
	if (dat == SLFF_ForceU_SlipSlnsV)
	{
		name = "ForceU_SlipSlns";
		return;
	}
	if (dat == SLFF_InterfacialDamageModeV)
	{
		name = "InterfacialDamageMode";
		return;
	}
	if (dat == SLFF_InterfacialDamageTractionV)
	{
		name = "InterfacialDamageTraction";
		return;
	}
//	if (dat == NormTanCompV)
//	{
//		name = "NormTanComp";
//		return;
//	}
	if (dat == SLFF_ContactV)
	{
		name = "Contact";
		return;
	}
	if (dat == SLFF_SlipV)
	{
		name = "Slip";
		return;
	}
	if (dat == SLFF_FrictionV)
	{
		name = "Friction";
		return;
	}
	if (dat == SLFF_FrictionPropertiesAreIso)
	{
		name = "FrictionPropertiesAreIso";
		return;
	}
	if (dat == SLFF_ContactSeparationRegularizationV)
	{
		name = "ContactSeparationRegularizationV";
		return;
	}
	if (dat == SLFF_funHV)
	{
		name = "funH";
		return;
	}
	if (dat == SLFF_RBH_VelDDotV)
	{
		name = "RBH_VelDDotV";
		return;
	}
	if (dat == SLFF_applyFactor_k_2YReductionB)
	{
		name = "applyFactor_k_2YReduction";
		return;
	}
	if (dat == SLFF_applyYReductionBy_wPrime_dB)
	{
		name = "applyYReductionBy_wPrime_d";
		return;
	}
	if (dat == SLFF_DamageTarget_boundedBy1B)
	{
		name = "DamageTarget_boundedBy1";
		return;
	}
	if (dat == SLFF_td_pn)
	{
		name = "td_pn";
		return;
	}
	if (dat == SLFF_td_pc)
	{
		name = "td_pc";
		return;
	}
	name = "none";
//	cout << "dat\t" << (int)dat << '\n';
}

ostream& operator<<(ostream& out, SLFracFlagT  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFracFlagT& dat)
{
	string name2;
	for (int i = 0; i < SLFracFlagT_SIZE; ++i)
	{
		getName(SLFracFlagT(i), name2);
		if (name2 == name)
		{
			dat = SLFracFlagT(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFracFlagT& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_Eeffective_StsType dat, string& name)
{
	if (dat == effSs_Off)
	{
		name = "effSs_Off";
		return;
	}
	if (dat == effSs_sqrt)
	{
		name = "effSs_sqrt";
		return;
	}
	if (dat == effSs_mc)
	{
		name = "effSs_mc";
		return;
	}
	if (dat == effSs_mct)
	{
		name = "effSs_mct";
		return;
	}
	if (dat == effSs_mctc)
	{
		name = "effSs_mctc";
		return;
	}
	if (dat == effSs_HoekBrown)
	{
		name = "effSs_HoekBrown";
		return;
	}
	if (dat == effSs_Rankin)
	{
		name = "effSs_Rankin";
		return;
	}
	if (dat == effSs_Tresca)
	{
		name = "effSs_Tresca";
		return;
	}
	if (dat == effSsCompressive)
	{
		name = "effSsCompressive";
		return;
	}
	if (dat == effSs_vm)
	{
		name = "effSs_vm";
		return;
	}
	if (dat == effSs_vmm)
	{
		name = "effSs_vmm";
		return;
	}
	if (dat == effSs_DP)
	{
		name = "effSs_DP";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_Eeffective_StsType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_Eeffective_StsType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_Eeffective_StsType_SIZE; ++i)
	{
		getName(SLFF_Eeffective_StsType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_Eeffective_StsType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_Eeffective_StsType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_Eeffective_StnType dat, string& name)
{
	if (dat == effSn_Off)
	{
		name = "effSn_Off";
		return;
	}
	if (dat == effSn_Mazars)
	{
		name = "effSn_Mazars";
		return;
	}
	if (dat == effSn_maxEps1)
	{
		name = "effSn_maxEps1";
		return;
	}
	if (dat == effSn_vmm_Peerling)
	{
		name = "effSn_vmm_Peerling";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_Eeffective_StnType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_Eeffective_StnType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_Eeffective_StnType_SIZE; ++i)
	{
		getName(SLFF_Eeffective_StnType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_Eeffective_StnType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_Eeffective_StnType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_EpsEl2SigmaType dat, string& name)
{
	if (dat == eps2sig_noBreakDown)
	{
		name = "noBreakDown";
		return;
	}
	if (dat == eps2sig_epsmp)
	{
		name = "epsmp";
		return;
	}
	if (dat == eps2sig_sigmp)
	{
		name = "sigmp";
		return;
	}
	if (dat == eps2sig_epshd)
	{
		name = "epshd";
		return;
	}
	if (dat == eps2sig_sighd)
	{
		name = "sighd";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_EpsEl2SigmaType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_EpsEl2SigmaType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_EpsEl2SigmaType_SIZE; ++i)
	{
		getName(SLFF_EpsEl2SigmaType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_EpsEl2SigmaType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_EpsEl2SigmaType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_EpsE2YbarType dat, string& name)
{
	if (dat == Ybar_energy)
	{
		name = "energy";
		return;
	}
	if (dat == Ybar_sqrtEnergy)
	{
		name = "sqrtEnergy";
		return;
	}
	if (dat == Ybar_stress)
	{
		name = "stress";
		return;
	}
	if (dat == Ybar_strain)
	{
		name = "strain";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_EpsE2YbarType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_EpsE2YbarType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_EpsE2YbarType_SIZE; ++i)
	{
		getName(SLFF_EpsE2YbarType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_EpsE2YbarType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_EpsE2YbarType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_YbarEnergyType dat, string& name)
{
	if (dat == Ybar_energyOff)
	{
		name = "energyOff";
		return;
	}
	if (dat == YbarStrnStrsBDE)
	{
		name = "YbarStrnStrsBDE";
		return;
	}
	if (dat == Ybar_energy_TotalStress)
	{
		name = "energy_TotalStress";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_YbarEnergyType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_YbarEnergyType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_YbarEnergyType_SIZE; ++i)
	{
		getName(SLFF_YbarEnergyType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_YbarEnergyType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_YbarEnergyType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_OmegaDType dat, string& name)
{
	if (dat == omegad_1md)
	{
		name = "1md";
		return;
	}
	if (dat == omegad_1md2)
	{
		name = "1md2";
		return;
	}
	if (dat == omegad_rational)
	{
		name = "rational";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_OmegaDType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_OmegaDType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_OmegaDType_SIZE; ++i)
	{
		getName(SLFF_OmegaDType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_OmegaDType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_OmegaDType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_AlphaDType dat, string& name)
{
	if (dat == alphad_xi_0_AT2)
	{
		name = "xi_0_AT2";
		return;
	}
	if (dat == alphad_xi_1_AT1)
	{
		name = "xi_1_AT1";
		return;
	}
	if (dat == alphad_xi_2_PZ_CZM)
	{
		name = "xi_2_PZ_CZM";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_AlphaDType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_AlphaDType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_AlphaDType_SIZE; ++i)
	{
		getName(SLFF_AlphaDType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_AlphaDType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_AlphaDType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_TSRType dat, string& name)
{
	if (dat == tsr_Zero)
	{
		name = "tsr_Zero";
		return;
	}
	if (dat == tsr_Xu_Needleman)
	{
		name = "tsr_Xu_Needleman";
		return;
	}
	if (dat == tsr_Ortiz)
	{
		name = "tsr_Ortiz";
		return;
	}
	if (dat == tsr_linear)
	{
		name = "tsr_linear";
		return;
	}
	if (dat == tsr_bilinear)
	{
		name = "tsr_bilinear";
		return;
	}
	if (dat == tsr_exponential)
	{
		name = "tsr_exponential";
		return;
	}
	if (dat == tsr_hyperbolic)
	{
		name = "tsr_hyperbolic";
		return;
	}
	if (dat == tsr_concrete)
	{
		name = "tsr_concrete";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_TSRType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_TSRType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_TSRType_SIZE; ++i)
	{
		getName(SLFF_TSRType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_TSRType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_TSRType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

bool IsExtrinsic(SLFF_TSRType dat)
{
	if ((dat == tsr_Ortiz) || (dat == tsr_linear))
		return true;
	return false;
}

double GetEnergyConstantFactor(SLFF_TSRType dat)
{
	if ((dat == tsr_Ortiz) || (dat == tsr_linear))
		return 0.5;
	if (dat == tsr_Xu_Needleman)
		return exp(1.0);
	THROW("Option not implemented\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_PF_velOptionType dat, string& name)
{
	if (dat == pfv_cd)
	{
		name = "pfv_cd";
		return;
	}
	if (dat == pfv_cs)
	{
		name = "pfv_cs";
		return;
	}
	if (dat == pfv_cr)
	{
		name = "pfv_cr";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_PF_velOptionType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_PF_velOptionType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_PF_velOptionType_SIZE; ++i)
	{
		getName(SLFF_PF_velOptionType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_PF_velOptionType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_PF_velOptionType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_PF_fieldOptionSrcLType dat, string& name)
{
	if (dat == pfsl_history)
	{
		name = "history";
		return;
	}
	if (dat == pfsl_pf)
	{
		name = "pf";
		return;
	}
	if (dat == pfsl_pf_bounded)
	{
		name = "pf_bounded";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_PF_fieldOptionSrcLType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_PF_fieldOptionSrcLType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_PF_fieldOptionSrcLType_SIZE; ++i)
	{
		getName(SLFF_PF_fieldOptionSrcLType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_PF_fieldOptionSrcLType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_PF_fieldOptionSrcLType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(NormTanCompType dat, string& name)
{
	if (dat == ntc_n)
	{
		name = "norm";
		return;
	}
	if (dat == ntc_t)
	{
		name = "tang";
		return;
	}
	if (dat == ntc_c)
	{
		name = "comp";
		return;
	}
	if (dat == ntc_none)
	{
		name = "none";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, NormTanCompType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, NormTanCompType& dat)
{
	string name2;
	for (int i = 0; i < NormTanCompType_SIZE; ++i)
	{
		getName(NormTanCompType(i), name2);
		if (name2 == name)
		{
			dat = NormTanCompType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, NormTanCompType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_ForceU_SlipSlnsType dat, string& name)
{
	if (dat == sl_forceU0_forceSlip0)
	{
		name = "forceU0_forceSlip0";
		return;
	}
	if (dat == sl_forceU1_forceSlip0)
	{
		name = "forceU1_forceSlip0";
		return;
	}
	if (dat == sl_forceU1_forceSlip1)
	{
		name = "forceU1_forceSlip1";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_ForceU_SlipSlnsType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_ForceU_SlipSlnsType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_ForceU_SlipSlnsType_SIZE; ++i)
	{
		getName(SLFF_ForceU_SlipSlnsType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_ForceU_SlipSlnsType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_ForceU_SlipSlnsType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_ForceU_SlipSlnsType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_InterfacialDamageModeType dat, string& name)
{
	if (dat == sl_interfacial_damage_off)
	{
		name = "interfacial_damage_off";
		return;
	}
	if (dat == sl_interfacial_damage_on)
	{
		name = "interfacial_damage_on";
		return;
	}
	if (dat == sl_interfacial_damage_mixed)
	{
		name = "interfacial_damage_mixed";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_InterfacialDamageModeType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_InterfacialDamageModeType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_InterfacialDamageModeType_SIZE; ++i)
	{
		getName(SLFF_InterfacialDamageModeType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_InterfacialDamageModeType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_InterfacialDamageModeType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_InterfacialDamageModeType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_InterfacialDamageTractionType dat, string& name)
{
	if (dat == sl_interfacial_damage_traction_zero)
	{
		name = "interfacial_damage_traction_zero";
		return;
	}
	if (dat == sl_interfacial_damage_traction_TSR)
	{
		name = "interfacial_damage_traction_TSR";
		return;
	}
	if (dat == sl_interfacial_damage_traction_HF)
	{
		name = "interfacial_damage_traction_HF";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_InterfacialDamageTractionType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_InterfacialDamageTractionType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_InterfacialDamageTractionType_SIZE; ++i)
	{
		getName(SLFF_InterfacialDamageTractionType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_InterfacialDamageTractionType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_InterfacialDamageTractionType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_InterfacialDamageTractionType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_ContactType dat, string& name)
{
	if (dat == sl_contact_off)
	{
		name = "contact_off";
		return;
	}
	if (dat == sl_contact_on)
	{
		name = "contact_on";
		return;
	}
	if (dat == sl_contact_mixed)
	{
		name = "contact_mixed";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_ContactType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_ContactType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_ContactType_SIZE; ++i)
	{
		getName(SLFF_ContactType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_ContactType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_ContactType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_ContactType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_SlipType dat, string& name)
{
	if (dat == sl_slip_off)
	{
		name = "slip_off";
		return;
	}
	if (dat == sl_slip_on)
	{
		name = "slip_on";
		return;
	}
	if (dat == sl_slip_mixed)
	{
		name = "slip_mixed";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_SlipType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_SlipType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_SlipType_SIZE; ++i)
	{
		getName(SLFF_SlipType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_SlipType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_SlipType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_SlipType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_FrictionType dat, string& name)
{
	if (dat == sl_friction_constant)
	{
		name = "friction_constant";
		return;
	}
	if (dat == sl_friction_slip_weakeningLinear)
	{
		name = "friction_slip_weakeningLinear";
		return;
	}
	if (dat == sl_friction_slip_weakeningPower)
	{
		name = "friction_slip_weakeningPower";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_FrictionType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_FrictionType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_FrictionType_SIZE; ++i)
	{
		getName(SLFF_FrictionType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_FrictionType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_FrictionType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_FrictionType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

bool IsSlipWeakening(SLFF_FrictionType tp)
{
	return ((tp == sl_friction_slip_weakeningLinear) || (tp == sl_friction_slip_weakeningPower));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_ContactSeparationRegularizationType dat, string& name)
{
	if (dat == csreg_linear)
	{
		name = "reg_linear";
		return;
	}
	if (dat == csreg_sqrt1)
	{
		name = "reg_sqrt1";
		return;
	}
	if (dat == csreg_sqrt2)
	{
		name = "reg_sqrt2";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_ContactSeparationRegularizationType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_ContactSeparationRegularizationType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_ContactSeparationRegularizationType_SIZE; ++i)
	{
		getName(SLFF_ContactSeparationRegularizationType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_ContactSeparationRegularizationType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_ContactSeparationRegularizationType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_ContactSeparationRegularizationType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_funHType dat, string& name)
{
	if (dat == AllixHExp)
	{
		name = "AllixHExp";
		return;
	}
	name = "none";
	cout << "dat\t" << (int)dat << '\n';
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_funHType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_funHType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_funHType_SIZE; ++i)
	{
		getName(SLFF_funHType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_funHType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_funHType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_funHType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getName(SLFF_RBH_VelDDotType dat, string& name)
{
	if (dat == RBH_VelDDot_none)
	{
		name = "RBH_VelDDot_none";
		return;
	}
	if (dat == RBH_VelDDot_vnp_only)
	{
		name = "RBH_VelDDot_vnp_only";
		return;
	}
	if (dat == RBH_VelDDot_vtabs_only)
	{
		name = "RBH_VelDDot_vtabs_only";
		return;
	}
	if (dat == RBH_VelDDot_vnp_vtabs)
	{
		name = "RBH_VelDDot_vnp_vtabs";
		return;
	}
	if (dat == RBH_VelDDot_vel_effective)
	{
		name = "RBH_VelDDot_vel_effective";
		return;
	}
	name = "none";
	THROW("Invalid value\n");
}

ostream& operator<<(ostream& out, SLFF_RBH_VelDDotType  dat)
{
	string name;
	getName(dat, name);
	out << name;
	return out;
}

bool string2Type(string name, SLFF_RBH_VelDDotType& dat)
{
	string name2;
	for (int i = 0; i < SLFF_RBH_VelDDotType_SIZE; ++i)
	{
		getName(SLFF_RBH_VelDDotType(i), name2);
		if (name2 == name)
		{
			dat = SLFF_RBH_VelDDotType(i);
			return true;
		}
	}
	return false;
}

istream& operator>>(istream& in, SLFF_RBH_VelDDotType& dat)
{
	string bufstr, name;
	READ_NSTRING(in, bufstr, name);
	int tmpi;
	if (fromString(name, tmpi) == true)
	{
		dat = SLFF_RBH_VelDDotType(tmpi);
	}
	if (string2Type(name, dat) == false)
	{
		THROW("invalid input parameter\n");
	}
	return in;
}
