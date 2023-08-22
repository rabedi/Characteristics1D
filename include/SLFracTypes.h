#ifndef SLFRAC_TYPES__H
#define SLFRAC_TYPES__H

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <map>
using namespace std;

// types for fracture mechanics

typedef enum {
	SLFF_rho = 0, SLFF_E = 1, SLFF_nu = 2, SLFF_Cn = 3, SLFF_Ct = 4, SLFF_Zn = 5, SLFF_Zt = 6, SLFF_cn = 7, SLFF_ct = 8, SLFF_cR = 9,
	// strengths s = general strength (sn, st, sc = normal, shear, compressive strengths), c = cohesion, beta = sn / st, ... in square model
	// k for friction is put in the friction category
	// eta = sc / sn
	SLF_s = 10, SLF_sn = 11, SLF_st = 12, SLF_sc = 13, SLF_c = 14, SLF_beta = 15, SLF_eta = 16, SLF_brittleness = 17,
	// fracture toughness, work of seperation, even possibly KIc if needed later
	// YScale: Elastic = sn^2/E/2; phase field = Gn/bn
	SLF_Gn = 20, SLF_Gt = 21, SLF_YScaleElastic = 22, SLF_YScalePF = 23,
	// displacement scales (dn, dt are displacement scales for TSR and other interfacial models), 
	// bn, bt: length scales of phase field or gradient damage models
	// dvn, dvt are displacement scales for RBH model and my displacement-base interfacial damage model
	// slip weakening length scale and contact regularization length scale are put in contact/friction group
	// SLF_lOverbn and , SLF_lOverbt ( = a1blcFactor) are multiplied by a1Factor (often 4/pi) in the form of a1 (a1 = a1blcFactor * a1Factor)
	SLF_dn = 30, SLF_dt = 31, SLF_bn = 32, SLF_bt = 33, SLF_dvn = 34, SLF_dvt = 35, SLF_lOverbn = 36, SLF_lOverbt = 37,

	// time scales
	SLF_taun = 40, SLF_taut = 41, //SLF_tauRateModel = 42,

	// length scales (Process Zone size)
	SLF_ln = 50, SLF_lt = 51,

	// velocity scales
	SLF_vn = 60, SLF_vt = 61,

	// contact/separation: deltas for contact regularization (can be relative to SLF_dn (default) or not)  
	SLF_dcont = 70, SLF_dsep = 71,

	// friction, phi is friction coefficient (SLF_dt_slip slip weakening length scale)
	SLF_k = 80, SLF_phiDegree = 81, SLF_phi = 82, SLF_dt_slip = 83,
	// for slip weakening model two parameters of Dp and epsilon are needed (additional parameters)
	// deltaUEff = (<delu>+^2 + betaDelU^2 * delU_t^2)
	SLF_beta_delU = 84, 
	SLF_damaged_kUpper = 85, SLF_damaged_kLower = 86, SLF_damaged_para0_epsSW = 87, SLF_damaged_para1_n = 88, SLF_damaged_para2 = 87,
	SLF_damaged_c = 90, SLF_damaged_sn = 91,
	// for anisotropic k we can have additional parameters to form anisotropic k, for 2D para1's are values for negative direction
	SLF_damaged_kUpper_para1 = 92, SLF_damaged_kLower_para1 = 93,
	SLF_damaged_kUpper_para2 = 94, SLF_damaged_kLower_para2 = 95,
	SLF_damaged_kUpper_para3 = 96, SLF_damaged_kLower_para3 = 97,

	// bulk damage, phase field parameters
	// instead of D in sigma = (1 - D)C epsilon, ((1 - k)(1 - D) + k) C epsilon is used
	SLF_factor_k = 100,
	// see (*1) below on how ybar and c can change versus d in computing target damage
	SLF_td_kn = 101,
	SLF_td_kc = 102,

	// Some other parameters
	// H(<Dtarget - D>+) needs some parameters, first parameter is below
	SLF_Hpara0_Allix_a = 120,

	// phase-field computed parameters
	// rho D'' + d D' + E nabla D = ... (parameters are with tilde) is the equation for phase field
	SLF_pf_timescale = 130, SLF_pf_rhoTilde = 131, SLF_pf_dTilde = 132, SLF_pf_ETilde = 133, 
	// interfacial models
	// in-situ normal, shear0, and shear2 stresses
	SLF_InSitu_sn = 200, SLF_InSitu_st0 = 201, SLF_InSitu_st1 = 202,

	SLFracPropT_SIZE
} SLFracPropT;

typedef enum
{
	// if a random field value is fed to parameters how strengths, Gs, ... are affected
	SLFF_scalarRandomMultiplierV = 10,
	// what type of effSsective (or effective) is used (MC, sqrt, VM, Tresca, ...)
	// type of bulk damage
	SLFF_isBulkDamage = 11,
	SLFF_BulkDamagePFV = 12,

	SLFF_Eeffective_StsV = 20,
	SLFF_Eeffective_StnV = 21,
	// breakdown of epsilon or sigma in going from eps to sig
	SLFF_EpsEl2SigmaV = 22,
	SLFF_EpsE2YbarV = 23, // Ene, sqrt(Ene), strin_eff, strs_eff for Y
	// how the damage force is defined (energy, energySqrt, stress effective, strain effective)
	SLFF_YbarEnergyV = 24, // Ybar (off, strnspecific, total)
	// w(D) = (1 - D), (1 - D)^2, ...
	SLFF_OmegaDV = 25,
	SLFF_AlphaDV = 26,
	// what TSR is used
	SLFF_TSRV = 27,
	//	NormTanCompV = 28,
	SLFF_PF_velOptionV = 28,
	// what field (phase field or its history) to use in [xi/2 + (1 - xi) field]
	SLFF_PF_fieldOptionSrcLV = 29,

	//	interface values
	SLFF_InterfaceTypeV = 30, // it takes one of these values from BndryFT
	SLFF_InterfaceHasInSituV = 31, // computed from parameters: SLF_InSitu_sn, SLF_InSitu_st0, SLF_InSitu_st1, user does not generally provide it
	// if the interface is undamaged, fully damaged, or partially damaged
	SLFF_ForceU_SlipSlnsV = 32,
	SLFF_InterfacialDamageModeV = 33,
	// what solution to be used for damaged part (zero, TSR, HF)
	SLFF_InterfacialDamageTractionV = 34,
	SLFF_ContactV = 35,
	SLFF_SlipV = 36,
	SLFF_FrictionV = 37,
	SLFF_FrictionPropertiesAreIso = 38,
	SLFF_ContactSeparationRegularizationV = 39,

	// bulk damage for eps -> sigma, we use kfactor (k) as sigma = ((1 - k)(1 - D) + k) C epsilon, that is omega(D) is replaced by (1 - k) omega(D) + k
	// Y = -omega(d)' Ybar for bulk and phase field, if the change above is used, this will change to -(1 - k)omega(D)' Ybar, 
	// it is suggested to keep this off
	SLFF_applyFactor_k_2YReductionB = 41,
	// in the expression above, Y = -omega(d)' Ybar, we ignore -omega(d)' which results in Y = Ybar (as if omega(d) = 1 - d is used for computing Y)
	SLFF_applyYReductionBy_wPrime_dB = 42,
	// in Allix's model damage target can be greater than one, the boolean below makes it below or equal to 1
	SLFF_DamageTarget_boundedBy1B = 43,
	// (*1) Dtarget = (y/ybar - c) / (1 - c), c = ((1 - kc) * (1 - d)^pc + kc) c0 , ybar = ((1 - kn) * (1 - d)^pn + kn) kbar0, 
	// these determine how lower (c0) and upper / lower (ybar) decrease as damage increases, parameters kn, and kc are provided in paramaters
	SLFF_td_pn = 44,
	SLFF_td_pc = 45,
// 	SLFF_IDV = 50
	SLFF_funHV = 51,
	// how to add velocity based DDOt, see SLFF_RBH_VelDDotType
	SLFF_RBH_VelDDotV = 52,
	SLFracFlagT_SIZE
} SLFracFlagT;

// s, G, d, tau stand for strength, work of separation / fracture toughness, displacement, and time scales
// x means that value is multiplied by the scalar
// const means the secondary parameter remains constant
// slff_sclrm_ICDamagex means factor is used for initial field damage
typedef enum
{
	slff_sclrm_none, slff_sclrm_sx_dconst, slff_sclrm_sx_Gconst, slff_sclrm_sx_tauconst,
	slff_sclrm_Gx_sconst, slff_sclrm_Gx_dconst, slff_sclrm_Gx_tauconst, slff_sclrm_ICDamagex,
	SLFF_scalarRandomMultiplierType_SIZE
} SLFF_scalarRandomMultiplierType;

// effSs stands for effective STRESS
// sqrt -> sqrt(<sn>^2 + beta^ st^2)
// mc -> st - k sn + c (mct, mctc add separate tensile and compressive failure too)
// rankin -> <s1>
// vm (von Mises) - J2, plasticity
// vmm (modified von Mises) -> nonsymmetric w.r.t. tensile compressive, Rogge_1983_The use of NONSAP to compare the Von Mises and a modified Von Mises yield criteria.pdf, equation 2
// has the form sqrt(a I1 + b J2)
// Drucker–Prager, a I1 + b sqrt(J2), specicif form depends on eta = sigma_c / sigma_n
// eta = sigma_c / sigma_n

typedef enum
{
	effSs_Off, effSs_sqrt, effSs_mc, effSs_mct, effSs_mctc, effSs_HoekBrown, effSs_Rankin, effSs_Tresca, effSsCompressive, effSs_vm, effSs_vmm, effSs_DP, SLFF_Eeffective_StsType_SIZE
} SLFF_Eeffective_StsType;

// effSn stands for effective STRAIN
// effSn_Mazars = sqrt(sum(<eps_i>+^2)
// effSn_maxEps1: Max{<eps_i>+)
// effSn_vmm_Peerling: function of eta and nu
typedef enum
{
	effSn_Off, effSn_Mazars, effSn_maxEps1, effSn_vmm_Peerling, SLFF_Eeffective_StnType_SIZE
} SLFF_Eeffective_StnType;

// mp -> breakdown to minus, plus type
// hd -> break down to hydraulic plus/minus and deviatoric
// eps or sig before mp, hd means that the breakdown is on eps or sigma
typedef enum
{
	eps2sig_noBreakDown, eps2sig_epsmp, eps2sig_sigmp, eps2sig_epshd, eps2sig_sighd, SLFF_EpsEl2SigmaType_SIZE
} SLFF_EpsEl2SigmaType;

// how the damage force is defined (energy, energySqrt, stress effective, strain effective)
typedef enum
{
	Ybar_energy, Ybar_sqrtEnergy, Ybar_stress, Ybar_strain, SLFF_EpsE2YbarType_SIZE
} SLFF_EpsE2YbarType;

// BDE in the second entry means that the energy is computed from the breakdown of the strain stress
typedef enum
{
	Ybar_energyOff, YbarStrnStrsBDE, Ybar_energy_TotalStress, SLFF_YbarEnergyType_SIZE
} SLFF_YbarEnergyType;

// omega(d) form (1 - d), (1 - d)^2, rational using p, a1, a2, a3 Q(d) = a1 d P(d), P(d) = (1 + a2 d + a2.a3 d^2)
typedef enum
{
	omegad_1md, omegad_1md2, omegad_rational, SLFF_OmegaDType_SIZE
} SLFF_OmegaDType;


// alpha(d) = xi d + (1 - xi) d2
// xi = 0 -> d2 -> AT2
// xi = 1 -> d -> AT1
// xi = 2 -> 2d - d2 -> PF_CZM
typedef enum
{
	alphad_xi_0_AT2 = 0, alphad_xi_1_AT1 = 1, alphad_xi_2_PZ_CZM = 2, SLFF_AlphaDType_SIZE
} SLFF_AlphaDType;

// EFM_2018_1021_Revision 1_V0.pdf 
// Phase field and gradient enhanced damage models for quasi-brittle failure: a numerical comparative study
// Mandal, Nguyen, Heidarpour, Table 1 -> entries linear to concrete are from that table
// concrete: Cornelissen et al. (1986)
// tsr_Zero is used for damage model, etc.
typedef enum
{
	tsr_Zero, tsr_Xu_Needleman, tsr_Ortiz, tsr_linear, tsr_bilinear, tsr_exponential, tsr_hyperbolic, tsr_concrete, SLFF_TSRType_SIZE
} SLFF_TSRType;

typedef enum
{
	pfv_cd, pfv_cr, pfv_cs, SLFF_PF_velOptionType_SIZE
} SLFF_PF_velOptionType;

typedef enum
{
	pfsl_history, pfsl_pf, pfsl_pf_bounded, SLFF_PF_fieldOptionSrcLType_SIZE
} SLFF_PF_fieldOptionSrcLType;

// normal and tangent directions
typedef enum
{
	ntc_n, ntc_t, ntc_c, ntc_none, NormTanCompType_SIZE
} NormTanCompType;

#define NormTanSIZE	2
#define NormTanCompSIZE	3

// off -> interface always in separation mode
// on -> always in contact mode
// mixed -> can switch between separation and contact modes

typedef enum
{
	sl_forceU0_forceSlip0, sl_forceU1_forceSlip0, sl_forceU1_forceSlip1, SLFF_ForceU_SlipSlnsType_SIZE
} SLFF_ForceU_SlipSlnsType;

typedef enum
{
	sl_interfacial_damage_off, sl_interfacial_damage_on, sl_interfacial_damage_mixed, SLFF_InterfacialDamageModeType_SIZE
} SLFF_InterfacialDamageModeType;

typedef enum
{
	sl_interfacial_damage_traction_zero, sl_interfacial_damage_traction_TSR, sl_interfacial_damage_traction_HF, SLFF_InterfacialDamageTractionType_SIZE
} SLFF_InterfacialDamageTractionType;

typedef enum
{
	sl_contact_off, sl_contact_on, sl_contact_mixed, SLFF_ContactType_SIZE
} SLFF_ContactType;

typedef enum
{
	sl_slip_off, sl_slip_on, sl_slip_mixed, SLFF_SlipType_SIZE
} SLFF_SlipType;

// type of the friction model
typedef enum
{
	sl_friction_constant, sl_friction_slip_weakeningLinear, sl_friction_slip_weakeningPower, SLFF_FrictionType_SIZE
} SLFF_FrictionType;

typedef enum
{
	csreg_linear, csreg_sqrt1, csreg_sqrt2, SLFF_ContactSeparationRegularizationType_SIZE
} SLFF_ContactSeparationRegularizationType;


/////////////////////////////////////////////////////
// enumeration for Allix's model

// H is a function that takes the value of 0 at 0 and 1 at 1 (the latter may not always be satisfied). 
// It's also an increasing function
// AllixHExp: (1 - exp(-ax))
typedef enum {AllixHExp, SLFF_funHType_SIZE} SLFF_funHType;

// DDot = DDotAllix(effTraction) + DDotRBH(delv)
// This type specifies DDotRBH(delv)
//										DDotRBH
//		RBH_VelDDot_none:				none
//		RBH_VelDDot_vnp_only			<vn>+/delun
//		RBH_VelDDot_vtabs_only			|vt|/delut
//		RBH_VelDDot_vnp_vtabs			<vn>+/delun + |vt|/delut
//		RBH_VelDDot_vel_effective		sqrt( [<vn>+/delun]^2 + betav^2 * [|vt|/delut]^2	)
typedef enum { RBH_VelDDot_none, RBH_VelDDot_vnp_only, RBH_VelDDot_vtabs_only, RBH_VelDDot_vnp_vtabs, RBH_VelDDot_vel_effective, SLFF_RBH_VelDDotType_SIZE } SLFF_RBH_VelDDotType;

typedef enum {noLengthBulkDamage, AllixTimeDelay, SLFF_BulkDamagePFType_SIZE} SLFF_BulkDamagePFType;

// other fracture types
// these consider the relation between damage and damage target in the form
// f(D_target, D) = <D_target - D>+
// NOT read in initial set of flags
typedef enum
{
	damage_target_zero, damage_target_lesst_damage, damage_target_gtt_damage, DamageDTargetRelT_SIZE
} DamageDTargetRelT;


bool IsMC(SLFF_Eeffective_StsType effS);

/////
void getName(SLFF_scalarRandomMultiplierType dat, string& name);
ostream& operator<<(ostream& out, SLFF_scalarRandomMultiplierType dat);
bool string2Type(string name, SLFF_scalarRandomMultiplierType& dat);
istream& operator>>(istream& in, SLFF_scalarRandomMultiplierType& dat);

/////
void getName(SLFF_BulkDamagePFType dat, string& name);
ostream& operator<<(ostream& out, SLFF_BulkDamagePFType dat);
bool string2Type(string name, SLFF_BulkDamagePFType& dat);
istream& operator>>(istream& in, SLFF_BulkDamagePFType& dat);

void getName(SLFracPropT dat, string& name);
ostream& operator<<(ostream& out, SLFracPropT dat);
bool string2Type(string name, SLFracPropT& dat);
istream& operator>>(istream& in, SLFracPropT& dat);

/////
void getName(SLFracFlagT dat, string& name);
ostream& operator<<(ostream& out, SLFracFlagT dat);
bool string2Type(string name, SLFracFlagT& dat);
istream& operator>>(istream& in, SLFracFlagT& dat);

/////
void getName(SLFF_Eeffective_StsType dat, string& name);
ostream& operator<<(ostream& out, SLFF_Eeffective_StsType dat);
bool string2Type(string name, SLFF_Eeffective_StsType& dat);
istream& operator>>(istream& in, SLFF_Eeffective_StsType& dat);

/////
void getName(SLFF_Eeffective_StnType dat, string& name);
ostream& operator<<(ostream& out, SLFF_Eeffective_StnType dat);
bool string2Type(string name, SLFF_Eeffective_StnType& dat);
istream& operator>>(istream& in, SLFF_Eeffective_StnType& dat);

/////
void getName(SLFF_EpsEl2SigmaType dat, string& name);
ostream& operator<<(ostream& out, SLFF_EpsEl2SigmaType dat);
bool string2Type(string name, SLFF_EpsEl2SigmaType& dat);
istream& operator>>(istream& in, SLFF_EpsEl2SigmaType& dat);

/////
void getName(SLFF_EpsE2YbarType dat, string& name);
ostream& operator<<(ostream& out, SLFF_EpsE2YbarType dat);
bool string2Type(string name, SLFF_EpsE2YbarType& dat);
istream& operator>>(istream& in, SLFF_EpsE2YbarType& dat);

/////
void getName(SLFF_YbarEnergyType dat, string& name);
ostream& operator<<(ostream& out, SLFF_YbarEnergyType dat);
bool string2Type(string name, SLFF_YbarEnergyType& dat);
istream& operator>>(istream& in, SLFF_YbarEnergyType& dat);

/////
void getName(SLFF_OmegaDType dat, string& name);
ostream& operator<<(ostream& out, SLFF_OmegaDType dat);
bool string2Type(string name, SLFF_OmegaDType& dat);
istream& operator>>(istream& in, SLFF_OmegaDType& dat);

/////
void getName(SLFF_AlphaDType dat, string& name);
ostream& operator<<(ostream& out, SLFF_AlphaDType dat);
bool string2Type(string name, SLFF_AlphaDType& dat);
istream& operator>>(istream& in, SLFF_AlphaDType& dat);

/////
void getName(SLFF_TSRType dat, string& name);
ostream& operator<<(ostream& out, SLFF_TSRType dat);
bool string2Type(string name, SLFF_TSRType& dat);
istream& operator>>(istream& in, SLFF_TSRType& dat);
bool IsExtrinsic(SLFF_TSRType dat);
// energy = returnValue * deltaC * sigmaC
double GetEnergyConstantFactor(SLFF_TSRType dat);
/////
void getName(SLFF_PF_velOptionType dat, string& name);
ostream& operator<<(ostream& out, SLFF_PF_velOptionType dat);
bool string2Type(string name, SLFF_PF_velOptionType& dat);
istream& operator>>(istream& in, SLFF_PF_velOptionType& dat);

/////
void getName(SLFF_PF_fieldOptionSrcLType dat, string& name);
ostream& operator<<(ostream& out, SLFF_PF_fieldOptionSrcLType dat);
bool string2Type(string name, SLFF_PF_fieldOptionSrcLType& dat);
istream& operator>>(istream& in, SLFF_PF_fieldOptionSrcLType& dat);

/////
void getName(NormTanCompType dat, string& name);
ostream& operator<<(ostream& out, NormTanCompType dat);
bool string2Type(string name, NormTanCompType& dat);
istream& operator>>(istream& in, NormTanCompType& dat);


/////
void getName(SLFF_ForceU_SlipSlnsType dat, string& name);
ostream& operator<<(ostream& out, SLFF_ForceU_SlipSlnsType dat);
bool string2Type(string name, SLFF_ForceU_SlipSlnsType& dat);
istream& operator>>(istream& in, SLFF_ForceU_SlipSlnsType& dat);

/////
void getName(SLFF_InterfacialDamageModeType dat, string& name);
ostream& operator<<(ostream& out, SLFF_InterfacialDamageModeType dat);
bool string2Type(string name, SLFF_InterfacialDamageModeType& dat);
istream& operator>>(istream& in, SLFF_InterfacialDamageModeType& dat);

/////
void getName(SLFF_InterfacialDamageTractionType dat, string& name);
ostream& operator<<(ostream& out, SLFF_InterfacialDamageTractionType dat);
bool string2Type(string name, SLFF_InterfacialDamageTractionType& dat);
istream& operator>>(istream& in, SLFF_InterfacialDamageTractionType& dat);

/////
void getName(SLFF_ContactType dat, string& name);
ostream& operator<<(ostream& out, SLFF_ContactType dat);
bool string2Type(string name, SLFF_ContactType& dat);
istream& operator>>(istream& in, SLFF_ContactType& dat);

/////
void getName(SLFF_SlipType dat, string& name);
ostream& operator<<(ostream& out, SLFF_SlipType dat);
bool string2Type(string name, SLFF_SlipType& dat);
istream& operator>>(istream& in, SLFF_SlipType& dat);

/////
void getName(SLFF_FrictionType dat, string& name);
ostream& operator<<(ostream& out, SLFF_FrictionType dat);
bool string2Type(string name, SLFF_FrictionType& dat);
istream& operator>>(istream& in, SLFF_FrictionType& dat);
bool IsSlipWeakening(SLFF_FrictionType tp);

void getName(SLFF_ContactSeparationRegularizationType dat, string& name);
ostream& operator<<(ostream& out, SLFF_ContactSeparationRegularizationType dat);
bool string2Type(string name, SLFF_ContactSeparationRegularizationType& dat);
istream& operator>>(istream& in, SLFF_ContactSeparationRegularizationType& dat);

/////
void getName(SLFF_funHType dat, string& name);
ostream& operator<<(ostream& out, SLFF_funHType dat);
bool string2Type(string name, SLFF_funHType& dat);
istream& operator>>(istream& in, SLFF_funHType& dat);

/////
void getName(SLFF_RBH_VelDDotType dat, string& name);
ostream& operator<<(ostream& out, SLFF_RBH_VelDDotType dat);
bool string2Type(string name, SLFF_RBH_VelDDotType& dat);
istream& operator>>(istream& in, SLFF_RBH_VelDDotType& dat);


#endif