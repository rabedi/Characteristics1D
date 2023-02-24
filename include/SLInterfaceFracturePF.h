#ifndef SL_INTERFACE_FRACTURE_PF__H
#define SL_INTERFACE_FRACTURE_PF__H

#include "globalMacros.h"
#include "globalTypesClasses.h"
#include "SLFracTypes.h"

// these are parameters for global model, e.g. ring problem r, p, ..., error tolerances, ...
class SLFractureGlobal_Configuration
{
public:
	SLFractureGlobal_Configuration();
	void Initialize_SLFractureGlobal_Configuration_After_Reading();
	void Read(string configNameIn);// does not read computed parts
	void Print(ostream& out, bool printComputedVals) const;
	// this function proposes a new time step size, from the current size and maximum error to tolerance ratio
	// acceptFlag
	void Update_AdaptivityFlag_From_current_delT_max_err2tol(double current_delT, double error_2_tol_max, AdaptivityS& as, bool& accept_point);
	// if any time scale entry is < 0, it's scaled relative to the minimum time scale of the domain
	void UpdateTimeScales(double min_domain_del_t, double max_domain_del_t);
	// parameters that can be used in general time marching and computing values of step (n + 1) from step n and possibly earlier
	// in general only one earlier step (i.e. n) is used to form values of step (n + 1)
	// at initial time, etc. only 1 is used
	int max_num_earlier_steps;
	// consider when there is only one earlier value for the ODE
	// fDot = source_f. At time step n + 1, with only one earlier value, we'll update f_n+1 as:
	// -	iteration 0:	alpha = one_step_update_alpha_it0			f_n+1 = f_n + source_f(n) * alpha * del_t
	//		that is, alpha = 0, simply uses earlier value
	double one_step_update_alpha_it0;
	// -	iteration > 0	alpha = one_step_update_alpha_it_gt0		f_n+1 = f_n + (source_f(n) * alpha + (1 - alpha) * source_f(n+1)) * del_t
	double one_step_update_alpha_it_gt0;

	// for solution from step n -> n + 1 there are two sources of error
	// 1.	within iterations from n -> n + 1: continuously previous assumption of traction and velocities are compared to newly derived ones (star values)
	//		it any is above corresponding tolerance, the in-between n -> n + 1 iteration continues
	//		max_n2np1_num_subiterations is the maximum number of in-between iterations n -> n + 1 

	// 2.  from time step n to n + 1 solution. If the change from step n to n + 1 is large, the time step will be adjusted and the step will be retaken
	// tolerances are given below. 
	//					For del_v and del_s tolerances		del_s_type (v) == 
	//									ect_nonDimensional ->  error is nondimensional (e.g. stress/sigmaC and velocity/velC)
	//									ect_Dimensional	   -> absolute value is used
	//									ect_notActive	   ->error is not active
	//					see enumeration errorCheckT

	// 1. iterations for going from time step n -> n + 1
	double within_step_iter_del_s_tol;	
	double within_step_iter_del_v_tol;

	// does not allow indefinite iterations, after reaching the number below one of the two things happen - 
	int	   within_step_iter_max_num_iter;
	//  refine_after_max_num_iter:
	//		true:
	//				in between_steps_adaptivity_s sets			(proposed delT = current_delT * refine_after_max_num_iter_del_t_factor)
	//						refine: if proposed delT is < uniform_del_t
	//						terminate (a_terminate_run_prematurely) otherwise
	//		false:	assumes convergence of iterations and goes to 
	bool refine_after_max_num_iter;
	double uniform_del_t;
	double min_del_t;
	double max_del_t;
	// del_t reduction factor
	double refine_after_max_num_iter_del_t_factor;

	// 2. from time step n -> n + 1, whether time step should be decreased. It also suggests a new time step.

	// this boolean specifies if between step adaptivity check is done or not
	bool   between_steps_adaptivity;
	double between_steps_del_s_tol;
	double between_steps_del_v_tol;
	// damage is nondimensional -> so no normalization is needed
	double between_steps_del_damage_tol;
	// this is a tolerance for a nondimensional quantity as damage_source is already multiplied by tau
	double between_steps_damage_source_tau_tol;
	// change of c from separation to contact from earlier to current step (contact to separation is not included as that one can be rapid)
	double between_steps_del_sep2cont_c;

	errorCheckT del_s_type;
	errorCheckT del_v_type;

	errorCheckT del_s_type_within_step;
	errorCheckT del_v_type_within_step;

	/// values used for adaptivity limit decisions and change of time step
	// 1.2 & 0.8
	double coarsening_error_ratio_lim, refinement_error_ratio_lim;
	// 2.0 & 0.5
	double coarsening_delt_factor, refinement_delt_factor;
	// computed
	double inv_uniform_del_t;

	double 	g_fp_tol;
	// for comparing time values
	double g_fp_tol_time;

	/// printing options
	string rootFolder;
	// these are the booleans that determine if certain things should be output
	// depending on whether an interface is contact-bonded only / BC, dealing with adaptive run, etc. some of these flags are turned to false
	// scalars are contact, ... area fractions calculated for fracture problems
	bool print_Scalars;
	// adaptivity log
	bool print_Adaptivity;
	// iterations for fracture interface
	bool print_outIterationConv;

	// run terminating criteria
	double terminate_run_target_time;	// simulation goes until all spatial points are above this time
	double terminate_run_target_max_damage;	// simulation continues until damage at one point reaches this value
	// others (e.g. contact, etc.) can be added later, especially for contact problems

	AdaptivityF a_flag_default;
private:
	void Read(istream& in);// does not read computed parts
};

extern SLFractureGlobal_Configuration* g_slf_conf;

// this class stores fracture / contact properties that will be used to compute star values
class SL_Interface_Fracture_PF
{
public:
	SL_Interface_Fracture_PF();
	void Read_SL_Interface_Fracture_PF(istream& in, int interfaceFlag);// does not read computed parts
	void Print(ostream& out, bool printComputedVals) const;
	void PrintShort(ostream& out) const;

	void Finalize_Fracture_PF_AfterReading();
	void CreateBonded_PF();
	inline bool HasAnyContactFriction() const {		return (damageOffOnMix != sl_interfacial_damage_off);	}
	// This is similar to the function above in return value, but has a different purpose. 
	// It decides whether any convergence check within_step is needed (from iteration k -> k+1)
	// For purely linear problem, i.e. Bonded solution, no (estimated) traces of velocity, traction, and displacement is used. So, no need to make sure those estimates are accurate.
	inline bool Need_BetweenIteration_ConvergenceCheck() const { return (damageOffOnMix != sl_interfacial_damage_off); }
	// initial damage and source are inputs that can be modified
	void GetIC_ScalarValues(const VEC& sigma, double sigmaCFactor, double deltaCFactor, double& initial_damage, double& maxEffDelU, double& initial_damage_source, double* scalar_Vals) const;
	inline double Get_Const_friction_Coef() const {		return damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax;	}
	// beta values used to compute effective scalar stress and displacement/velocities
	inline bool FullyDebondedInterfaceInTensile() const {	return ((damageOffOnMix == sl_interfacial_damage_on) && ((damageTractionModel == sl_interfacial_damage_traction_zero) || (tsrModel == tsr_Zero)));	}
	// an implicit (backward Euler) solution procedue for Ortiz modle is given in TSR1D.* files. If this boolean is true, such method is going to be use
	inline bool IsSimpleOrtizModel() const { return ((damageOffOnMix == sl_interfacial_damage_on) && (contactOffOnMix == sl_contact_off) && (damageTractionModel == sl_interfacial_damage_traction_TSR) && (tsrModel == tsr_Ortiz)); }
	// initial damage becomes initial damage for damage evolution models and is set to some nonzero max displacement for cohesive zone models
//	void Update_InitialDamage_MaxDelUEffective(double iniDamageIn, double deltaC, double& iniDamageOut, double& maxDelU);
	inline bool TRS_Tensile_Mode() const { return ((damageOffOnMix == sl_interfacial_damage_on) && (damageTractionModel == sl_interfacial_damage_traction_TSR) && (tsrModel != tsr_Zero)); }
	void Get_sigmaC_deltaC_phiC_scales(double& sigmaCScale, double& deltaCScale, double& energyCScale) const;

	double beta_traction;
	double beta_delU;


	////////////////////////////////////////////////////////////
	// Bonded-debonded transition
	SLFF_InterfacialDamageModeType damageOffOnMix;

	////////////////////////////////////////////////////////////
	// contact-separation transition
	SLFF_ContactType			   contactOffOnMix;
	// values for separation contact interpolation (these are normalized w.r.t. deltaC)
	double rel_dcont, rel_dsep, rel_dcont_or_dsep;
	SLFF_ContactSeparationRegularizationType sep2cont_regType;

	////////////////////////////////////////////////////////////
	// stick-slip transition
	SLFF_SlipType				   slipOffOnMix;
	SLFF_FrictionType			   frictionModel;
	// parameters for friction coefficient
	bool isIsoFrictionParameters;

	// slip-weakening upper and lower k's
	// shape of decay parameters
	double damaged_slip_dt, damaged_para0_epsSW, damaged_para1_n;
	// friction coefficient parameters
	// 2D	positive y direction OR 3D for the angle of max value
	double damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax, damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax;
	// 2D	negative y direction OR 2D for the angle of inx value
	double damaged_kUpper_2D_neg_y_dir_or_3D_angularMin, damaged_kLower_2D_neg_y_dir_or_3D_angularMin;
	// 3D
	//	3D: the minimum values is at 180 (opposite direction) apart - kUpper and lower will be angle dependent
	double theta4MaxValDeg;
	double theta4MaxValRad; // angle (in radians) for which k_upper and lower both take their angular max value

	////////////////////////////////////////////////////////////
	// separation mode parameters
	SLFF_InterfacialDamageTractionType damageTractionModel;
	SLFF_TSRType				   tsrModel;

	//////////////////////////// Damage evolution source
	// (*1) Dtarget = (y / ybar - c)/(1 - c), ybar: normalizer, c = brittleness
	//	ybar = ((1 - kn) (1 - D) ^ pn + kn) normalizer0
	// c = ((1 - kc) (1 - D) ^ pc + kc) c0
		
	double c0;
	double tau;
	double Hpara0;

	SLFF_funHType hfunT;

	int pn, pc;
	bool cConst, nConst;
	double factor_k;
	double kn;
	double kc;
	// double maxDamageValue; 	// maximum value damage can take in REGULARIZATION process. This  is commented out as no regularization is done here
	// if srcZeroFldGtmxV == true, for damage > mxV4SrcZero, src becomes zero
	bool srcZeroFldGtmxV;
	double mxV4SrcZero;
	// if damage gets close to 1 within this tolerance, full damage is considered
	// see Finalize_Fracture_PF_AfterReading() on how it's set
	double zeroTol4DamageEq1;
	bool fx_boundedBy1;

	// RBH model
	SLFF_RBH_VelDDotType rbhvDTT;
	// these are the relative displacement scales of the RBH model to displacement scale(s) of Allix model
	double dvnRatio, dvtRatio;

	//////////////////////////// strengths
	// general strength
	double gen_strength;
	// tensile strength (for MC if sn < 0, it's a factor of MC given strength (sigmaC/k)
	double sn;
	// cohesion for MC models
	double c;
	// friction coefficient
//	double k;
	// effective stress model
	SLFF_Eeffective_StsType effType;

	// displacement scales
	double deltaC;


	//////////////////////////// in-situ stresses
	bool has_in_situ;
	VEC in_situ_stress;
};

#endif

