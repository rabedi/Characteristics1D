#include "SLInterfaceFracturePF.h"
#include "LAfuncsFinalStep.h"
#include "SL_Interface_PtData.h"
#include "SLInterfaceCalculator.h"

extern SLFractureGlobal_Configuration* g_slf_conf = NULL;

SLFractureGlobal_Configuration::SLFractureGlobal_Configuration()
{
	max_num_earlier_steps = 1;
	one_step_update_alpha_it0 = 1.0;
	one_step_update_alpha_it_gt0 = 0.5;

	within_step_iter_del_s_tol = 0.01;
	within_step_iter_del_v_tol = 0.01;
	within_step_iter_max_num_iter = 10;
	refine_after_max_num_iter = false;
	uniform_del_t = 0.001;
	min_del_t = 0.00001;
	max_del_t = 10.0;
	inv_uniform_del_t = 1.0 / uniform_del_t;
	refine_after_max_num_iter_del_t_factor = 0.5;

	between_steps_adaptivity = true;
	between_steps_del_s_tol = 0.04;
	between_steps_del_v_tol = 0.04;
	between_steps_del_damage_tol = 0.1;
	between_steps_damage_source_tau_tol = 0.1;
	between_steps_del_sep2cont_c = -1.0; // inactive

	del_s_type = ect_nonDimensional;
	del_v_type = ect_nonDimensional;

	del_s_type_within_step = ect_nonDimensional;
	del_v_type_within_step = ect_nonDimensional;

	coarsening_error_ratio_lim = 1.2, refinement_error_ratio_lim = 0.8;
	coarsening_delt_factor = 2.0, refinement_delt_factor = 0.5;

	g_fp_tol = 1e-5;

	print_Scalars = true;
	print_Adaptivity = true;
	print_outIterationConv = true;
	rootFolder = "";

	terminate_run_target_time = 100.0;
	terminate_run_target_max_damage = 1.0 - 1e-7;

	a_flag_default = a_none;
}

void SLFractureGlobal_Configuration::Initialize_SLFractureGlobal_Configuration_After_Reading() 
{
	if (between_steps_del_s_tol < 0)
		del_s_type = ect_notActive;
	if (between_steps_del_v_tol < 0)
		del_v_type = ect_notActive;

	if (within_step_iter_del_s_tol < 0)
		del_s_type_within_step = ect_notActive;
	if (within_step_iter_del_v_tol < 0)
		del_v_type_within_step = ect_notActive;

	inv_uniform_del_t = 1.0 / uniform_del_t;
	g_fp_tol_time = g_fp_tol * uniform_del_t;

	a_flag_default = a_coarsen;
	if (!between_steps_adaptivity)
	{
		print_Adaptivity = false;
		a_flag_default = a_none;
	}

	if (rootFolder.size() > 0)
	{
		do_mkdir(rootFolder.c_str());
		rootFolder += "/";
	}
}

void SLFractureGlobal_Configuration::Read(istream& in)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
		else
		{
			while ((buf != "infile_ts_adaptivity") && (!in.eof()))
				READ_NSTRING(in, buf, buf);
			if (in.eof())
				THROW("Reached end of file looking for infile_ts_adaptivity\n");
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
		if (buf == "max_num_earlier_steps")
		{
			READ_NINTEGER(in, buf, max_num_earlier_steps);
		}
		else if (buf == "one_step_update_alpha_it0")
		{
			READ_NDOUBLE(in, buf, one_step_update_alpha_it0);
		}
		else if (buf == "one_step_update_alpha_it_gt0")
		{
			READ_NDOUBLE(in, buf, one_step_update_alpha_it_gt0);
		}
		else if (buf == "within_step_iter_del_s_tol")
		{
			READ_NDOUBLE(in, buf, within_step_iter_del_s_tol);
		}
		else if (buf == "within_step_iter_del_v_tol")
		{
			READ_NDOUBLE(in, buf, within_step_iter_del_v_tol);
		}
		else if (buf == "within_step_iter_max_num_iter")
		{
			READ_NINTEGER(in, buf, within_step_iter_max_num_iter);
		}
		else if (buf == "refine_after_max_num_iter")
		{
			READ_NBOOL(in, buf, refine_after_max_num_iter);
		}
		else if (buf == "del_s_type_within_step")
			in >> del_s_type_within_step;
		else if (buf == "del_v_type_within_step")
			in >> del_v_type_within_step;
		else if (buf == "uniform_del_t")
		{
			READ_NDOUBLE(in, buf, uniform_del_t);
			if (sfcm.success && (sfcm.cfl_factor > 0))
				uniform_del_t = -sfcm.cfl_factor;

			string key = "uniform_del_t";	map<string, string>* mpPtr;
			double value = -1;
			bool found = Find_Version_Value(key, value, mpPtr);
			if (!found)
			{
				string key = "cfl";
				found = Find_Version_Value(key, value, mpPtr);
				value = -value;
			}
			if (found)
				uniform_del_t = value;
			// whether to keep delt fixed for all spatial mesh resolutions or not?
			key = "dt_resap";
			double resolutionPower;
			found = Find_Version_Value(key, resolutionPower, mpPtr);
			if (found && (fabs(resolutionPower) > 1e-4))
			{
				key = "resolutionFactor"; // see if the resolution factor is something other than 0 and 1
				double resolutionFactor;
				found = Find_Version_Value(key, resolutionFactor, mpPtr);
				if ((found) && ((resolutionFactor - 1.0) > 1e-2)) // means mesh is coarsened
					uniform_del_t /= pow(resolutionFactor, resolutionPower);
			}

			
			g_logout << "\tuniform_del_t\t" << uniform_del_t;
		}
		else if (buf == "max_del_t")
		{
			READ_NDOUBLE(in, buf, max_del_t);
		}
		else if (buf == "min_del_t")
		{
			READ_NDOUBLE(in, buf, min_del_t);
		}
		else if (buf == "refine_after_max_num_iter_del_t_factor")
		{
			READ_NDOUBLE(in, buf, refine_after_max_num_iter_del_t_factor);
		}
		else if (buf == "between_steps_adaptivity")
		{
			READ_NBOOL(in, buf, between_steps_adaptivity);
		}
		else if (buf == "between_steps_del_s_tol")
		{
			READ_NDOUBLE(in, buf, between_steps_del_s_tol);
		}
		else if (buf == "between_steps_del_v_tol")
		{
			READ_NDOUBLE(in, buf, between_steps_del_v_tol);
		}
		else if (buf == "between_steps_del_damage_tol")
		{
			READ_NDOUBLE(in, buf, between_steps_del_damage_tol);
		}
		else if (buf == "between_steps_damage_source_tau_tol")
		{
			READ_NDOUBLE(in, buf, between_steps_damage_source_tau_tol);
		}
		else if (buf == "between_steps_del_sep2cont_c")
		{
			READ_NDOUBLE(in, buf, between_steps_del_sep2cont_c);
		}
		else if (buf == "del_s_type")
			in >> del_s_type;
		else if (buf == "del_v_type")
			in >> del_v_type;
		else if (buf == "coarsening_error_ratio_lim")
		{
			READ_NDOUBLE(in, buf, coarsening_error_ratio_lim);
		}
		else if (buf == "refinement_error_ratio_lim")
		{
			READ_NDOUBLE(in, buf, refinement_error_ratio_lim);
		}
		else if (buf == "coarsening_delt_factor")
		{
			READ_NDOUBLE(in, buf, coarsening_delt_factor);
		}
		else if (buf == "refinement_delt_factor")
		{
			READ_NDOUBLE(in, buf, refinement_delt_factor);
		}
		else if (buf == "print_Scalars")
		{
			READ_NBOOL(in, buf, print_Scalars);
		}
		else if (buf == "print_Adaptivity")
		{
			READ_NBOOL(in, buf, print_Adaptivity);
		}
		else if (buf == "print_outIterationConv")
		{
			READ_NBOOL(in, buf, print_outIterationConv);
		}
		else if (buf == "terminate_run_target_time")
		{
			READ_NDOUBLE(in, buf, terminate_run_target_time);
			if (sfcm.success && (sfcm.tFinal > 0))
				terminate_run_target_time = sfcm.tFinal;

			string key = "tF";	map<string, string>* mpPtr;
			double value = -1;
			bool found = Find_Version_Value(key, value, mpPtr);
			if (found)
				terminate_run_target_time = value;

			// Final time already adjusted in Domain_AllInterfacesAllTimes
			//key = "delc_Tf_fact";
			//if (Find_Version_Value(key, value, mpPtr) == true)
			//	terminate_run_target_time *= value;

			g_logout << "\tterminate_run_target_time\t" << terminate_run_target_time;
		}
		else if (buf == "terminate_run_target_max_damage")
		{
			READ_NDOUBLE(in, buf, terminate_run_target_max_damage);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Initialize_SLFractureGlobal_Configuration_After_Reading();
}

void SLFractureGlobal_Configuration::Read(string configNameIn)
{
	fstream in(configNameIn.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << "configNameIn\t" << configNameIn << '\n';
		THROW("Cannot open file\n");
	}
	Read(in);
}

void SLFractureGlobal_Configuration::Print(ostream& out, bool printComputedVals) const
{
	out << "{\n";
	out << "max_num_earlier_steps\t" << max_num_earlier_steps << '\n';
	out << "one_step_update_alpha_it0\t" << one_step_update_alpha_it0 << '\n';
	out << "one_step_update_alpha_it_gt0\t" << one_step_update_alpha_it_gt0 << '\n';

	out << "within_step_iter_del_s_tol\t" << within_step_iter_del_s_tol << '\n';
	out << "within_step_iter_del_v_tol\t" << within_step_iter_del_v_tol << '\n';
	out << "within_step_iter_max_num_iter\t" << within_step_iter_max_num_iter << '\n';

	out << "refine_after_max_num_iter\t" << refine_after_max_num_iter << '\n';
	out << "uniform_del_t\t" << uniform_del_t << '\n';
	out << "refine_after_max_num_iter_del_t_factor\t" << refine_after_max_num_iter_del_t_factor << '\n';

	out << "between_steps_adaptivity\t" << between_steps_adaptivity << '\n';
	out << "between_steps_del_s_tol\t" << between_steps_del_s_tol << '\n';
	out << "between_steps_del_v_tol\t" << between_steps_del_v_tol << '\n';
	out << "between_steps_del_damage_tol\t" << between_steps_del_damage_tol << '\n';
	out << "between_steps_damage_source_tau_tol\t" << between_steps_damage_source_tau_tol << '\n';
	
	out << "del_s_type\t" << del_s_type << '\n';
	out << "del_v_type\t" << del_v_type << '\n';

	out << "coarsening_error_ratio_lim\t" << coarsening_error_ratio_lim << '\n';
	out << "refinement_error_ratio_lim\t" << refinement_error_ratio_lim << '\n';

	out << "coarsening_delt_factor\t" << coarsening_delt_factor << '\n';
	out << "refinement_delt_factor\t" << refinement_delt_factor << '\n';

	out << "}\n";
	if (!printComputedVals)
		return;
//	out << "[\n";
//	out << "]\n";
}

void SLFractureGlobal_Configuration::Update_AdaptivityFlag_From_current_delT_max_err2tol(double current_delT, double error_2_tol_max, AdaptivityS& as, bool& accept_point)
{
	double inv_error_2_tol_max = 1.0 / (error_2_tol_max + 1e-50);
	if (inv_error_2_tol_max > coarsening_error_ratio_lim) // coarsening range
	{
		accept_point = true;
		as.Update_AdaptFlag(a_coarsen);
		as.a_delt = MIN(coarsening_delt_factor * current_delT, max_del_t);
	}
	else if (inv_error_2_tol_max < refinement_error_ratio_lim) // coarsening range
	{
		accept_point = false;
		as.a_delt = refinement_delt_factor * current_delT;
		if (as.a_delt < min_del_t)
			as.Update_AdaptFlag(a_terminate_run_prematurely);
		else
			as.Update_AdaptFlag(a_refine);
	}
	else
	{
		accept_point = true;
		as.Update_AdaptFlag(a_none);
		as.a_delt = current_delT;
	}
}

void SLFractureGlobal_Configuration::UpdateTimeScales(double min_domain_del_t, double max_domain_del_t)
{
	CFL_provided = (uniform_del_t < 0);
	if (CFL_provided)
		uniform_del_t *= -min_domain_del_t;
	else
	{
		if (uniform_del_t > min_domain_del_t)
		{
			cout << "uniform_del_t\t" << uniform_del_t << '\n';
			cout << "min_domain_del_t\t" << min_domain_del_t << '\n';
			THROW("reduce uniform_del_t in the config file to at least min_domain_del_t for CFL = 1\n");
		}
	}
	if (min_del_t < 0)
		min_del_t *= -min_domain_del_t;
	if (max_del_t < 0)
		max_del_t *= -min_domain_del_t;

	inv_uniform_del_t = 1.0 / uniform_del_t;
	if (terminate_run_target_time < 0)
		terminate_run_target_time *= -max_del_t;
	Initialize_SLFractureGlobal_Configuration_After_Reading();
	if (!between_steps_adaptivity)
	{
		long long maxNumStepsNeeded = (long long)ceil(max_domain_del_t / uniform_del_t);
		if ((maxNumStepsNeeded + 2) > SIZE_PT_TIME_SEQUENCE)
		{
			cout << "min_domain_del_t\t" << min_domain_del_t << '\n';
			cout << "max_domain_del_t\t" << max_domain_del_t << '\n';
			cout << "In file globalMacros.h, #define SIZE_PT_TIME_SEQUENCE must change\n";
			cout << "current value SIZE_PT_TIME_SEQUENCE = " << SIZE_PT_TIME_SEQUENCE << '\n';
			cout << "change it at least to maxNumStepsNeeded + 2 = " << (maxNumStepsNeeded + 2) << '\n';
			THROW("too small of a time history storage\n");
		}
	}
}

SL_Interface_Fracture_PF::SL_Interface_Fracture_PF()
{
	beta_traction = 1.0;
	beta_delU = 1.0;
	deltaC = -1.0;

	damageOffOnMix = sl_interfacial_damage_mixed;

	contactOffOnMix = sl_contact_mixed;

	rel_dcont_or_dsep = -0.01;
	rel_dcont = -0.03;
	rel_dsep = 0.0;
	sep2cont_regType = csreg_sqrt1;


	slipOffOnMix = sl_slip_mixed;
	frictionModel = sl_friction_constant; // sl_friction_slip_weakeningLinear;

	isIsoFrictionParameters = true;

	damaged_slip_dt = 0.0;
	damaged_para0_epsSW = 0.0;
	damaged_para1_n = 0.0;

	damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax = 0.3, damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax = 0.2;
	damaged_kUpper_2D_neg_y_dir_or_3D_angularMin = 0.4, damaged_kLower_2D_neg_y_dir_or_3D_angularMin = 0.35;

	theta4MaxValDeg = 30.0;

	damageTractionModel = sl_interfacial_damage_traction_zero; //sl_interfacial_damage_traction_TSR;
	tsrModel = tsr_Xu_Needleman; // tsr_Zero, tsr_Xu_Needleman, tsr_Ortiz, tsr_linear;

	sn = 1.0;
	c = 1.0;
//	k = 0.3;
	gen_strength = -1.0;

	c0 = 0.7;
	tau = 1.0;
	Hpara0 = 10.0;

	hfunT = AllixHExp;
	pn = 0, pc = 0;
	factor_k = 0.0; // 0.05 
	kn = 0.0; //  0.05;
	kc = 0.0;
	fx_boundedBy1 = true;

	// RBH
	rbhvDTT = RBH_VelDDot_none;
	dvnRatio = 10.0, dvtRatio = 10.0;

//	maxDamageValue = 1.0; // 0.95
	srcZeroFldGtmxV = true;
	mxV4SrcZero = 1.0;
	zeroTol4DamageEq1 = -1e-2;

	// in-situ
	has_in_situ = false;
	setValue(in_situ_stress, 0.0);
}

void SL_Interface_Fracture_PF::Read_SL_Interface_Fracture_PF(istream& in, int interfaceFlag)
{
	bool versionChange = ((sfcm.success) && (sfcm_gen.bulkFlag == interfaceFlag));
	double value;
	string key;	map<string, string>* mpPtr;

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
		if (buf == "beta_traction")
		{
			READ_NDOUBLE(in, buf, beta_traction);
		}
		else if (buf == "beta_delU")
		{
			READ_NDOUBLE(in, buf, beta_delU);
		}
		else if (buf == "damageOffOnMix")
			in >> damageOffOnMix;
		else if (buf == "contactOffOnMix")
			in >> contactOffOnMix;
		else if (buf == "rel_dcont")
		{
			READ_NDOUBLE(in, buf, rel_dcont);
		}
		else if (buf == "rel_dsep")
		{
			READ_NDOUBLE(in, buf, rel_dsep);
		}
		else if (buf == "rel_dcont_or_dsep")
		{
			READ_NDOUBLE(in, buf, rel_dcont_or_dsep);
		}
		else if (buf == "sep2cont_regType")
			in >> sep2cont_regType;
		else if (buf == "slipOffOnMix")
			in >> slipOffOnMix;
		else if (buf == "frictionModel")
			in >> frictionModel;
		else if (buf == "isIsoFrictionParameters")
		{
			READ_NBOOL(in, buf, isIsoFrictionParameters);
		}
		else if (buf == "damaged_slip_dt")
		{
			READ_NDOUBLE(in, buf, damaged_slip_dt);
		}
		else if (buf == "damaged_para0_epsSW")
		{
			READ_NDOUBLE(in, buf, damaged_para0_epsSW);
		}
		else if (buf == "damaged_para1_n")
		{
			READ_NDOUBLE(in, buf, damaged_para1_n);
		}
		else if (buf == "damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax")
		{
			READ_NDOUBLE(in, buf, damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax);
		}
		else if (buf == "damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax")
		{
			READ_NDOUBLE(in, buf, damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax);
		}
		else if (buf == "damaged_kUpper_2D_neg_y_dir_or_3D_angularMin")
		{
			READ_NDOUBLE(in, buf, damaged_kUpper_2D_neg_y_dir_or_3D_angularMin);
		}
		else if (buf == "damaged_kLower_2D_neg_y_dir_or_3D_angularMin")
		{
			READ_NDOUBLE(in, buf, damaged_kLower_2D_neg_y_dir_or_3D_angularMin);
		}
		else if (buf == "theta4MaxValDeg")
		{
			READ_NDOUBLE(in, buf, theta4MaxValDeg);
		}
		else if (buf == "damageTractionModel")
			in >> damageTractionModel;
		else if (buf == "tsrModel")
			in >> tsrModel;
		else if (buf == "c0")
		{
			READ_NDOUBLE(in, buf, c0);
			key = "c0";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
				c0 = value;
		}
		else if (buf == "tau")
		{
			READ_NDOUBLE(in, buf, tau);
		}
		else if (buf == "Hpara0")
		{
			READ_NDOUBLE(in, buf, Hpara0);
		}
		else if (buf == "hfunT")
			in >> hfunT;
		else if (buf == "pn")
		{
			READ_NINTEGER(in, buf, pn);
		}
		else if (buf == "pc")
		{
			READ_NINTEGER(in, buf, pc);
		}
		else if (buf == "factor_k")
		{
			READ_NDOUBLE(in, buf, factor_k);
		}
		else if (buf == "kn")
		{
			READ_NDOUBLE(in, buf, kn);
		}
		else if (buf == "kc")
		{
			READ_NDOUBLE(in, buf, kc);
		}
		else if (buf == "srcZeroFldGtmxV")
		{
			READ_NBOOL(in, buf, srcZeroFldGtmxV);
		}
		else if (buf == "mxV4SrcZero")
		{
			READ_NDOUBLE(in, buf, mxV4SrcZero);
		}
		else if (buf == "zeroTol4DamageEq1")
		{
			READ_NDOUBLE(in, buf, zeroTol4DamageEq1);
		}
		else if (buf == "fx_boundedBy1")
		{
			READ_NBOOL(in, buf, fx_boundedBy1);
		}
		else if (buf == "rbhvDTT")
			in >> rbhvDTT;
		else if (buf == "dvnRatio")
		{
			READ_NDOUBLE(in, buf, dvnRatio);
			key = "dvnRatio";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
				dvnRatio = value;
		}
		else if (buf == "dvtRatio")
		{
			READ_NDOUBLE(in, buf, dvtRatio);
		}
		else if (buf == "gen_strength")
		{
			READ_NDOUBLE(in, buf, gen_strength);
			key = "sigc";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
				gen_strength = value;
		}
		else if (buf == "sn")
		{
			READ_NDOUBLE(in, buf, sn);
		}
		else if (buf == "c")
		{
			READ_NDOUBLE(in, buf, c);
		}
//		else if (buf == "k")
//		{
//			READ_NDOUBLE(in, buf, k);
//		}
		else if (buf == "effType")
			in >> effType;
		else if (buf == "deltaC")
		{
			READ_NDOUBLE(in, buf, deltaC);
			if (versionChange)
			{
				key = "ldelc";
				if (Find_Version_Value(key, value, mpPtr) == true)
					deltaC = pow(10.0, value);
				else
				{
					key = "delc";
					if (Find_Version_Value(key, value, mpPtr) == true)
						deltaC = value;
				}
				key = "delc_Tf_fact";
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					deltaC *= value;
					deltaC /= sfcm.sigmaCFactor;
				}
			}
			g_logout << "\tdeltaC\t" << deltaC;


			g_logout << "\tdeltaC\t" << deltaC;
		}
		else if (buf == "has_in_situ")
		{
			READ_NBOOL(in, buf, has_in_situ);
		}
		else if (buf == "in_situ_stress")
		{
			ReadV(in_situ_stress, in);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Finalize_Fracture_PF_AfterReading();
}

void SL_Interface_Fracture_PF::Print(ostream& out, bool printComputedVals) const
{
	out << "{\n";
	out << "beta_traction\t" << beta_traction << '\n';
	out << "beta_delU\t" << beta_delU << '\n';

	out << "damageOffOnMix\t" << damageOffOnMix << '\n';
	out << "contactOffOnMix\t" << contactOffOnMix << '\n';
	out << "rel_dcont\t" << rel_dcont << '\n';
	out << "rel_dsep\t" << rel_dsep << '\n';

	out << "sep2cont_regType\t" << sep2cont_regType << '\n';

	out << "slipOffOnMix\t" << slipOffOnMix << '\n';
	out << "frictionModel\t" << frictionModel << '\n';
	out << "isIsoFrictionParameters\t" << isIsoFrictionParameters << '\n';

	out << "damaged_slip_dt\t" << damaged_slip_dt << '\n';
	out << "damaged_para0_epsSW\t" << damaged_para0_epsSW << '\n';
	out << "damaged_para1_n\t" << damaged_para1_n << '\n';

	out << "damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax\t" << damaged_kUpper_iso_2D_pos_y_dir_or_3D_angularMax << '\n';
	out << "damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax\t" << damaged_kLower_iso_2D_pos_y_dir_or_3D_angularMax << '\n';
	out << "damaged_kUpper_2D_neg_y_dir_or_3D_angularMin\t" << damaged_kUpper_2D_neg_y_dir_or_3D_angularMin << '\n';
	out << "damaged_kLower_2D_neg_y_dir_or_3D_angularMin\t" << damaged_kLower_2D_neg_y_dir_or_3D_angularMin << '\n';

	out << "theta4MaxValDeg\t" << theta4MaxValDeg << '\n';

	out << "damageTractionModel\t" << damageTractionModel << '\n';
	out << "tsrModel\t" << tsrModel << '\n';

	out << "c0\t" << c0 << '\n';
	out << "tau\t" << tau << '\n';
	out << "Hpara0\t" << Hpara0 << '\n';


	out << "hfunT\t" << hfunT << '\n';
	out << "pn\t" << pn << '\n';
	out << "pc\t" << pc << '\n';
	out << "factor_k\t" << factor_k << '\n';
	out << "kn\t" << kn << '\n';
	out << "kc\t" << kc << '\n';

	out << "srcZeroFldGtmxV\t" << srcZeroFldGtmxV << '\n';
	out << "mxV4SrcZero\t" << mxV4SrcZero << '\n';

	out << "zeroTol4DamageEq1\t" << zeroTol4DamageEq1 << '\n';
	out << "fx_boundedBy1\t" << fx_boundedBy1 << '\n';

	out << "rbhvDTT\t" << rbhvDTT << '\n';
	out << "dvnRatio\t" << dvnRatio << '\n';
	out << "dvtRatio\t" << dvtRatio << '\n';

	out << "gen_strength\t" << gen_strength << '\n';
	out << "sn\t" << sn << '\n';
	out << "c\t" << c << '\n';
//	out << "k\t" << k << '\n';
	out << "effType\t" << effType << '\n';
	out << "deltaC\t" << deltaC << '\n';
	out << "has_in_situ\t" << has_in_situ << '\n';
	out << "in_situ_stress\n"; PrintV(in_situ_stress, out);	out << '\n';
	out << "}\n";
	if (!printComputedVals)
		return;
	out << "[\n";
	out << "theta4MaxValRad\t" << theta4MaxValRad << '\n';
	out << "cConst\t" << cConst << '\n';
	out << "nConst\t" << nConst << '\n';
	out << "]\n";
}

void SL_Interface_Fracture_PF::PrintShort(ostream& out) const
{
	out << "damageOffOnMix\t" << damageOffOnMix << '\t';
	out << "contactOffOnMix\t" << contactOffOnMix << '\t';
	out << "slipOffOnMix\t" << slipOffOnMix << '\t';
	out << "frictionModel\t" << frictionModel << '\t';
	out << "gen_strength\t" << gen_strength << '\t';
	out << "sn\t" << sn << '\t';
	out << "deltaC\t" << deltaC << '\t';
	out << "tsrModel\t" << tsrModel << '\n';
}

void SL_Interface_Fracture_PF::CreateBonded_PF()
{
	damageOffOnMix = sl_interfacial_damage_off;
	contactOffOnMix = sl_contact_on;
	slipOffOnMix = sl_slip_off;
	Finalize_Fracture_PF_AfterReading();
}

void SL_Interface_Fracture_PF::Finalize_Fracture_PF_AfterReading()
{
	if (TRS_Tensile_Mode())
		c0 = 1.0;
	theta4MaxValRad = PI / 180.0 * theta4MaxValDeg;
	if (gen_strength < 0)
		gen_strength = sn;

	if (!srcZeroFldGtmxV)
		mxV4SrcZero = 1000.0;

	if (rel_dcont_or_dsep > 0)
	{
		// typical contact problems, where the interface is under compressive stress and should be treated to be in contact (that's why delu = 0 -> contact = 1)
		if (FullyDebondedInterfaceInTensile())
		{
			rel_dcont = 0.0;
			rel_dsep = rel_dcont_or_dsep;
		}
		else // other problems for contact mode some interpenetration happens
		{
			rel_dcont = -rel_dcont_or_dsep;
			rel_dsep = 0.0;
		}
	}

	cConst = (pc == 0);
	nConst = (pn == 0);

#if DiM1
	if (rbhvDTT != RBH_VelDDot_none)
		rbhvDTT = RBH_VelDDot_vnp_only;
	isIsoFrictionParameters = true;
	slipOffOnMix = sl_slip_off;
#else
#if USE_ISO_ASSUMPTION
	isIsoFrictionParameters = true;
#endif
#endif
	if (zeroTol4DamageEq1 < 0.0)
	{
		if (tau > 0.0)
			zeroTol4DamageEq1 *= 1.0 / tau;
		else
			zeroTol4DamageEq1 *= -1.0;
	}
}

void SL_Interface_Fracture_PF::GetIC_ScalarValues(const VEC& sigma, double sigmaCFactor, double deltaCFactor, double& initial_damage, double& maxEffDelU, double& initial_damage_source, double* scalar_Vals) const
{
	for (int i = 0; i < s1i_SIZE; ++i)
		scalar_Vals[i] = 0.0;

	scalar_Vals[s1i_fricCoef] = Get_Const_friction_Coef();

	double sigmaCModified = gen_strength * sigmaCFactor;
	Vec_nt sigma_I_nt_parts;
	sigma_I_nt_parts.Compute_derivedValues_From_vec(sigma, beta_traction);
	double strengthNormalizer;
	double effStress = Compute_Interface_EffectiveStress(this, sigmaCModified, sigma_I_nt_parts, strengthNormalizer);
	scalar_Vals[s1i_effectiveStrs] = effStress;

	scalar_Vals[s1i_contactRel] = 1.0;
	if (contactOffOnMix == sl_contact_off)
		scalar_Vals[s1i_contactRel] = 0.0;
	scalar_Vals[s1i_stickRel] = 1.0;
	if (slipOffOnMix == sl_slip_on)
		scalar_Vals[s1i_stickRel] = 0.0;

	maxEffDelU = 0.0;

	if (damageOffOnMix == sl_interfacial_damage_off)
	{
		initial_damage = 0.0;
		initial_damage_source = 0.0;
		scalar_Vals[s1i_absStick] = 1.0;
	}
	else
	{
		if (damageOffOnMix == sl_interfacial_damage_on)
		{
			if (damageTractionModel == sl_interfacial_damage_traction_TSR)
			{
				double deltaCTmp = deltaCFactor * deltaC;
				if ((tsrModel == tsr_Ortiz) || (tsrModel == tsr_linear))
					maxEffDelU = deltaCTmp * initial_damage;
				else if (tsrModel == tsr_Xu_Needleman)
				{
					maxEffDelU = 7.6383520679940 * deltaCTmp * initial_damage; // 0.01 sigmaC
		//			maxEffDelU = 10.2334134764515 * deltaCTmp * initial_damage; // 0.001 sigmaC
				}
				else if (fabs(initial_damage) > 1e-4)
				{
					THROW("Update function for TSR model so initial damage is mapped to some nonzero maxDelU\n");
				}
			}
			initial_damage = 1.0;
			initial_damage_source = 0.0;
		}
		scalar_Vals[s1i_D] = initial_damage;
		scalar_Vals[s1i_DsrcAllix] = initial_damage_source;
		scalar_Vals[s1i_Dsrc] = initial_damage_source;
		double relCont = scalar_Vals[s1i_contactRel];
		double relStick = scalar_Vals[s1i_stickRel];
		double omd = 1.0 - initial_damage;
		double drelCont = initial_damage * relCont;
		double aI = omd + drelCont * relStick;
		double aII = drelCont * (1.0 - relStick);
		double aIII = 1.0 - aI - aII; // 1 - omd - initial_damage * relCont = initial_damage * (1 - relCont) 
		scalar_Vals[s1i_absStick] = aI;
		scalar_Vals[s1i_absSlip] = aII;
		scalar_Vals[s1i_absSep] = aIII;
	}
}

void SL_Interface_Fracture_PF::Get_sigmaC_deltaC_phiC_scales(double& sigmaCScale, double& deltaCScale, double& energyCScale) const
{
	sigmaCScale = gen_strength;
	deltaCScale = deltaC;
	energyCScale = sigmaCScale * deltaCScale;
	if (TRS_Tensile_Mode())
	{
		if (tsrModel == tsr_Ortiz)
			energyCScale *= 0.5;
		else if (tsrModel == tsr_Xu_Needleman)
			energyCScale *= exp(1.0);
		else
		{
			THROW("Imeplement the option here or leave the factoring out\n");
		}
	}
}


