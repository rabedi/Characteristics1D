#include "SLBulk_Properties.h"
#include "LAfuncsFinalStep.h"
#include "SL_OneInterfaceAllTimes.h"
#include "SLDescriptorData.h"
#include "Domain_AllInterfacesAllTimes.h"

OnePoint_Data_Upstream_w_lOr_r::OnePoint_Data_Upstream_w_lOr_r()
{
#if HAVE_SOURCE
	setValue(source_v, 0.0);
	setValue(source_sigma, 0.0);
#endif
}

OnePoint_Data_Downstream_w_lOr_r::OnePoint_Data_Downstream_w_lOr_r()
{
#if HAVE_SOURCE
	setValue(source_v, 0.0);
	setValue(source_sigma, 0.0);
#endif
}

OnePoint_inBulk_Fields::OnePoint_inBulk_Fields()
{
	damage = 0.0; // 1e40;
	delu0 = 0.0;
}

void OnePoint_inBulk_Fields::OnePoint_inBulk_Fields_Print_Data(ostream & out, bool print_x, bool print_time)
{
	if (print_time)
		out << time << '\t';
	if (print_x)
		out << x << '\t';
#if DiM1
	if (!g_domain->b_ring_opened1D)
		out << u[0] << '\t' << v[0] << '\t' << eps[0] << '\t' << sigma[0];
	else
	{
		double ax = g_SL_desc_data.load_parameters[0] * x, axt = ax * time;
		double tmp = u[0] - axt;
		out << tmp << '\t';
		tmp = v[0] - ax;
		out << tmp << '\t';
		out << eps[0] << '\t' << sigma[0];
	}
#if HAVE_SOURCE
	out << '\t' << source_v[0] << '\t' << source_sigma[0];
#endif
#else
	out << u[0];
	for (int i = 1; i < DiM; ++i)
		out << '\t' << u[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << v[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << eps[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << sigma[i];

#if HAVE_SOURCE
	for (int i = 0; i < DiM; ++i)
		out << '\t' << source_v[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << source_sigma[i];
#endif
#endif
#if RING_PROBLEM
	out << '\t' << v_r;
#endif
	out << '\t' << damage;
	out << '\t' << delu0;
	out << '\n';
}

void OnePoint_inBulk_Fields::OnePoint_inBulk_Fields_Print_Header(ostream& out, bool print_segi_pi, bool print_x, bool print_time)
{
	vector<HeaderLabels> hLabels;
	OnePoint_inBulk_Fields_Data_Get_Header_Labels(hLabels, print_segi_pi, print_x, print_time);
	PrintHeader(hLabels, out);
}




void OnePoint_inBulk_Fields::OnePoint_inBulk_Fields_Data_Get_Header_Labels(vector<HeaderLabels>& hLabels, bool print_segi_pi, bool print_x, bool print_time)
{
	HeaderLabels hLabel;
	string group = "index";
	hLabel.group = group;

	if (print_segi_pi)
	{
		hLabel.subgroup = "segi";
		hLabel.textLabel = "segi";
		hLabel.latexLabel = "i_{\\mathrm{seg}}";
		hLabels.push_back(hLabel);

		hLabel.subgroup = "pti";
		hLabel.textLabel = "pti";
		hLabel.latexLabel = "i_{\\mathrm{pt}}";
		hLabels.push_back(hLabel);
	}
	if (print_time)
	{
		hLabel.subgroup = "time";
		hLabel.textLabel = "time";
		hLabel.latexLabel = "t";
		hLabels.push_back(hLabel);
	}
	if (print_x)
	{
		hLabel.subgroup = "x";
		hLabel.textLabel = "x";
		hLabel.latexLabel = "x";
		hLabels.push_back(hLabel);
	}

	group = "fld";
	hLabel.group = group;

#if DiM1
	hLabel.subgroup = "u";
	hLabel.textLabel = "u";
	hLabel.latexLabel = "u";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "v";
	hLabel.textLabel = "v";
	hLabel.latexLabel = "v";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "eps";
	hLabel.textLabel = "eps";
	hLabel.latexLabel = "\\epsilon";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "sigma";
	hLabel.textLabel = "sigma";
	hLabel.latexLabel = "\\sigma";
	hLabels.push_back(hLabel);
#if HAVE_SOURCE
	hLabel.subgroup = "source_v";
	hLabel.textLabel = "source_v";
	hLabel.latexLabel = "v_{\\mathrm{src}}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "source_sigma";
	hLabel.textLabel = "source_sigma";
	hLabel.latexLabel = "\\sigma_{\\mathrm{src}}";
	hLabels.push_back(hLabel);
#endif
#else
	for (int i = 0; i < DiM; ++i)
	{
			string tmps = "u" + pINDs[i];
			hLabel.textLabel = tmps;
			tmps = "u_" + pINDs[i];
			hLabel.latexLabel = tmps;
			hLabels.push_back(hLabel);
	}
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "v" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "v_" + pINDs[i];
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "eps" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\epsilon_" + pINDs[i];
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "sigma" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\sigma_" + pINDs[i];
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
#endif
#if RING_PROBLEM
	hLabel.subgroup = "v_r";
	hLabel.textLabel = "v_r";
	hLabel.latexLabel = "v_r";
	hLabels.push_back(hLabel);
#endif
	hLabel.group = "damage";
	hLabel.subgroup = "damage";
	hLabel.textLabel = "damage";
	hLabel.latexLabel = "D";
	hLabels.push_back(hLabel);


	hLabel.group = "fld";
	hLabel.subgroup = "delu";
	hLabel.textLabel = "delu0";
	hLabel.latexLabel = "\\Delta{u}_0";
	hLabels.push_back(hLabel);
}

double	getCdIsotropic(double rho, double nu, double E)
{
#if DiM1
	return sqrt(E / rho);
#endif
	double A;
#if DiM2
#if USE_PLANE_STRAIN
	A = E * (1.0 - nu) / (1.0 + nu) / (1 - 2.0 * nu);
#else
	A = E / (1.0 - nu * nu);
#endif
#endif
#if DiM3
	A = E * (1.0 - nu) / (1.0 + nu) / (1 - 2.0 * nu);
#endif
	return sqrt(A / rho);
}

double	getCsIsotropic(double rho, double nu, double E)
{
#if DiM1
	return sqrt(E / rho);
#endif
	double C = E / (1.0 + nu);
	return sqrt(C / 2.0 / rho);
}

SL_Bulk_Properties::SL_Bulk_Properties()
{
	E_iso = -1.0;
	nu_iso = -1.0;
	rho = 1.0;
	is_iso = true;
	flag = 1;

#if HAVE_SOURCE_ORDER0_q
	D_vv = 0.0, D_vsigma = 0.0, D_sigmav = 0.0, D_sigmasigma = 0.0;
	has_DTerms = false;
#endif
	scaling_isDone = false;
}

void SL_Bulk_Properties::Read_SL_Bulk_Properties(istream& in, int bulkFlag)
{
	flag = bulkFlag;
	string str_flag;
	toString(flag, str_flag);
	bool versionChange = ((sfcm.success) && (sfcm_gen.bulkFlag == flag));
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
		if (buf == "rho")
		{
			READ_NDOUBLE(in, buf, rho);
			key = "rho";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
			{
				rho = value;
			}
			else
			{
				key = "rho" + str_flag;
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					rho = value;
				}
				else if (g_SL_desc_data.tdLoadComputer != NULL)
				{
					if ((flag == g_bulkFlags_IncImp[SDL]) || (flag == g_bulkFlags_IncImp[SDR]))
					{
						key = "rhoAmbient";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							rho = value;
						}
					}
				}
			}
		}
		else if (buf == "flag")
		{
			READ_NINTEGER(in, buf, flag);
		}
		else if (buf == "nu_iso")
		{
			READ_NDOUBLE(in, buf, nu_iso);
		}
		else if ((buf == "E_iso") || (buf == "C_iso"))
		{
			READ_NDOUBLE(in, buf, E_iso);
			key = "E";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
			{
				E_iso = value;
			}
			else
			{
				key = "E" + str_flag;
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					E_iso = value;
				}
				else if (g_SL_desc_data.tdLoadComputer != NULL)
				{
					if ((flag == g_bulkFlags_IncImp[SDL]) || (flag == g_bulkFlags_IncImp[SDR]))
					{
						key = "EAmbient";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							E_iso = value;
						}
					}
				}
			}
		}
		else if (buf == "E_iso_4_plainStrain_sym_y")
		{
			double E;
			READ_NDOUBLE(in, buf, E);
			E_iso = E * (1.0 - nu_iso) / (1.0 + nu_iso) / (1 - 2.0 * nu_iso);
		}
#if HAVE_SOURCE_ORDER0_q
		else if (buf == "D_vv")
		{
			READ_NDOUBLE(in, buf, D_vv);
			key = "D_vv";
			if (versionChange && (Find_Version_Value(key, value, mpPtr) == true))
			{
				D_vv = value;
			}
			else
			{
				key = "D_vv" + str_flag;
				if (Find_Version_Value(key, value, mpPtr) == true)
				{
					D_vv = value;
				}
				else if (g_SL_desc_data.tdLoadComputer != NULL)
				{
					if ((flag == g_bulkFlags_IncImp[SDL]) || (flag == g_bulkFlags_IncImp[SDR]))
					{
						key = "D_vvAmbient";
						if (Find_Version_Value(key, value, mpPtr) == true)
						{
							D_vv = value;
						}
					}
				}
			}
		}
		else if (buf == "D_vsigma")
		{
			READ_NDOUBLE(in, buf, D_vsigma);
		}
		else if (buf == "D_sigmav")
		{
			READ_NDOUBLE(in, buf, D_sigmav);
		}
		else if (buf == "D_sigmasigma")
		{
			READ_NDOUBLE(in, buf, D_sigmasigma);
		}
#endif
#if COMPUTE_ANISO_BULK
		else if (buf == "C")
			ReadM(C, in);
#endif
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Initialize_FromInputParas();
}

void SL_Bulk_Properties::Print(ostream &out, bool printComputedVals) const
{
	out << "{\n";
	out << "flag\t" << flag << '\n';
	out << "rho\t" << rho << '\n';
	out << "E_iso\t" << E_iso << '\n';
	out << "nu_iso\t" << nu_iso << '\n';
#if HAVE_SOURCE_ORDER0_q
	out << "has_DTerms\t" << has_DTerms << '\n';
	out << "D_vv\t" << D_vv << '\n';
	out << "D_vsigma\t" << D_vsigma << '\n';
	out << "D_sigmav\t" << D_sigmav << '\n';
	out << "D_sigmasigma\t" << D_sigmasigma << '\n';
#endif
#if COMPUTE_ANISO_BULK
	out << "C\n";
	PrintM(C, out);
	out << "\n";
#endif
	out << "}\n";
	if (!printComputedVals)
		return;
	out << "[\n";
	out << "is_iso\t" << is_iso << '\n';
	out << "cd_iso\t" << cd_iso << '\n';
	out << "ws\n";		PrintV(ws, out);	out << "\n";
	out << "c_rhos\n";		PrintV(c_rhos, out);	out << "\n";
	out << "inv_c_rhos\n";		PrintV(inv_c_rhos, out);	out << "\n";
#if HAVE_SOURCE_ORDER0_q
	out << "allDiagonal_Ds_are_zero\t" << allDiagonal_Ds_are_zero << '\n';
	out << "D_wlwl\n";		PrintV(D_wlwl, out);	out << "\n";
	out << "D_wrwr\n";		PrintV(D_wrwr, out);	out << "\n";
	out << "D_wlwr\n";		PrintV(D_wlwr, out);	out << "\n";
	out << "D_wrwl\n";		PrintV(D_wrwl, out);	out << "\n";
#endif
#if COMPUTE_ANISO_BULK
	out << "C_over_rho\n";		PrintM(C_over_rho, out);	out << "\n";
	out << "E\n";		PrintM(E, out);	out << "\n";
	out << "C_over_rho\n";		PrintM(C_over_rho, out);	out << "\n";
	out << "F\n";		PrintM(F, out);	out << "\n";
	out << "H\n";		PrintM(H, out);	out << "\n";
	out << "Z\n";		PrintM(Z, out);	out << "\n";
	out << "Y\n";		PrintM(Y, out);	out << "\n";
#endif
	out << "]\n";
}

void SL_Bulk_Properties::Initialize_FromInputParas()
{
	is_iso = (E_iso > 0.0);
#if USE_ISO_ASSUMPTION
	if (!is_iso)
		THROW("For iso assumption, E_iso must be > 0 to enter iso block below\n");
#endif
	double rhoInv = 1.0 / rho;
	if (is_iso) // material is isotropic 
	{
		cd_iso = getCdIsotropic(rho, nu_iso, E_iso);
		cs_iso = cd_iso;
		c_rhos[0] = cd_iso * rho;
		ws[0] = cd_iso;
#if DiM2a3_F		
		cs_iso = getCsIsotropic(rho, nu_iso, E_iso);
		ws[1] = cs_iso;
		c_rhos[1] = cs_iso * rho;
#if DiM3
		ws[2] = cs_iso;
		c_rhos[2] = c_rhos[1];
#endif
#endif
#if COMPUTE_ANISO_BULK
		double C00 = c_rhos[0] * cd_iso;
		C[0][0] = C00;
		double z0 = c_rhos[0];
		E[0][0] = 1.0;
		F[0][0] = z0;
		Z[0][0] = z0;
		double zInv0 = 1.0 / z0;
		H[0][0] = zInv0;
		Y[0][0] = zInv0;
#if DiM2a3_F		
		double G = c_rhos[1] * cs_iso;
		C[0][1] = 0.0;
		C[1][0] = 0.0;
		C[1][1] = G;
#if DiM3
		C[1][2] = 0.0;
		C[2][1] = 0.0;
		C[0][2] = 0.0;
		C[2][0] = 0.0;
		C[2][2] = G;
		ws[2] = cs_iso;
		c_rhos[2] = c_rhos[1];
#endif
		double z, zInv;
		for (int i = 0; i < DiM; ++i)
		{
			for (int j = 0; j < DiM; ++j)
			{
				C_over_rho[i][j] = C[i][j] * rhoInv;
				E[i][j] = 0.0;
				F[i][j] = 0.0;
				H[i][j] = 0.0;
				Z[i][j] = 0.0;
				Y[i][j] = 0.0;
			}
			E[i][i] = 1.0;
			z = c_rhos[i];
			F[i][i] = z;
			Z[i][i] = z;
			zInv = 1.0 / z;
			H[i][i] = zInv;
			Y[i][i] = zInv;
		}
#endif
		Inverse(C, inv_C);
#endif
		for (int i = 0; i < DiM; ++i)
			effective_CsIso[i] = ws[i] * ws[i] * rho;
	}
#if COMPUTE_ANISO_BULK
	else
	{
		// anisotropic case
		// C must have been provided
		for (int i = 0; i < DiM; ++i)
			for (int j = 0; j < DiM; ++j)
				C_over_rho[i][j] = C[i][j] * rhoInv;

		VEC eigenValues;
		if (!EigenValuesEigenVectors_PositiveSymMatrix(C_over_rho, E, eigenValues) == false)
		{
			THROW("Cannot solve eigenvalues of C_over_rho\n");
		}
		for (int i = 0; i < DiM; ++i)
			ws[i] = sqrt(eigenValues[i]);

		double crho;
		for (int i = 0; i < DiM; ++i)
		{
			crho = rho * ws[i];
			c_rhos[i] = crho;
			for (int j = 0; j < DiM; ++j)
			{
				F[i][j] = crho * E[i][j];
			}
		}
		Inverse(F, H);

		for (int i = 0; i < DiM; ++i)
		{
			for (int j = 0; j < DiM; ++j)
			{
				Z[i][j] = 0.0;
				for (int k = 0; k < DiM; ++k)
					Z[i][j] += E[k][i] * F[k][j];
			}
		}
		Inverse(Z, Y);
	}
#endif
	for (int i = 0; i < DiM; ++i)
	{
		inv_c_rhos[i] = 1.0 / c_rhos[i];
		inv_ws[i] = 1.0 / ws[i];
	}
#if HAVE_SOURCE_ORDER0_q
	// D_wlwl = 0.5 * [(D_sigmasigma + D_vv) + (Dprime_sigmav + Dprime_vsigma)]
	// D_wrwr = 0.5 * [(D_sigmasigma + D_vv) - (Dprime_sigmav + Dprime_vsigma)]
	// D_wlwr = 0.5 * [(D_sigmasigma - D_vv) - (Dprime_sigmav - Dprime_vsigma)]
	// D_wrwl = 0.5 * [(D_sigmasigma - D_vv) + (Dprime_sigmav - Dprime_vsigma)]
	//	Dprime_sigmav = D_sigmav / crho
	//	Dprime_vsigma = D_vsigma * crho
	VEC half_Dprime_sigmav, half_Dprime_vsigma, half_Dprime_sigmav_p_vsigma, half_Dprime_sigmav_m_vsigma;
	double half_ss_p_vv = 0.5 * (D_sigmasigma + D_vv), half_ss_m_vv = 0.5 * (D_sigmasigma - D_vv);
	static double tol = 1e-18;
	if (fabs(D_vv) > tol)
		has_DTerms = true;
	else if (fabs(D_vsigma) > tol)
		has_DTerms = true;
	else if (fabs(D_sigmav) > tol)
		has_DTerms = true;
	else if (fabs(D_sigmasigma) > tol)
		has_DTerms = true;

	allDiagonal_Ds_are_zero = true;
	double tmp; tol = MAX(fabs(half_ss_p_vv), 1e-8);
	for (int i = 0; i < DiM; ++i)
	{
		half_Dprime_sigmav[i] = D_sigmav * inv_c_rhos[i];
		half_Dprime_vsigma[i] = D_vsigma *     c_rhos[i];
		half_Dprime_sigmav_p_vsigma[i] = 0.5 * (half_Dprime_sigmav[i] + half_Dprime_vsigma[i]);
		half_Dprime_sigmav_m_vsigma[i] = 0.5 * (half_Dprime_sigmav[i] - half_Dprime_vsigma[i]);

		D_wlwl[i] = half_ss_p_vv + half_Dprime_sigmav_p_vsigma[i];
		D_wrwr[i] = half_ss_p_vv - half_Dprime_sigmav_p_vsigma[i];
		tmp = half_ss_m_vv - half_Dprime_sigmav_m_vsigma[i];
		if (fabs(tmp) > tol)
			allDiagonal_Ds_are_zero = false;
		D_wlwr[i] = tmp;

		tmp = half_ss_m_vv + half_Dprime_sigmav_m_vsigma[i];
		if (fabs(tmp) > tol)
			allDiagonal_Ds_are_zero = false;
		D_wrwl[i] = tmp;
	}
#endif	
}

void SL_Bulk_Properties::Initialize_FromOther_withFactors(const SL_Bulk_Properties& other, double factorC, double factorRho, double factor_damping)
{
	(*this) = other;
	if (scaling_isDone)
		return;
	double factor_c2 = factorC / factorRho;
	double factor_c = sqrt(factor_c2);
	double factor_Z = factorRho * factor_c;
	double factor_Y = 1.0 / factor_Z;
	rho *= factorRho;
#if HAVE_SOURCE_ORDER0_q
	// the following is not technically correct. Damping factor is introduced loosely and in practice is intended to scale all damping factors in characteristics system ...
	if (has_DTerms)
	{
		D_vv *= factor_damping;
		D_vsigma *= factor_damping;
		D_sigmav *= factor_damping;
		D_sigmasigma *= factor_damping;
	}
#endif
	E_iso *= factorC;
	cd_iso *= factor_c;
	cs_iso *= factor_c;

#if COMPUTE_ANISO_BULK
	FactorMat(C, factorC);
	FactorMat(C_over_rho, factor_c2);
	FactorMat(Z, factor_Z);
	FactorMat(Y, factor_Y);
	FactorMat(F, factor_Z);
	FactorMat(H, factor_Y);
#endif
	FactorVec(ws, factor_c);
	FactorVec(c_rhos, factor_Z);
	FactorVec(effective_CsIso, factorC);
	FactorVec(inv_c_rhos, factor_Y);
#if HAVE_SOURCE_ORDER0_q
	// damping factors are individually scaled such that the whole damping matrix in w system is scaled ...
	FactorVec(D_wlwl, factor_damping);
	FactorVec(D_wrwr, factor_damping);
	FactorVec(D_wlwr, factor_damping);
	FactorVec(D_wrwl, factor_damping);
#endif
	scaling_isDone = true;
}

void SL_Bulk_Properties::Form_Iso_E1Rho1Default()
{
	E_iso = 1.0;
	nu_iso = 0.3;
	rho = 1.0;
	is_iso = true;
	Initialize_FromInputParas();
}

double SL_Bulk_Properties::q_to_w(const VEC & v, const VEC & sigma, int w_index, bool is_wr)
{
	if (is_iso)
	{
		if (is_wr)
			return sigma[w_index] - c_rhos[w_index] * v[w_index];
		else
			return sigma[w_index] + c_rhos[w_index] * v[w_index];
	}
#if	!USE_ISO_ASSUMPTION
	VEC tmp;
	for (int i = 0; i < DiM; ++i)
		tmp[i] = sigma[i];
	if (is_wr)
	{
		for (int i = 0; i < DiM; ++i)
			for (int j = 0; j < DiM; ++j)
				tmp[i] -= Z(i, j) * v[i];
	}
	else
	{
		for (int i = 0; i < DiM; ++i)
			for (int j = 0; j < DiM; ++j)
				tmp[i] += Z(i, j) * v[i];
	}
/*	if (is_wr)
		for (int i = 0; i < DiM; ++i)
			tmp[i] = sigma[i] - c_rhos[w_index] * v[i];
	else
		for (int i = 0; i < DiM; ++i)
			tmp[i] = sigma[i] + c_rhos[w_index] * v[i];
*/
	double ret = 0.0;
	for (int i = 0; i < DiM; ++i)
		ret += E[w_index][i] * tmp[i];
	return  ret;
#endif
	THROW("Invalid iso/aniso assumption\n");
}

void SL_Bulk_Properties::q_to_w_vec(const VEC& v, const VEC& sigma, VEC& w_rOrlGoing, bool is_wr)
{
	if (is_iso)
	{
		if (is_wr)
			for (int i = 0; i < DiM; ++i)
				w_rOrlGoing[i] = sigma[i] - c_rhos[i] * v[i];
		else
			for (int i = 0; i < DiM; ++i)
				w_rOrlGoing[i] = sigma[i] + c_rhos[i] * v[i];
	}
#if	!USE_ISO_ASSUMPTION
	VEC Fv, Esigma;
	ProductMatVec(F, v, Fv);
	ProductMatVec(E, sigma, Esigma);
	if (is_wr)
		SubtractVec(Esigma, Fv, w_rOrlGoing);
	else
		AddVec(Esigma, Fv, w_rOrlGoing);
#endif
}

void SL_Bulk_Properties::q_to_w_2vecs(const VEC& v, const VEC& sigma, VEC& w_rGoing, VEC& w_lGoing)
{
	if (is_iso)
	{
		for (int i = 0; i < DiM; ++i)
			w_rGoing[i] = sigma[i] - c_rhos[i] * v[i];
		for (int i = 0; i < DiM; ++i)
			w_lGoing[i] = sigma[i] + c_rhos[i] * v[i];
	}
#if	!USE_ISO_ASSUMPTION
	VEC Fv, Esigma;
	ProductMatVec(F, v, Fv);
	ProductMatVec(E, sigma, Esigma);
	SubtractVec(Esigma, Fv, w_rGoing);
	AddVec(Esigma, Fv, w_lGoing);
#endif
}

void SL_Bulk_Properties::w_to_q_2vecs(const VEC& w_rGoing, const VEC& w_lGoing, VEC& v, VEC& sigma)
{
	if (is_iso)
	{
		for (int i = 0; i < DiM; ++i)
			sigma[i] = 0.5 * (w_lGoing[i] + w_rGoing[i]);
		for (int i = 0; i < DiM; ++i)
			v[i] = 0.5 * inv_c_rhos[i] * (w_lGoing[i] - w_rGoing[i]);
	}
#if !USE_ISO_ASSUMPTION
	// v = 0.5 H (wl - wr);		sigma = 0.5 E^t (wl + wr)
	VEC half_wl_minus_wr, half_wl_plus_wr;
	SubtractVec(w_lGoing, w_rGoing, half_wl_minus_wr);
	AddVec(w_lGoing, w_rGoing, half_wl_plus_wr);
	FactorVec(half_wl_plus_wr, 0.5);
	FactorVec(half_wl_minus_wr, 0.5);
	ProductMatVec(H, half_wl_minus_wr, v);
	ProductTransposeMatVec(E, half_wl_plus_wr, sigma);
#endif
}

void SL_Bulk_Properties::Compute_stress_from_v_and_characteristics(const VEC& w, const VEC& v, bool is_wr, VEC& sigma)
{
	if (is_iso)
	{
		if (is_wr) // left side
		{
			for (int i = 0; i < DiM; ++i)
				sigma[i] = w[i] + c_rhos[i] * v[i];
		}
		else
		{
			for (int i = 0; i < DiM; ++i)
				sigma[i] = w[i] - c_rhos[i] * v[i];
		}
	}
#if !USE_ISO_ASSUMPTION
	// there are two parts
	VEC Etw, Zv;
	ProductTransposeMatVec(E, w, Etw);
	ProductMatVec(Z, v, Zv);
	if (is_wr) // left side: sigma = EtW + Zv
	{
		for (int i = 0; i < DiM; ++i)
			sigma[i] = Etw[i] + Zv[i];
	}
	else
		if (is_wr) // left side: sigma = EtW + Zv
		{
			for (int i = 0; i < DiM; ++i)
				sigma[i] = Etw[i] - Zv[i];
		}
#endif
}

void SL_Bulk_Properties::Compute_Downstream_characteristic_From_Upstream_qs_and_Forces(AllPoints_Data_Upstream_w_lOr_r& upstream_pts, OnePoint_Data_Downstream_w_lOr_r& downstream_pt, bool is_wr, double x, bool useDSCharsAsUSChars, double* ring_opened1D_alPtr)
{
	bool hasVelJump = (ring_opened1D_alPtr != NULL);
	//! Compute upstream characteristics from q (v, sigma) or all upstream points
	if (!hasVelJump)
	{
		for (int i = 0; i < DiM; ++i)
			upstream_pts.wUpstream[i] = q_to_w(upstream_pts.up_Pts[i].v, upstream_pts.up_Pts[i].sigma, i, is_wr);
	}
	else
	{
		VEC vel;
		CopyVec(upstream_pts.up_Pts[0].v, vel);
		vel[0] += *ring_opened1D_alPtr;
		upstream_pts.wUpstream[0] = q_to_w(vel, upstream_pts.up_Pts[0].sigma, 0, is_wr);
	}
	//! If the problem does not have any source term, this is directly translated to downstream w:
#if !HAVE_SOURCE
	CopyVec(upstream_pts.wUpstream, downstream_pt.wDownstream);
	return;
#else
	if (useDSCharsAsUSChars)
	{
		CopyVec(upstream_pts.wUpstream, downstream_pt.wDownstream);
		return;
	}

	// problem has source term -> upstream w changes as it gets to downstream w
	//! Calculating upstream force from individual q's
	for (int i = 0; i < DiM; ++i)
		upstream_pts.source_w_Upstream[i] = q_to_w(upstream_pts.up_Pts[i].source_v, upstream_pts.up_Pts[i].source_sigma, i, is_wr);
	//! Calculating downstream force again from source_v, source_q (they may be zero)
	VEC sources_w_Downstream, souces_w_ave_UpDownStreams;
	q_to_w_vec(downstream_pt.source_v, downstream_pt.source_sigma, sources_w_Downstream, is_wr);
	for (int i = 0; i < DiM; ++i)
		souces_w_ave_UpDownStreams[i] = 0.5 * (upstream_pts.source_w_Upstream[i] + sources_w_Downstream[i]);

	// now compute downstream w
	// simple case: no 0th order q's in the problem
	#if !HAVE_SOURCE_ORDER0_q
	for (int i = 0; i < DiM; ++i)
		downstream_pt.wDownstream[i] = upstream_pts.wUpstream[i] + souces_w_ave_UpDownStreams[i] * upstream_pts.up_Pts[i].delT;
	#else
	// more difficult case, as there are diagonal d terms from 0th order q
	// see if an exact solution can be devised, depending on whether off-diagonal components of D_omega matrix are zero or not
	if (allDiagonal_Ds_are_zero && 	(!g_domain->b_ring_opened1D)) // better case, with analytical solution
	{
		VEC ds; // diagoinal q's for characteristics
		if (is_wr)
			CopyVec(D_wrwr, ds);
		else
			CopyVec(D_wlwl, ds);
		// components of exp(-d del_t)
		double exp_m_delT;
		double w_upstream, source_upstream, source_downstream;
		for (int i = 0; i < DiM; ++i)
		{
			double delT = upstream_pts.up_Pts[i].delT;
			exp_m_delT = exp(-ds[i] * delT);
			w_upstream = upstream_pts.wUpstream[i];
			source_upstream = upstream_pts.source_w_Upstream[i];
			source_downstream = sources_w_Downstream[i];
			// w_up * exp(-d delt) + 0.5 delt * (source_up exp(-d delt) + source_down)
			downstream_pt.wDownstream[i] = w_upstream * exp_m_delT + 0.5 * delT * (source_upstream * exp_m_delT + source_downstream);
		}
	}
	else // this is the case where estimate of q for downstream should be used
	{
		VEC added_source_w;
		for (int i = 0; i < DiM; ++i)
		{
			// computing mid-point q (v, sigma) along each characteristics
			VEC ave_v, ave_sigma;
			for (int j = 0; j < DiM; ++j)
			{
				ave_v[j] = 0.5 * (downstream_pt.last_estimate_v[j] + upstream_pts.up_Pts[i].v[j]);
				ave_sigma[j] = 0.5 * (downstream_pt.last_estimate_sigma[j] + upstream_pts.up_Pts[i].sigma[j]);
			}
			if (g_domain->b_ring_opened1D)
			{
				if (g_domain->b_ring_opened1D_damping_on_full_vTheta) // computing true v_theta
					ave_v[0] -= g_SL_desc_data.load_parameters[0] * x;
				else // damping computing from v vTheta + ax which does not really make sense
				{
					if (hasVelJump) // left boundary, velocity is incremented by aL
						ave_v[0] += *ring_opened1D_alPtr;
				}
			}
			// computing the corresponding source (-Dq)
			VEC ave_src_v, ave_src_sigma;
			for (int j = 0; j < DiM; ++j)
			{
				ave_src_v[j] = -(D_vv	  * ave_v[j] + D_vsigma * ave_sigma[j]);
				ave_src_sigma[j] = -(D_sigmav * ave_v[j] + D_sigmasigma * ave_sigma[j]);
			}
			// turning the source in q to source in w
			VEC source_vec_w;
			q_to_w_vec(ave_src_v, ave_src_sigma, source_vec_w, is_wr);
			// for this characeristics direction i, only component i is relevant and will be used below
			added_source_w[i] = source_vec_w[i];
		}
// This is the obsolete version that uses D coefficients in w parameter. It's replaced with D's in q parameter above
#if 0
		VEC w_downstream_same_side, w_downstream_other_side;
		q_to_w_vec(downstream_pt.last_estimate_v, downstream_pt.last_estimate_sigma, w_downstream_same_side, is_wr);
		q_to_w_vec(downstream_pt.last_estimate_v, downstream_pt.last_estimate_sigma, w_downstream_other_side, !is_wr);
		VEC w_same_side;
		for (int i = 0; i < DiM; ++i)
			w_same_side[i] = 0.5 * (w_downstream_same_side[i] + upstream_pts.wUpstream[i]);
		VEC *d_same_side, *d_other_side;
		if (is_wr)
		{
			d_same_side = &D_wrwr;
			d_other_side = &D_wrwl;
		}
		else
		{
			d_same_side = &D_wlwl;
			d_other_side = &D_wlwr;
		}
#endif
		double prev_src, added_src, corrected_src;
		for (int i = 0; i < DiM; ++i)
		{
			prev_src = souces_w_ave_UpDownStreams[i];
			added_src = added_source_w[i];
//			added_src = -((*d_same_side)[i] * w_same_side[i] + (*d_other_side)[i] * w_downstream_other_side[i]);
			corrected_src = prev_src + added_src;
//			for (int i = 0; i < DiM; ++i)
			downstream_pt.wDownstream[i] = upstream_pts.wUpstream[i] + corrected_src * upstream_pts.up_Pts[i].delT;
		}
	}
	#endif
#endif
}

void SL_Bulk_Properties::Compute_Downstream_characteristic_From_in_situ_Stress(const VEC& in_situ_stress, const VEC& wDownstream_WO_in_situ, VEC& wDownstream_W_in_situ)
{
	if (is_iso)
	{
		for (int i = 0; i < DiM; ++i)
			wDownstream_W_in_situ[i] = wDownstream_WO_in_situ[i] + in_situ_stress[i];
	}
#if !USE_ISO_ASSUMPTION
	double tmp;
	for (int i = 0; i < DiM; ++i)
	{
		tmp = wDownstream_WO_in_situ[i];
		// w = E sigma +/- F v, so added contribution from in-situ stress is E in-situ-stress
		for (int j = 0; j < DiM; ++j)
			tmp += E[i][j] * in_situ_stress[j];
		wDownstream_W_in_situ[i] += tmp;
	}
#endif
}

void SL_Bulk_Properties::Compute_Stress_from_Strain(const VEC& eps, VEC& sigma)
{
	if (is_iso)
	{
		for (int i = 0; i < DiM; ++i)
			sigma[i] = eps[i] * effective_CsIso[i];
	}
#if !USE_ISO_ASSUMPTION
	else
		ProductMatVec(C, eps, sigma);
#endif
}

void SL_Bulk_Properties::Compute_Strain_from_Stress(const VEC& sigma, VEC& eps)
{
	if (is_iso)
	{
		for (int i = 0; i < DiM; ++i)
			eps[i] = sigma[i] / effective_CsIso[i];
	}
#if !USE_ISO_ASSUMPTION
	else
		ProductMatVec(inv_C, sigma, eps);
}
#endif
}

void SL_Bulk_Properties::Compute_Bulk_vseps_IC(OnePoint_inBulk_Fields& point_fieldsNT)
{
	g_SL_desc_data.Get_IC(this, flag, point_fieldsNT.x, point_fieldsNT.u, point_fieldsNT.v, point_fieldsNT.eps
#if RING_PROBLEM
		, point_fieldsNT.v_r
#endif
	);
	Compute_Stress_from_Strain(point_fieldsNT.eps, point_fieldsNT.sigma);
}

void SL_Bulk_Properties::Compute_Bulk_vseps_NonIC(double x, SL_OneInterfaceAllTimes *interfaceLeftOfBulk, SL_OneInterfaceAllTimes* interfaceRightOfBulk, OnePoint_inBulk_Fields& point_fieldsNT, OnePoint_inBulk_Fields* point_fieldsPT_or_NTPIPtr)
{
	VEC wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ;
	setValue(wlSide_rGoing_WO_in_situ, 0.0);
	setValue(wrSide_lGoing_WO_in_situ, 0.0);
	VEC *wDownstream;
	double *ring_opened1D_alPtr = NULL;

#if RING_PROBLEM
	double v_r;
#endif
	for (int eside = 0; eside < NUM_SIDES; ++eside)
	{
		bool is_wr;
		/////////////////////////////////////////////////////////////////
		// A. Upstream values
		int epside;
		SL_OneInterfaceAllTimes *interfaceSide;
		double delx;
		if (eside == SDL)
		{
			wDownstream = &wlSide_rGoing_WO_in_situ;
			epside = SDR;
			interfaceSide = interfaceLeftOfBulk;
			delx = point_fieldsNT.x - interfaceSide->interface_x;
			if (g_domain->isPeriodic && (delx < 0.0))
				delx += g_domain->L;
			is_wr = true;
		}
		else
		{
			wDownstream = &wrSide_lGoing_WO_in_situ;
			epside = SDL;
			interfaceSide = interfaceRightOfBulk;
			delx = interfaceSide->interface_x - point_fieldsNT.x;
			if (g_domain->isPeriodic && (delx < 0.0))
				delx += g_domain->L;
			is_wr = false;
		}
		AllPoints_Data_Upstream_w_lOr_r upstream_pts;
		for (int i = 0; i < DiM; ++i)
		{
			double delt = delx * inv_ws[i];
			double tAbsolute = point_fieldsNT.time - delt;
			if (tAbsolute > 0) // hit the neighbor time sequence
			{
				upstream_pts.up_Pts[i].delT = delt;
				SL_interfacePPtData* ptSlnPtr;
				bool ptSlnPtr_Deletable;
				bool found = interfaceSide->timeSeqData.Interpolate_Pt_Solution_At_time(tAbsolute, ptSlnPtr, ptSlnPtr_Deletable);
				if (found == false)
				{
					THROW("Cannot find the point\n");
				}

				CopyVec(ptSlnPtr->sl_side_ptData[epside].sigma_downstream_final, upstream_pts.up_Pts[i].sigma);
				CopyVec(ptSlnPtr->sl_side_ptData[epside].v_downstream_final, upstream_pts.up_Pts[i].v);
#if HAVE_SOURCE
				g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(interfaceSide->interface_x, tAbsolute, flag, upstream_pts.up_Pts[i].source_v, upstream_pts.up_Pts[i].source_sigma);
#if RING_PROBLEM
				v_r = ptSlnPtr->v_r_final;
#endif
#endif
				if (ptSlnPtr_Deletable)
					delete ptSlnPtr;
			}
			else // characteristics hits IC
			{
				double delT = point_fieldsNT.time;
				double del_x = ws[i] * delT;
				double xIC;
				if (eside == SDL)
					xIC = point_fieldsNT.x - del_x;
				else
					xIC = point_fieldsNT.x + del_x;

				upstream_pts.up_Pts[i].delT = delT;

				// getting ICs
				VEC u_0, eps_0;
				g_SL_desc_data.Get_IC(this, flag, xIC, u_0, upstream_pts.up_Pts[i].v, eps_0
#if RING_PROBLEM
					, v_r
#endif
				);
				Compute_Stress_from_Strain(eps_0, upstream_pts.up_Pts[i].sigma);

#if HAVE_SOURCE
				g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(xIC, 0.0, flag, upstream_pts.up_Pts[i].source_v, upstream_pts.up_Pts[i].source_sigma);
#endif
			}
#if HAVE_SOURCE
#if RING_PROBLEM
			double sigma_theta_source_final = E_iso * v_r / g_domain->ring_R; // E v_r / R (same on both sides)
			upstream_pts.up_Pts[i].source_sigma[0] += sigma_theta_source_final;
#endif
#endif
		}
		/////////////////////////////////////////////////////////////////
		// B. Downstream values, parts needed
		OnePoint_Data_Downstream_w_lOr_r downstream_pt;
		bool useDSCharsAsUSChars = false;
#if HAVE_SOURCE
		if (point_fieldsPT_or_NTPIPtr != NULL)
		{
			CopyVec(point_fieldsPT_or_NTPIPtr->v, downstream_pt.last_estimate_v);
			CopyVec(point_fieldsPT_or_NTPIPtr->sigma, downstream_pt.last_estimate_sigma);
#if RING_PROBLEM
			v_r = point_fieldsPT_or_NTPIPtr->v_r;
#endif
		}
		else
		{
			// we don't have downstream sigma and v so simply use downstream charateristics
			useDSCharsAsUSChars = true;
			CopyVec(upstream_pts.up_Pts[0].v, downstream_pt.last_estimate_v);
			CopyVec(upstream_pts.up_Pts[0].sigma, downstream_pt.last_estimate_sigma);
			//#if RING_PROBLEM
						// nothing is needed for v_r is has the latest value already
			//#endif
		}
		g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(point_fieldsNT.x, point_fieldsNT.time, flag, downstream_pt.source_v, downstream_pt.source_sigma);
#if RING_PROBLEM
		double sigma_theta_source_final = E_iso * v_r / g_domain->ring_R; // E v_r / R (same on both sides)
		downstream_pt.source_sigma[0] += sigma_theta_source_final;
#endif
#else
		setValue(downstream_pt.last_estimate_v, 0.0);
		setValue(downstream_pt.last_estimate_sigma, 0.0);
#endif
		Compute_Downstream_characteristic_From_Upstream_qs_and_Forces(upstream_pts, downstream_pt, is_wr, x, false, ring_opened1D_alPtr);
		CopyVec(downstream_pt.wDownstream, *wDownstream);
	}
	// C: now having upstream and downstream values calculate the star value for the bulk
	w_to_q_2vecs(wlSide_rGoing_WO_in_situ, wrSide_lGoing_WO_in_situ, point_fieldsNT.v, point_fieldsNT.sigma);
	// D: Compute Strain
	Compute_Strain_from_Stress(point_fieldsNT.sigma, point_fieldsNT.eps);
}

void SL_Bulk_Properties::Compute_Bulk_Values(double x, SL_OneInterfaceAllTimes *interfaceLeftOfBulk, SL_OneInterfaceAllTimes* interfaceRightOfBulk, OnePoint_inBulk_Fields& point_fieldsNT, bool isIC, OnePoint_inBulk_Fields* point_fieldsPT_or_NTPIPtr, unsigned int maxIter, double relTol4Conv)
{
	if (!isIC)
	{
		bool treatWOSourceUpdate = true;
#if HAVE_SOURCE
		treatWOSourceUpdate = (point_fieldsPT_or_NTPIPtr == NULL);
#endif
		if (treatWOSourceUpdate)
			Compute_Bulk_vseps_NonIC(x, interfaceLeftOfBulk, interfaceRightOfBulk, point_fieldsNT, NULL);
		else
		{
			double Zscalar = c_rhos[0];
			unsigned int iter = 0;

			OnePoint_inBulk_Fields prevItVal = *point_fieldsPT_or_NTPIPtr;
			double normValP = MAX(Norm2(prevItVal.sigma), Zscalar * Norm2(prevItVal.v));
			while (iter++ < maxIter)
			{
				Compute_Bulk_vseps_NonIC(x, interfaceLeftOfBulk, interfaceRightOfBulk, point_fieldsNT, &prevItVal);
				VEC delv, delsigma;
				SubtractVec(point_fieldsNT.v, prevItVal.v, delv);
				SubtractVec(point_fieldsNT.sigma, prevItVal.sigma, delsigma);
				double delNorm = MAX(Norm2(delsigma), Zscalar * Norm2(delv));
				double normValN = MAX(Norm2(point_fieldsNT.sigma), Zscalar * Norm2(point_fieldsNT.v));
				double normVal = MAX(normValN, normValP);
				if (delNorm <= (normVal + 1e-50)* relTol4Conv)
					break;
				prevItVal = point_fieldsNT;
				normValP = normValN;
			}
		}
		if (point_fieldsPT_or_NTPIPtr != NULL)
		{
			double delt = point_fieldsNT.time - point_fieldsPT_or_NTPIPtr->time;
			double hdelt = 0.5 * delt;
			// updating u
			for (int i = 0; i < DiM; ++i)
				point_fieldsNT.u[i] = point_fieldsPT_or_NTPIPtr->u[i] + hdelt * (point_fieldsNT.v[i] + point_fieldsPT_or_NTPIPtr->v[i]);
#if RING_PROBLEM
			double ring_p = 0.5 * (g_SL_desc_data.GetRing_p(point_fieldsNT.x, point_fieldsNT.time) +
				g_SL_desc_data.GetRing_p(point_fieldsPT_or_NTPIPtr->x, point_fieldsPT_or_NTPIPtr->time));
			double sigma_theta = 0.5 * (point_fieldsNT.sigma[0] + point_fieldsPT_or_NTPIPtr->sigma[0]);
			point_fieldsNT.v_r = point_fieldsPT_or_NTPIPtr->v_r + delt * (ring_p - sigma_theta / rho / g_domain->ring_R);
#endif
		}
		else
		{
			setValue(point_fieldsNT.u, 0.0);
#if RING_PROBLEM
			point_fieldsNT.v_r = 0.0;
#endif
		}
	}
	else
	{
		Compute_Bulk_vseps_IC(point_fieldsNT);
	}
	// source updated
#if HAVE_SOURCE
	g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(point_fieldsNT.x, point_fieldsNT.time, flag, point_fieldsNT.source_v, point_fieldsNT.source_sigma);
#if RING_PROBLEM
	double sigma_theta_source_final = E_iso * point_fieldsNT.v_r / g_domain->ring_R; // E v_r / R (same on both sides)
	point_fieldsNT.source_sigma[0] += sigma_theta_source_final;
#endif
#endif
}

SL_Elastic_InterfaceProperties::SL_Elastic_InterfaceProperties()
{
	bulk_leftPtr = NULL;
	bulk_rightPtr= NULL;
	bulk_inside_domain4BoundaryPtr = NULL;
	for (int i = 0; i < DiM; ++i)
		directionalBCType[i] = bct_Unspecified;
	bulk_leftPtr_deletable = true;
	bulk_rightPtr_deletable = true;
}

SL_Elastic_InterfaceProperties::~SL_Elastic_InterfaceProperties()
{
	if ((bulk_leftPtr != NULL) && (bulk_leftPtr_deletable))
		delete bulk_leftPtr;
	if ((bulk_rightPtr != NULL) && (bulk_rightPtr_deletable))
		delete bulk_rightPtr;
}

SL_Elastic_InterfaceProperties::SL_Elastic_InterfaceProperties(const SL_Elastic_InterfaceProperties& other)
{
	bulk_leftPtr = NULL;
	bulk_rightPtr = NULL;
	bulk_inside_domain4BoundaryPtr = NULL;
	for (int i = 0; i < DiM; ++i)
		directionalBCType[i] = bct_Unspecified;
	bulk_leftPtr_deletable = true;
	bulk_rightPtr_deletable = true;
	(*this) = other;
}

SL_Elastic_InterfaceProperties& SL_Elastic_InterfaceProperties::operator=(const SL_Elastic_InterfaceProperties& other)
{
	interfaceLoc = other.interfaceLoc;
	isIsoInterface = other.isIsoInterface;

	for (unsigned int i = 0; i < DiM; ++i)
		directionalBCType[i] = other.directionalBCType[i];
	
	if (bulk_leftPtr_deletable)
	{
		if (bulk_leftPtr != NULL)
			delete bulk_leftPtr;
	}
	bulk_leftPtr = NULL;
	if (bulk_leftPtr_deletable)
	{
		if (bulk_leftPtr != NULL)
			delete bulk_leftPtr;
	}
	bulk_leftPtr = NULL;
	bulk_leftPtr_deletable = other.bulk_leftPtr_deletable;
	if (bulk_leftPtr_deletable)
	{
		if (other.bulk_leftPtr != NULL)
		{
			bulk_leftPtr = new SL_Bulk_Properties();
			(*bulk_leftPtr) = *other.bulk_leftPtr;
		}
	}
	else
		bulk_leftPtr= other.bulk_leftPtr;

	bulk_rightPtr_deletable = other.bulk_rightPtr_deletable;
	if (bulk_rightPtr_deletable)
	{
		if (other.bulk_rightPtr != NULL)
		{
			bulk_rightPtr = new SL_Bulk_Properties();
			(*bulk_rightPtr) = *other.bulk_rightPtr;
		}
	}
	else
		bulk_rightPtr = other.bulk_rightPtr;

	bulk_inside_domain4BoundaryPtr = NULL;
	if (bulk_leftPtr == NULL) // left of the domain
	{
		interfaceLoc = ilt_left;
		for (int i = 0; i < DiM; ++i)
			if (directionalBCType[i] == bct_Unspecified)
				directionalBCType[i] = bct_Neumann;
		bulk_inside_domain4BoundaryPtr = bulk_rightPtr;
	}
	if (bulk_rightPtr == NULL) // right of the domain
	{
		interfaceLoc = ilt_right;
		for (int i = 0; i < DiM; ++i)
			if (directionalBCType[i] == bct_Unspecified)
				directionalBCType[i] = bct_Neumann;
		bulk_inside_domain4BoundaryPtr = bulk_leftPtr;
	}

	CopyVec(other.iso_wlSide_rGoing_2_sigmaI, iso_wlSide_rGoing_2_sigmaI);
	CopyVec(other.iso_wrSide_lGoing_2_sigmaI, iso_wrSide_lGoing_2_sigmaI);

	CopyVec(other.iso_wlSide_rGoing_2_vI, iso_wlSide_rGoing_2_vI);
	CopyVec(other.iso_wrSide_lGoing_2_vI, iso_wrSide_lGoing_2_vI);

	CopyVec(other.iso_sigmaStar_2_vStar_lSide, iso_sigmaStar_2_vStar_lSide);
	CopyVec(other.iso_sigmaStar_2_vStar_rSide, iso_sigmaStar_2_vStar_rSide);

	CopyVec(other.iso_wlSide_rGoing_2_vStar_lSide, iso_wlSide_rGoing_2_vStar_lSide);
	CopyVec(other.iso_wrSide_lGoing_2_vStar_rSide, iso_wrSide_lGoing_2_vStar_rSide);

	YAve_n = other.YAve_n;
	YAve_t = other.YAve_t;

#if !USE_ISO_ASSUMPTION
	CopyMat(other.wlSide_rGoing_2_sigmaI, wlSide_rGoing_2_sigmaI);
	CopyMat(other.wrSide_lGoing_2_sigmaI, wrSide_lGoing_2_sigmaI);

	CopyMat(other.wlSide_rGoing_2_vI, wlSide_rGoing_2_vI);
	CopyMat(other.wrSide_lGoing_2_vI, wrSide_lGoing_2_vI);

	CopyMat(other.sigmaStar_2_vStar_lSide, sigmaStar_2_vStar_lSide);
	CopyMat(other.sigmaStar_2_vStar_rSide, sigmaStar_2_vStar_rSide);
	CopyMat(other.wlSide_rGoing_2_vStar_lSide, wlSide_rGoing_2_vStar_lSide);
	CopyMat(other.wrSide_lGoing_2_vStar_rSide, wrSide_lGoing_2_vStar_rSide);

	CopyMat(other.Yl_plus_Yr, Yl_plus_Yr);
	CopyMat(other.Zl_plus_Zr, Zl_plus_Zr);
	CopyMat(other.inv_Yl_plus_Yr, inv_Yl_plus_Yr);
	CopyMat(other.inv_Zl_plus_Zr, inv_Zl_plus_Zr);

#if DiM2a3_F
	sumYnn = other.sumYnn;
	CopyVec(other.sumYnt, sumYnt);
	CopyVec(other.sumYtn, sumYtn);

	CopyVec(other.gammaVec, gammaVec);
	CopyVec(other.betaVec, betaVec);

	CopyMat(other.sumYtt, sumYtt);
	CopyMat(other.sumYtt_eff, sumYtt_eff);
	CopyMat(other.alphaMat, alphaMat);
#endif

#endif
	return *this;
}

void SL_Elastic_InterfaceProperties::Read(istream& in)
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
		if (buf == "directionalBCType")
		{
			for (int i = 0; i < DiM; ++i)
				in >> directionalBCType[i];
		}
		else if (buf == "bulk_leftPtr")
		{
			bulk_leftPtr_deletable = true;
			bulk_leftPtr = new SL_Bulk_Properties();
			bulk_leftPtr->Read_SL_Bulk_Properties(in, 1);
		}
		else if (buf == "bulk_rightPtr")
		{
			bulk_rightPtr_deletable = true;
			bulk_rightPtr = new SL_Bulk_Properties();
			bulk_rightPtr->Read_SL_Bulk_Properties(in, 1);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	FormMaps();
}

void SL_Elastic_InterfaceProperties::Print(ostream& out, bool printComputedVals) const
{
	out << "{\n";
	out << "directionalBCType\n";	
	for (int i = 0; i < DiM; ++i)
		out << directionalBCType[i] << '\t';
	out << '\n';
	if ((bulk_leftPtr_deletable) && (bulk_leftPtr != NULL))
	{
		out << "bulk_leftPtr\n";
		bulk_leftPtr->Print(out, printComputedVals);
		out << '\n';
	}
	if ((bulk_rightPtr_deletable) && (bulk_rightPtr != NULL))
	{
		out << "bulk_rightPtr\n";
		bulk_rightPtr->Print(out, printComputedVals);
		out << '\n';
	}
	out << "}\n";
	if (!printComputedVals)
		return;
	out << "[\n";
	out << "isIsoInterface\t" << isIsoInterface << '\n';
	out << "interfaceLoc\t" << interfaceLoc << '\n';
	out << "bulk_leftPtr_deletable\t" << bulk_leftPtr_deletable << '\n';
	out << "bulk_rightPtr_deletable\t" << bulk_rightPtr_deletable << '\n';

	out << "YAve_n\t" << YAve_n << '\n';
	out << "YAve_t\t" << YAve_t << '\n';

	out << "iso_wlSide_rGoing_2_sigmaI\n";	PrintV(iso_wlSide_rGoing_2_sigmaI, out); out << '\n';
	out << "iso_wrSide_lGoing_2_sigmaI\n";	PrintV(iso_wrSide_lGoing_2_sigmaI, out); out << '\n';
	out << "iso_wlSide_rGoing_2_vI\n";	PrintV(iso_wlSide_rGoing_2_vI, out); out << '\n';
	out << "iso_wrSide_lGoing_2_vI\n";	PrintV(iso_wrSide_lGoing_2_vI, out); out << '\n';

	out << "iso_sigmaStar_2_vStar_lSide\n";	PrintV(iso_sigmaStar_2_vStar_lSide, out); out << '\n';
	out << "iso_sigmaStar_2_vStar_rSide\n";	PrintV(iso_sigmaStar_2_vStar_rSide, out); out << '\n';
	out << "iso_wlSide_rGoing_2_vStar_lSide\n";	PrintV(iso_wlSide_rGoing_2_vStar_lSide, out); out << '\n';
	out << "iso_wrSide_lGoing_2_vStar_rSide\n";	PrintV(iso_wrSide_lGoing_2_vStar_rSide, out); out << '\n';
#if !USE_ISO_ASSUMPTION
	out << "wlSide_rGoing_2_sigmaI\n";	PrintM(wlSide_rGoing_2_sigmaI, out); out << '\n';
	out << "wrSide_lGoing_2_sigmaI\n";	PrintM(wrSide_lGoing_2_sigmaI, out); out << '\n';
	out << "wlSide_rGoing_2_vI\n";	PrintM(wlSide_rGoing_2_vI, out); out << '\n';
	out << "wrSide_lGoing_2_vI\n";	PrintM(wrSide_lGoing_2_vI, out); out << '\n';

	out << "sigmaStar_2_vStar_lSide\n";	PrintM(sigmaStar_2_vStar_lSide, out); out << '\n';
	out << "sigmaStar_2_vStar_rSide\n";	PrintM(sigmaStar_2_vStar_rSide, out); out << '\n';
	out << "wlSide_rGoing_2_vStar_lSide\n";	PrintM(wlSide_rGoing_2_vStar_lSide, out); out << '\n';
	out << "wrSide_lGoing_2_vStar_rSide\n";	PrintM(wrSide_lGoing_2_vStar_rSide, out); out << '\n';

	out << "Yl_plus_Yr\n";	PrintM(Yl_plus_Yr, out); out << '\n';
	out << "Zl_plus_Zr\n";	PrintM(Zl_plus_Zr, out); out << '\n';
	out << "inv_Yl_plus_Yr\n";	PrintM(inv_Yl_plus_Yr, out); out << '\n';
	out << "inv_Zl_plus_Zr\n";	PrintM(inv_Zl_plus_Zr, out); out << '\n';
#if DiM2a3_F
	out << "sumYnn\t" << sumYnn << '\n';
	out << "sumYnt\n";	PrintVm1(sumYnt, out); out << '\n';
	out << "sumYtn\n";	PrintVm1(sumYtn, out); out << '\n';
	out << "sumYtt\n";	PrintMm1(sumYtt, out); out << '\n';
	out << "sumYtt_eff\n";	PrintMm1(sumYtt_eff, out); out << '\n';
	out << "alphaMat\n";	PrintMm1(alphaMat, out); out << '\n';
	out << "gammaVec\n";	PrintVm1(gammaVec, out); out << '\n';
	out << "betaVec\n";	PrintVm1(betaVec, out); out << '\n';
#endif
#endif
	out << "]\n";
}

void SL_Elastic_InterfaceProperties::FormMaps()
{
	interfaceLoc = ilt_twoSided;
	if (bulk_leftPtr == NULL) // left of the domain
	{
		interfaceLoc = ilt_left;
		for (int i = 0; i < DiM; ++i)
			if (directionalBCType[i] == bct_Unspecified)
				directionalBCType[i] = bct_Neumann;
		bulk_inside_domain4BoundaryPtr = bulk_rightPtr;
	}
	if (bulk_rightPtr == NULL) // right of the domain
	{
		interfaceLoc = ilt_right;
		for (int i = 0; i < DiM; ++i)
			if (directionalBCType[i] == bct_Unspecified)
				directionalBCType[i] = bct_Neumann;
		bulk_inside_domain4BoundaryPtr = bulk_leftPtr;
	}
	if (interfaceLoc != ilt_twoSided)
	{
		isIsoInterface = bulk_inside_domain4BoundaryPtr->is_iso;
		YAve_n = bulk_inside_domain4BoundaryPtr->inv_c_rhos[0];
		YAve_t = YAve_n;
#if DiM2a3_F
		YAve_t = bulk_inside_domain4BoundaryPtr->inv_c_rhos[1];
#endif
		return;
	}

	// the rest is for two-sided interface
	isIsoInterface = bulk_leftPtr->is_iso && bulk_rightPtr->is_iso;
	if (isIsoInterface)
	{
		double Zl, Zr, Yl, Yr, inv_Zlpr, inv_Ylpr;
		for (int i = 0; i < DiM; ++i)
		{
			Zl = bulk_leftPtr->c_rhos[i];
			Zr = bulk_rightPtr->c_rhos[i];

			Yl = bulk_leftPtr->inv_c_rhos[i];
			Yr = bulk_rightPtr->inv_c_rhos[i];
			inv_Zlpr = 1.0 / (Zl + Zr);
			inv_Ylpr = 1.0 / (Yl + Yr);
			iso_wlSide_rGoing_2_sigmaI[i] = inv_Ylpr * Yl;
			iso_wrSide_lGoing_2_sigmaI[i] = inv_Ylpr * Yr;
			iso_wlSide_rGoing_2_vI[i] = -inv_Zlpr;
			iso_wrSide_lGoing_2_vI[i] =  inv_Zlpr;

			iso_sigmaStar_2_vStar_lSide[i] = Yl;
			iso_sigmaStar_2_vStar_rSide[i] = -Yr;
			iso_wlSide_rGoing_2_vStar_lSide[i] = -Yl;
			iso_wrSide_lGoing_2_vStar_rSide[i] = Yr;
		}
		YAve_n = 0.5 * (bulk_leftPtr->inv_c_rhos[0] + bulk_rightPtr->inv_c_rhos[0]);
		YAve_t = YAve_n;
#if DiM2a3_F
		YAve_t = 0.5 * (bulk_leftPtr->inv_c_rhos[1] + bulk_rightPtr->inv_c_rhos[1]);
#endif
	}
#if !USE_ISO_ASSUMPTION
	for (int i = 0; i < DiM; ++i)
		for (int j = 0; j < DiM; ++j)
		{
			Yl_plus_Yr[i][j] = bulk_leftPtr->Y[i][j] + bulk_rightPtr->Y[i][j];
			Zl_plus_Zr[i][j] = bulk_leftPtr->Z[i][j] + bulk_rightPtr->Z[i][j];
		}

	YAve_n = 0.5 * Yl_plus_Yr[0][0];
	YAve_t = YAve_n;
#if DiM2a3_F
	YAve_t = 0.5 * Yl_plus_Yr[1][1];
#endif
	Inverse(Yl_plus_Yr, inv_Yl_plus_Yr);
	Inverse(Zl_plus_Zr, inv_Zl_plus_Zr);
	ProductMatMat(inv_Yl_plus_Yr, bulk_leftPtr->H, wlSide_rGoing_2_sigmaI);
	ProductMatMat(inv_Yl_plus_Yr, bulk_rightPtr->H, wrSide_lGoing_2_sigmaI);
	ProductMatMat(inv_Zl_plus_Zr, bulk_leftPtr->E, wlSide_rGoing_2_vI);
	FactorMat(wlSide_rGoing_2_vI, -1.0);
	ProductMatMat(inv_Zl_plus_Zr, bulk_rightPtr->E, wrSide_lGoing_2_vI);

	CopyMat(bulk_leftPtr->Y, sigmaStar_2_vStar_lSide);
	CopyMat_withFactor(bulk_rightPtr->Y, sigmaStar_2_vStar_rSide, -1.0);
	CopyMat_withFactor(bulk_leftPtr->H, wlSide_rGoing_2_vStar_lSide, -1.0);
	CopyMat(bulk_rightPtr->H, wrSide_lGoing_2_vStar_rSide);

#if DiM2a3_F
	Break_nt_Matrix_n_t_components(true, Yl_plus_Yr, sumYnn, sumYnt, sumYtn, sumYtt, sumYtt_eff, alphaMat);
	double inv_sumYnn = 1.0 / sumYnn;
	for (int i = 0; i < DiMm1; ++i)
	{
		gammaVec[i] = inv_sumYnn * sumYnt[i];
		betaVec[i] = 0.0;
	}
	for (int i = 0; i < DiMm1; ++i)
		for (int j = 0; j < DiMm1; ++j)
			betaVec[j] += gammaVec[i] * alphaMat[i][j];
#endif
#endif
}

void SL_Elastic_InterfaceProperties::Compute_sigmaI_vI_from_ws(const VEC& wlSide_rGoing, const VEC& wrSide_lGoing, VEC& vI, VEC& sigmaI)
{
	if (interfaceLoc != ilt_twoSided)
	{
		Compute_sigmaI_vI_from_wsAndBC__BoundaryCase(wlSide_rGoing, wrSide_lGoing, vI, sigmaI);
		return;
	}
	if (isIsoInterface)
	{
		for (int i = 0; i < DiM; ++i)
		{
			vI[i] = wlSide_rGoing[i] * iso_wlSide_rGoing_2_vI[i] + wrSide_lGoing[i] * iso_wrSide_lGoing_2_vI[i];
			sigmaI[i] = wlSide_rGoing[i] * iso_wlSide_rGoing_2_sigmaI[i] + wrSide_lGoing[i] * iso_wrSide_lGoing_2_sigmaI[i];
		}
		return;
	}
#if !USE_ISO_ASSUMPTION
	double tmp;
	for (int i = 0; i < DiM; ++i)
	{
		tmp = 0.0;
		for (int j = 0; j < DiM; ++j)
			tmp += wlSide_rGoing_2_vI[i][j] * wlSide_rGoing[j] + wrSide_lGoing_2_vI[i][j] * wrSide_lGoing[j];
		vI[i] = tmp;
		tmp = 0.0;
		for (int j = 0; j < DiM; ++j)
			tmp += wlSide_rGoing_2_sigmaI[i][j] * wlSide_rGoing[j] + wrSide_lGoing_2_sigmaI[i][j] * wrSide_lGoing[j];
		sigmaI[i] = tmp;
	}
#endif
}

void SL_Elastic_InterfaceProperties::Compute_sigmaI_vI_from_wsAndBC__BoundaryCase(const VEC& w_OR_BC_lSide, const VEC& w_OR_BC_rSide, VEC& vI, VEC& sigmaI)
{
	if (isIsoInterface)
	{
		BoundaryConditionT bcT;
		double bcVal;

		if (interfaceLoc == ilt_left)
		{
			for (int i = 0; i < DiM; ++i)
			{
				bcT = directionalBCType[i];
				if (bcT == bct_Dirichlet)
				{
					bcVal = w_OR_BC_lSide[i];
					vI[i] = bcVal;
					sigmaI[i] = w_OR_BC_rSide[i] - bulk_inside_domain4BoundaryPtr->c_rhos[i] * bcVal;
				}
				else if (bcT == bct_Neumann)
				{
					bcVal = w_OR_BC_lSide[i];
					sigmaI[i] = bcVal;
					vI[i] = bulk_inside_domain4BoundaryPtr->inv_c_rhos[i] * (w_OR_BC_rSide[i] - bcVal);
				}
				else if (bcT == bct_Characteristics) // characteristics is comping in
				{
					sigmaI[i] = 0.5 * (w_OR_BC_rSide[i] + w_OR_BC_lSide[i]);
					vI[i] = 0.5 * bulk_inside_domain4BoundaryPtr->inv_c_rhos[i] * (w_OR_BC_rSide[i] - w_OR_BC_lSide[i]);
				}
				else
				{
					cout << bcT << '\n';
					THROW("Invalid boundary flag type\n");
				}
			}
		}
		else if (interfaceLoc == ilt_right)
		{
			for (int i = 0; i < DiM; ++i)
			{
				bcT = directionalBCType[i];
				if (bcT == bct_Dirichlet)
				{
					bcVal = w_OR_BC_rSide[i];
					vI[i] = bcVal;
					sigmaI[i] = w_OR_BC_lSide[i] + bulk_inside_domain4BoundaryPtr->c_rhos[i] * bcVal;
				}
				else if (bcT == bct_Neumann)
				{
					bcVal = w_OR_BC_rSide[i];
					sigmaI[i] = bcVal;
					vI[i] = bulk_inside_domain4BoundaryPtr->inv_c_rhos[i] * (bcVal - w_OR_BC_lSide[i]);
				}
				else if (bcT == bct_Characteristics) // characteristics is comping in
				{
					sigmaI[i] = 0.5 * (w_OR_BC_rSide[i] + w_OR_BC_lSide[i]);
					vI[i] = 0.5 * bulk_inside_domain4BoundaryPtr->inv_c_rhos[i] * (w_OR_BC_rSide[i] - w_OR_BC_lSide[i]);
				}
				else
				{
					cout << bcT << '\n';
					THROW("Invalid boundary flag type\n");
				}
			}
		}
		else
		{
			cout << "interfaceLoc\t" << interfaceLoc << '\n';
			THROW("Invalid interface side type\n");
		}
		return;
	}
#if !USE_ISO_ASSUMPTION
	// anisotropic case
	// for this case only pure Dirichlet, Neumann, or characteristic BC is supported, specified by flag 0
	BoundaryConditionT bcT = directionalBCType[0];
	if (bcT == bct_Characteristics)
	{
		VEC half_wr_plus_wl, half_wr_minus_wl;
		for (int i = 0; i < DiM; ++i)
		{
			half_wr_plus_wl[i] = 0.5 * (w_OR_BC_rSide[i] + w_OR_BC_lSide[i]);
			half_wr_minus_wl[i] = w_OR_BC_rSide[i] - half_wr_plus_wl[i];
		}
		setValue(sigmaI, 0.0);
		setValue(vI, 0.0);
		// sigma* = 0.5 E^t (wr + wl);		v* = 0.5 H (wr - wl)
		for (int i = 0; i < DiM; ++i)
		{
			for (int j = 0; j < DiM; ++j)
			{
				vI[i] += bulk_inside_domain4BoundaryPtr->H[i][j] * half_wr_minus_wl[j];
				sigmaI[i] += bulk_inside_domain4BoundaryPtr->E[j][i] * half_wr_plus_wl[j];
			}
		}
		return;
	}

	if (interfaceLoc == ilt_left)
	{
		if (bcT == bct_Dirichlet)
		{
			for (int i = 0; i < DiM; ++i)
			{
				vI[i] = w_OR_BC_lSide[i];
				// sigmaI = (E^t w)r - Zr v*
				sigmaI[i] = 0.0;
				for (int j = 0; j < DiM; ++j)
					sigmaI[i] += (bulk_inside_domain4BoundaryPtr->E[j][i] * w_OR_BC_rSide[j] - bulk_inside_domain4BoundaryPtr->Z[i][j] * vI[j]);
			}
		}
		else if (bcT == bct_Neumann)
		{
			for (int i = 0; i < DiM; ++i)
			{
				sigmaI[i] = w_OR_BC_lSide[i];
				// vI = (Hw)r - Yr sigma*
				vI[i] = 0.0;
				for (int j = 0; j < DiM; ++j)
					vI[i] += (bulk_inside_domain4BoundaryPtr->H[i][j] * w_OR_BC_rSide[j] - bulk_inside_domain4BoundaryPtr->Y[i][j] * sigmaI[j]);
			}
		}
		else
		{
			cout << bcT << '\n';
			THROW("Invalid boundary flag type\n");
		}
	}
	else if (interfaceLoc == ilt_right)
	{
		if (bcT == bct_Dirichlet)
		{
			for (int i = 0; i < DiM; ++i)
			{
				vI[i] = w_OR_BC_rSide[i];
				// sigmaI = (E^t w)l + Zl v*
				sigmaI[i] = 0.0;
				for (int j = 0; j < DiM; ++j)
					sigmaI[i] += (bulk_inside_domain4BoundaryPtr->E[j][i] * w_OR_BC_lSide[j] + bulk_inside_domain4BoundaryPtr->Z[i][j] * vI[j]);
			}
		}
		else if (bcT == bct_Neumann)
		{
			for (int i = 0; i < DiM; ++i)
			{
				sigmaI[i] = w_OR_BC_rSide[i];
				// vI = (Hw)r - Yr sigma*
				vI[i] = 0.0;
				for (int j = 0; j < DiM; ++j)
					vI[i] -= (bulk_inside_domain4BoundaryPtr->H[i][j] * w_OR_BC_lSide[j] - bulk_inside_domain4BoundaryPtr->Y[i][j] * sigmaI[j]);
			}
		}
		else
		{
			cout << bcT << '\n';
			THROW("Invalid boundary flag type\n");
		}
	}
	else
	{
		cout << "interfaceLoc\t" << interfaceLoc << '\n';
		THROW("Invalid interface side type\n");
	}
#endif
}

void SL_Elastic_InterfaceProperties::Compute_vStars_from_sigmaStar_ws(const VEC& wlSide_rGoing, const VEC& wrSide_lGoing, const VEC& sigmaStar, VEC& vStarLeft, VEC& vStarRight)
{
	if (isIsoInterface)
	{
		for (int i = 0; i < DiM; ++i)
		{
			vStarLeft[i] = iso_sigmaStar_2_vStar_lSide[i] * sigmaStar[i] + iso_wlSide_rGoing_2_vStar_lSide[i] * wlSide_rGoing[i];
			vStarRight[i] = iso_sigmaStar_2_vStar_rSide[i] * sigmaStar[i] + iso_wrSide_lGoing_2_vStar_rSide[i] * wrSide_lGoing[i];
		}
		return;
	}
#if !USE_ISO_ASSUMPTION
	double tmp;
	for (int i = 0; i < DiM; ++i)
	{
		tmp = 0.0;
		for (int j = 0; j < DiM; ++j)
			tmp += sigmaStar_2_vStar_lSide[i][j] * sigmaStar[j] + wlSide_rGoing_2_vStar_lSide[i][j] * wlSide_rGoing[j];
		vStarLeft[i] = tmp;
		tmp = 0.0;
		for (int j = 0; j < DiM; ++j)
			tmp += sigmaStar_2_vStar_rSide[i][j] * sigmaStar[j] + wrSide_lGoing_2_vStar_rSide[i][j] * wrSide_lGoing[j];
		vStarRight[i] = tmp;
	}
#endif
}

void SL_Elastic_InterfaceProperties::Compute_Downstream_characteristic_From_in_situ_Stress(bool has_in_situ, const VEC& in_situ_stress, const VEC& wlSide_rGoing_WO_in_situ, const VEC& wrSide_lGoing_WO_in_situ,
	VEC& wlSide_rGoing, VEC& wrSide_lGoing)
{
	if (!has_in_situ)
	{
		CopyVec(wlSide_rGoing_WO_in_situ, wlSide_rGoing);
		CopyVec(wrSide_lGoing_WO_in_situ, wrSide_lGoing);
		return;
	}
	bulk_leftPtr->Compute_Downstream_characteristic_From_in_situ_Stress(in_situ_stress, wlSide_rGoing_WO_in_situ, wlSide_rGoing);
	bulk_rightPtr->Compute_Downstream_characteristic_From_in_situ_Stress(in_situ_stress, wrSide_lGoing_WO_in_situ, wrSide_lGoing);
}

void SL_Elastic_InterfaceProperties::Compute_stresses_from_vs_and_characteristicss(const VEC& wlSide_rGoing, const VEC& wrSide_lGoing, const VEC* vlPtr, const VEC* vrPtr,
	VEC& sigmal, VEC& sigmar)
{
	if ((bulk_leftPtr != NULL) && (vlPtr != NULL))
		bulk_leftPtr->Compute_stress_from_v_and_characteristics(wlSide_rGoing, *vlPtr, true, sigmal);
	if ((bulk_rightPtr != NULL) && (vrPtr != NULL))
		bulk_rightPtr->Compute_stress_from_v_and_characteristics(wrSide_lGoing, *vrPtr, false, sigmar);
}

void SL_Elastic_InterfaceProperties::Computecharacteristics_from_Traces_vel_stress_etc(const VEC* traceCurrentStress_WO_in_situ_LeftIn, const VEC* traceCurrentStress_WO_in_situ_RightIn, const VEC* traceCurrentVelLeftIn, const VEC* traceCurrentVelRightIn, VEC& wlSide_rGoing, VEC& wrSide_lGoing)
{
	if ((bulk_leftPtr != NULL) && (traceCurrentStress_WO_in_situ_LeftIn != NULL))
		bulk_leftPtr->q_to_w_vec(*traceCurrentVelLeftIn, *traceCurrentStress_WO_in_situ_LeftIn, wlSide_rGoing, true);
	if ((bulk_rightPtr != NULL) && (traceCurrentStress_WO_in_situ_RightIn != NULL))
		bulk_rightPtr->q_to_w_vec(*traceCurrentVelRightIn, *traceCurrentStress_WO_in_situ_RightIn, wrSide_lGoing, false);
}

void SL_Elastic_InterfaceProperties::Get_Bulks2Side(SL_Bulk_Properties*& bulk_left_outPtr, SL_Bulk_Properties*& bulk_right_outPtr)
{
	if (bulk_leftPtr != NULL)
	{
		bulk_left_outPtr = bulk_leftPtr;
		if (bulk_rightPtr != NULL)
			bulk_right_outPtr = bulk_rightPtr;
		else
			bulk_right_outPtr = bulk_leftPtr;
	}
	else
	{
		bulk_left_outPtr = bulk_rightPtr;
		bulk_right_outPtr = bulk_rightPtr;
	}
}
