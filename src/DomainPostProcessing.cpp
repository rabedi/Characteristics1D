#include "DomainPostProcessing.h"
#include "globalFunctions.h"
#include "LAfuncsFinalStep.h"
#include "SLDescriptorData.h"
#include "Domain_AllInterfacesAllTimes.h"


Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum()
{
	Subdomain_spatial_ave_sum_set_zero();
}

void Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_set_zero()
{
	timeIndex = 0;
	timeVal = 0.0;
	setValue(eps_ave_bc, 0.0);
	setValue(sigma_ave_bc, 0.0);

	setValue(u_L, 0.0);
	setValue(u_R, 0.0);
	setValue(sigman_L, 0.0);
	setValue(sigman_R, 0.0);
	setValue(impulse_L, 0.0);
	setValue(impulse_R, 0.0);
	setValue(impulse_BC, 0.0);

	power_L = 0.0;
	power_R = 0.0;
	energy_L = 0.0;
	energy_R = 0.0;
	energy_BC = 0.0;

	setValue(powerIDiss_vec, 0.0);
	setValue(eneIDiss_vec, 0.0);
	energyIDiss_Recoverable = 0.0;
	powerIDiss = 0.0;
	eneIDiss = 0.0;
	min_interface_damage = 1.0; 
	mean_interface_damage = 0.0;
	sdiv_interface_damage = 0.0;
	max_interface_damage = 0.0;
	max_interface_damage_source = 0.0;
	max2EffDelU = 0.0;

	setValue(u_ave, 0.0);
	setValue(v_ave, 0.0);
	setValue(eps_ave_bulk_blk , 0.0);
	setValue(eps_ave_bulk_intfc, 0.0);
	setValue(sigma_ave_bulk, 0.0);

	setValue(linMomentum_dsum, 0.0);
	K_dsum = 0.0, U_dsum = 0.0, phi_dsum = 0.0;

	source_e_dtsum = 0.0;
#if HAVE_SOURCE
#if RING_PROBLEM
	v_r_ave = 0.0;
#endif
	setValue(source_linMomentum_dsum, 0.0);
	source_e_v_dsum = 0.0, source_e_sigma_dsum = 0.0, source_e_dsum = 0.0;

	setValue(source_linMomentum_dtsum, 0.0);
	source_e_v_dsum = 0.0, source_e_sigma_dtsum = 0.0;

#if HAVE_SOURCE_ORDER0_q
	setValue(damping_linMomentum_dsum, 0.0);
	damping_e_v_dsum = 0.0, damping_e_sigma_dsum = 0.0, damping_e_dsum = 0.0;

	setValue(damping_linMomentum_dtsum, 0.0);
	damping_e_v_dsum = 0.0, damping_e_sigma_dtsum = 0.0, damping_e_dtsum = 0.0;
#endif
#endif
	energy_eps_total_bc = 0.0, energy_eps_diss_bc = 0.0, energy_eps_recoverable_bc = 0.0;
	energy_eps_total_bulk = 0.0, energy_eps_diss_bulk = 0.0, energy_eps_recoverable_bulk = 0.0;
}

void Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_Read_Data(istream& in)
{
	in >> timeIndex;
	in >> timeVal;
	/////////////////////////////////////////////////////////////////
	// energies
	/////////////
	// current energy phi = K + U at current tme
	in >> phi_dsum;
	in >> K_dsum;
	in >> U_dsum;
	double rat;
	in >> rat;

	///////////// energy: INPUT
	// Input total = phi_time=0 + input BC + input source
	in >> input_energy;

	// phi0
	double phi0;
	in >> phi0;

	// input BC = BC_L + BC_R
	in >> energy_BC;
	in >> energy_L;
	in >> energy_R;

	// input source = source_v + source_sigma
	in >> source_e_dtsum;
#if HAVE_SOURCE
	in >> source_e_v_dtsum;
	in >> source_e_sigma_dtsum;
#endif

	///////////// energy: DISSIPATION
	// phys_diss_tot = eneIDiss + damping_e_dtsum
	in >> phys_diss_tot;
	in >> phys_diss_lost;
	in >> energyIDiss_Recoverable;
	in >> eneIDiss;
	in >> phys_diss_interface_lost;
	
	// components of energy dissipation
	DiM2a3(for (int i = 0; i < DiM; ++i) in >> eneIDiss_vec[i]; );
	double energy_phys_diss_2_input;
	double energy_phys_diss_lost_2_input;
	double energyIDiss_Recoverable_2_input;
	in >> energy_phys_diss_lost_2_input;
	in >> energy_phys_diss_2_input;
	in >> energyIDiss_Recoverable_2_input;

#if HAVE_SOURCE_ORDER0_q
	// damping_e_dtsum = damping_e_v_dtsum + damping_e_sigma_dtsum 
	in >> damping_e_dtsum;
	in >> damping_e_v_dtsum;
	in >> damping_e_sigma_dtsum;
#endif

	///////////// energy: Numerical
	// input - phi_current - phys_diss_tot
	in >> numerial_energy_diss;

	///////////// numerical dissipation energy: normalized ones
	double energy_numerical_diss_2_input;
	double energy_numerical_diss_2_phys_diss;

	in >> energy_numerical_diss_2_input;
	in >> energy_numerical_diss_2_phys_diss;

	////////////////////////////////////////////////////////////////////////////////////////////
	// average strain, stress (and stiffness)

	// strain
	// strain BC
	for (int i = 0; i < DiM; ++i)
		in >> eps_ave_bc[i];
	// strain body
	for (int i = 0; i < DiM; ++i)
		in >> eps_ave_bulk[i];
	for (int i = 0; i < DiM; ++i)
		in >> eps_ave_bulk_intfc[i];
	for (int i = 0; i < DiM; ++i)
		in >> eps_ave_bulk_blk[i];

	// stress
	// stress BC
	for (int i = 0; i < DiM; ++i)
		in >> sigma_ave_bc[i];
	// stress body
	for (int i = 0; i < DiM; ++i)
		in >> sigma_ave_bulk[i];

	// sigmas of left and right
	for (int i = 0; i < DiM; ++i)
		in >> sigman_L[i];
	for (int i = 0; i < DiM; ++i)
		in >> sigman_R[i];

	double psi_diss_2_totall_bc, psi_diss_bc_2_energy_phys_diss_lost;
	in >> energy_eps_total_bc;
	in >> energy_eps_diss_bc;
	in >> energy_eps_recoverable_bc;
	in >> psi_diss_2_totall_bc;
	in >> psi_diss_bc_2_energy_phys_diss_lost;

	double psi_diss_2_totall_bulk, psi_diss_bulk_2_energy_phys_diss_lost;
	in >> energy_eps_total_bulk;
	in >> energy_eps_diss_bulk;
	in >> energy_eps_recoverable_bulk;
	in >> psi_diss_2_totall_bulk;
	in >> psi_diss_bulk_2_energy_phys_diss_lost;

	// strain bc - body
	for (int i = 0; i < DiM; ++i)
		in >> eps_ave_bc_minus_bulk[i];
	// stress bc - body
	for (int i = 0; i < DiM; ++i)
		in >> sigma_ave_bc_minus_bulk[i];

	// ratio relevant under certain loading and before fracture for 1D
#if DiM1
	double stress2strain_bc;
	double stress2strain_bulk;
	in >> stress2strain_bc;
	in >> stress2strain_bulk;
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	// damage stats
	in >> mean_interface_damage;
	in >> sdiv_interface_damage;
	in >> max_interface_damage;
	in >> min_interface_damage;
	in >> max_interface_damage_source;
	in >> max2EffDelU;

	/////////////////////////////////////////////////////////////////
	// linear momenta
	/////////////
	// current linear momentum at current tme
	for (int i = 0; i < DiM; ++i)
		in >> linMomentum_dsum[i];

	// impulses BC = BC_L + BC_R
	for (int i = 0; i < DiM; ++i)
		in >> impulse_BC[i];
	for (int i = 0; i < DiM; ++i)
		in >> impulse_L[i];
	for (int i = 0; i < DiM; ++i)
		in >> impulse_R[i];

	// source 
#if HAVE_SOURCE
	for (int i = 0; i < DiM; ++i)
		in >> source_linMomentum_dtsum[i];
#endif

	///////////// DISSIPATION part
	// nothing from the interfaces all from damping
#if HAVE_SOURCE_ORDER0_q
	for (int i = 0; i < DiM; ++i)
		in >> damping_linMomentum_dtsum[i];
#endif

	///////////// Numerical error
	for (int i = 0; i < DiM; ++i)
		in >> numerical_linmom_error[i];


	/////////////////////////////////////////////////////////////////
	// powers
	// BC = BC_L + BC_R
	double power_BC;
	in >> power_BC;
	in >> power_L;
	in >> power_R;

#if HAVE_SOURCE
	in >> source_e_dsum;
	in >> source_e_v_dsum;
	in >> source_e_sigma_dsum;
#endif

	///////////// power: DISSIPATION
	in >> powerIDiss;
	// components of energy dissipation
	DiM2a3(for (int i = 0; i < DiM; ++i) in >> powerIDiss_vec[i]; );

#if HAVE_SOURCE_ORDER0_q
	in >> damping_e_dsum;
	in >> damping_e_v_dsum;
	in >> damping_e_sigma_dsum;
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	// u and v
	for (int i = 0; i < DiM; ++i)
		in >> u_L[i];
	for (int i = 0; i < DiM; ++i)
		in >> u_R[i];
	for (int i = 0; i < DiM; ++i)
		in >> u_ave[i];
	for (int i = 0; i < DiM; ++i)
		in >> v_ave[i];
#if HAVE_SOURCE
#if RING_PROBLEM
	in >> v_r_ave;
#endif
#endif
}

void Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_Write_Data(ostream& out, double invLength)
{
	double input_energy2Use = input_energy;
	if (g_SL_desc_data.bndryLoad_inputEnergy > 0.0)
		input_energy2Use = g_SL_desc_data.bndryLoad_inputEnergy;
	double tol_input_ene = input_energy2Use * 1e-7;
	double tol_input_ene_pl = tol_input_ene * invLength;

	out << timeIndex;
	out << '\t' << timeVal;
	/////////////////////////////////////////////////////////////////
	// energies
	/////////////
	// current energy phi = K + U at current tme
	out << '\t' << phi_dsum;
	out << '\t' << K_dsum;
	out << '\t' << U_dsum;
	double rat = 0.0;
	if (phi_dsum > tol_input_ene)
		rat = U_dsum / phi_dsum;
	out << "\t" << rat;

	///////////// energy: INPUT
	// Input total = phi_time=0 + input BC + input source
	out << '\t' << input_energy;

	// phi0
	double phi0 = input_energy - energy_BC;
//#if HAVE_SOURCE
	phi0 -= source_e_dtsum;
//#endif
	out << '\t' << phi0;

	// input BC = BC_L + BC_R
	out << '\t' << energy_BC;
	out << '\t' << energy_L;
	out << '\t' << energy_R;

	// input source = source_v + source_sigma
	out << '\t' << source_e_dtsum;
#if HAVE_SOURCE
	out << '\t' << source_e_v_dtsum;
	out << '\t' << source_e_sigma_dtsum;
#endif

	///////////// energy: DISSIPATION
	// phys_diss_tot = eneIDiss +  damping_e_dtsum
	out << '\t' << phys_diss_tot;
	out << '\t' << phys_diss_lost;
	out << '\t' << energyIDiss_Recoverable;
	out << '\t' << eneIDiss;
	out << '\t' << phys_diss_interface_lost;

	// components of energy dissipation
	DiM2a3(for (int i = 0; i < DiM; ++i) out << '\t' << eneIDiss_vec[i]; );

	double energy_phys_diss_2_input = computeRatio(phys_diss_tot, input_energy2Use);
	double energy_phys_diss_lost_2_input = computeRatio(phys_diss_lost, input_energy2Use);
	double energyIDiss_Recoverable_2_input = computeRatio(energyIDiss_Recoverable, input_energy2Use);
	out << '\t' << energy_phys_diss_lost_2_input;
	out << '\t' << energy_phys_diss_2_input;
	out << '\t' << energyIDiss_Recoverable_2_input;
#if HAVE_SOURCE_ORDER0_q
	// damping_e_dtsum  = damping_e_v_dtsum + damping_e_sigma_dtsum 
	out << '\t' << damping_e_dtsum;
	out << '\t' << damping_e_v_dtsum;
	out << '\t' << damping_e_sigma_dtsum;
#endif

	///////////// energy: Numerical
	// input - phi_current - phys_diss_tot
	out << '\t' << numerial_energy_diss;

	///////////// numerical dissipation energy: normalized ones
	double energy_numerical_diss_2_input = computeRatio(numerial_energy_diss, input_energy2Use);
	double energy_numerical_diss_2_phys_diss = 0.0;
	if (phys_diss_tot > tol_input_ene)
		energy_numerical_diss_2_phys_diss = computeRatio(numerial_energy_diss, phys_diss_tot);

	out << '\t' << energy_numerical_diss_2_input;
	out << '\t' << energy_numerical_diss_2_phys_diss;

	////////////////////////////////////////////////////////////////////////////////////////////
	// average strain, stress (and stiffness)

	// strain
	// strain BC
	for (int i = 0; i < DiM; ++i)
		out << '\t' << eps_ave_bc[i];
	// strain body
	for (int i = 0; i < DiM; ++i)
		out << '\t' << eps_ave_bulk[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << eps_ave_bulk_intfc[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << eps_ave_bulk_blk[i];

	// stress
	// stress BC
	for (int i = 0; i < DiM; ++i)
		out << '\t' << sigma_ave_bc[i];
	// stress body
	for (int i = 0; i < DiM; ++i)
		out << '\t' << sigma_ave_bulk[i];

	// sigmas of left and right
	for (int i = 0; i < DiM; ++i)
		out << '\t' << sigman_L[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << sigman_R[i];

	double psi_diss_2_totall_bc = 0.0;
	double psi_diss_bc_2_energy_phys_diss_lost = 0.0;
	double denom = invLength * phys_diss_lost;
	if (energy_eps_total_bc > tol_input_ene_pl)
		psi_diss_2_totall_bc = computeRatio(energy_eps_diss_bc, energy_eps_total_bc);
	if (denom > tol_input_ene_pl)
		psi_diss_bc_2_energy_phys_diss_lost = computeRatio(energy_eps_diss_bc, denom);
	out << '\t' << energy_eps_total_bc;
	out << '\t' << energy_eps_diss_bc;
	out << '\t' << energy_eps_recoverable_bc;
	out << '\t' << psi_diss_2_totall_bc;
	out << '\t' << psi_diss_bc_2_energy_phys_diss_lost;

	double psi_diss_2_totall_bulk = 0.0;
	double psi_diss_bulk_2_energy_phys_diss_lost = 0.0;
	denom = invLength * phys_diss_lost;
	if (energy_eps_total_bulk > tol_input_ene_pl)
		psi_diss_2_totall_bulk = computeRatio(energy_eps_diss_bulk, energy_eps_total_bulk);
	if (denom > tol_input_ene_pl)
		psi_diss_bulk_2_energy_phys_diss_lost = computeRatio(energy_eps_diss_bulk, denom);
	out << '\t' << energy_eps_total_bulk;
	out << '\t' << energy_eps_diss_bulk;
	out << '\t' << energy_eps_recoverable_bulk;
	out << '\t' << psi_diss_2_totall_bulk;
	out << '\t' << psi_diss_bulk_2_energy_phys_diss_lost;

	// strain bc - body
	for (int i = 0; i < DiM; ++i)
		out << '\t' << eps_ave_bc_minus_bulk[i];
	// stress bc - body
	for (int i = 0; i < DiM; ++i)
		out << '\t' << sigma_ave_bc_minus_bulk[i];

	// ratio relevant under certain loading and before fracture for 1D
#if DiM1
	double stress2strain_bc = computeRatio(sigma_ave_bc[0], eps_ave_bc[0]);
	double stress2strain_bulk = computeRatio(sigma_ave_bulk[0], eps_ave_bulk [0]);
	out << '\t' << stress2strain_bc;
	out << '\t' << stress2strain_bulk;
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	// damage stats
	out << '\t' << mean_interface_damage;
	out << '\t' << sdiv_interface_damage;
	out << '\t' << max_interface_damage;
	out << '\t' << min_interface_damage;
	out << '\t' << max_interface_damage_source;
	out << '\t' << max2EffDelU;

	/////////////////////////////////////////////////////////////////
	// linear momenta
	/////////////
	// current linear momentum at current tme
	for (int i = 0; i < DiM; ++i)
		out << '\t' << linMomentum_dsum[i];

	// impulses BC = BC_L + BC_R
	for (int i = 0; i < DiM; ++i)
		out << '\t' << impulse_BC[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << impulse_L[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << impulse_R[i];

	// source 
#if HAVE_SOURCE
	for (int i = 0; i < DiM; ++i)
		out << '\t' << source_linMomentum_dtsum[i];
#endif

	///////////// DISSIPATION part
	// nothing from the interfaces all from damping
#if HAVE_SOURCE_ORDER0_q
	for (int i = 0; i < DiM; ++i)
		out << '\t' << damping_linMomentum_dtsum[i];
#endif

	///////////// Numerical error
	for (int i = 0; i < DiM; ++i)
		out << '\t' << numerical_linmom_error[i];


	/////////////////////////////////////////////////////////////////
	// powers
	// BC = BC_L + BC_R
	double power_BC = power_L + power_R;
	out << '\t' << power_BC;
	out << '\t' << power_L;
	out << '\t' << power_R;

#if HAVE_SOURCE
	out << '\t' << source_e_dsum;
	out << '\t' << source_e_v_dsum;
	out << '\t' << source_e_sigma_dsum;
#endif

	///////////// power: DISSIPATION
	out << '\t' << powerIDiss;
	// components of energy dissipation
	DiM2a3(for (int i = 0; i < DiM; ++i) out << '\t' << powerIDiss_vec[i]; );

#if HAVE_SOURCE_ORDER0_q
	out << '\t' << damping_e_dsum;
	out << '\t' << damping_e_v_dsum;
	out << '\t' << damping_e_sigma_dsum;
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	// u and v
	for (int i = 0; i < DiM; ++i)
		out << '\t' << u_L[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << u_R[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << u_ave[i];
	for (int i = 0; i < DiM; ++i)
		out << '\t' << v_ave[i];
#if HAVE_SOURCE
#if RING_PROBLEM
	out << '\t' << v_r_ave;
#endif
#endif
	out << '\n';
}

void Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_Data_Write_Header(ostream& out)
{
	vector<HeaderLabels> hLabels;
	Subdomain_spatial_ave_sum_Data_Get_Header_Labels(hLabels);
	PrintHeader(hLabels, out);
}

void Subdomain_spatial_ave_sum::Size_interface_sigma_delu_vec(unsigned int numInterfaces)
{
	interface_sigma_delu_vec.resize(numInterfaces);
	for (unsigned int i = 0; i < numInterfaces; ++i)
	{
		setValue(interface_sigma_delu_vec[i].first, 0.0);
		setValue(interface_sigma_delu_vec[i].second, 0.0);
	}
}

void Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_Data_Get_Header_Labels(vector<HeaderLabels>& hLabels)
{
	HeaderLabels hLabel;
	hLabel.group = "index";
	hLabel.subgroup = "time";
	hLabel.textLabel = "timeIndex";
	hLabel.latexLabel = "ti";
	hLabels.push_back(hLabel);

	hLabel.group = "index";
	hLabel.subgroup = "time";
	hLabel.textLabel = "time";
	hLabel.latexLabel = "t";
	hLabels.push_back(hLabel);

	/////////////////////////////////////////////////////////////////
	// energies
	string group = "energy";
	hLabel.group = group;

	hLabel.subgroup = "tT";

	hLabel.textLabel = "phi";
	hLabel.latexLabel = "\\phi";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "K";
	hLabel.latexLabel = "K";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "U";
	hLabel.latexLabel = "U";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "U2phi";
	hLabel.latexLabel = "U2phi";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "input";
	hLabel.textLabel = "input_energy";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{inp}}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "t0";
	hLabel.textLabel = "phi0";
	hLabel.latexLabel = "\\phi_0";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "BC";
	hLabel.textLabel = "energy_BC";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{BC}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_L";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{L}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_R";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{R}}";
	hLabels.push_back(hLabel);

	// input source = source_v + source_sigma
	hLabel.subgroup = "source";
	hLabel.textLabel = "source_e";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{src}}";
	hLabels.push_back(hLabel);
#if HAVE_SOURCE
	hLabel.textLabel = "source_e_v";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{src},v}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "source_e_sigma";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{src},\\sigma}";
	hLabels.push_back(hLabel);
#endif

	///////////// energy: DISSIPATION
	hLabel.subgroup = "dissipation";
	hLabel.textLabel = "phys_diss_tot";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,tot}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "phys_diss_lost";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,lost}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "phys_diss_recoverable";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,recov}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "eneIDiss";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,I}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "phys_diss_interface_lost";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,lost,I}}";
	hLabels.push_back(hLabel);

#if !DiM1
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "eneIDiss" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\mathcal{E}_{\\mathrm{diss,I" + pINDs[i] + "}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
#endif

	hLabel.textLabel = "energy_phys_diss_2_input";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,tot}}/\\mathcal{E}_{\\mathrm{inp}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_phys_diss_lost_2_input";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,lost}}/\\mathcal{E}_{\\mathrm{inp}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energyIDiss_Recoverable_2_input";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{diss,recov}}/\\mathcal{E}_{\\mathrm{inp}}";
	hLabels.push_back(hLabel);

#if HAVE_SOURCE_ORDER0_q
	hLabel.subgroup = "damping";
	hLabel.textLabel = "damping_e";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{dmp}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "damping_e_v";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{dmp},v}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "damping_e_sigma";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{dmp},\\sigma}";
	hLabels.push_back(hLabel);
#endif

	hLabel.subgroup = "numerical";
	hLabel.textLabel = "numerial_energy_diss";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{num}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_numerical_diss_2_input";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{num}}/\\mathcal{E}_{\\mathrm{inp}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_numerical_diss_2_phys_diss";
	hLabel.latexLabel = "\\mathcal{E}_{\\mathrm{num}}/\\mathcal{E}_{\\mathrm{diss,tot}}";
	hLabels.push_back(hLabel);

	////////////////////////////////////////////////////////////////////////////////////////////
	// strain / stress
	group = "strs_strn";
	hLabel.group = group;
	hLabel.subgroup = "strn";

	// strain BC
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "eps_bc" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\bar{\\epsilon}_{" + pINDcommas[i] + "\\mathrm{bc}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// strain body
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "eps_bulk" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "\\bar{\\epsilon}";
#else
		tmps = "\\bar{\\epsilon}_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// strain body - interface part
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "eps_bulk_intfc" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "\\bar{\\epsilon}_{\\mathrm{I}}";
#else
		tmps = "\\bar{\\epsilon}_{\\mathrm{I}," + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// strain body - bulk part
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "eps_bulk_blk" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "\\bar{\\epsilon}_{\\mathrm{b}}";
#else
		tmps = "\\bar{\\epsilon}_{\\mathrm{b}," + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}


	hLabel.subgroup = "strs";
	// stress BC
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "sigma_bc" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\bar{\\sigma}_{" + pINDcommas[i] + "\\mathrm{bc}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// stress body
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "sigma_bulk" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "\\bar{\\sigma}";
#else
		tmps = "\\bar{\\sigma}_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// sigmas of left and right
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "sigma_L" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\bar{\\sigma}_{" + pINDcommas[i] + "\\mathrm{L}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "sigma_R" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\bar{\\sigma}_{" + pINDcommas[i] + "\\mathrm{R}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	hLabel.subgroup = "energy";
	hLabel.textLabel = "energy_eps_total_bc";
	hLabel.latexLabel = "\\bar{\\psi}_{\\mathrm{bc}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_eps_diss_bc";
	hLabel.latexLabel = "\\bar{\\psi}_{\\mathrm{diss,bc}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_eps_recoverable_bc";
	hLabel.latexLabel = "\\bar{\\psi}_{\\mathrm{recov,bc}}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "damage";
	hLabel.textLabel = "energy_eps_diss2total_bc";
	hLabel.latexLabel = "{\\bar{\\psi}_{\\mathrm{diss,bc}}}/\\bar{\\psi}_{\\mathrm{bc}}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "ratio";
	hLabel.textLabel = "energy_eps_dissL2phid_bc";
	hLabel.latexLabel = "{L\\bar{\\psi}_{\\mathrm{diss,bc}}}/\\mathcal{E}_{\\mathrm{diss,lost}}";
	hLabels.push_back(hLabel);

	////////////////////////////////////////
	hLabel.subgroup = "energy";
	hLabel.textLabel = "energy_eps_total_bulk";
	hLabel.latexLabel = "\\bar{\\psi}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_eps_diss_bulk";
	hLabel.latexLabel = "\\bar{\\psi}_{\\mathrm{diss}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "energy_eps_recoverable_bulk";
	hLabel.latexLabel = "\\bar{\\psi}_{\\mathrm{recov}}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "damage";
	hLabel.textLabel = "energy_eps_diss2total_bulk";
	hLabel.latexLabel = "{\\bar{\\psi}_{\\mathrm{diss}}}/\\bar{\\psi}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "ratio";
	hLabel.textLabel = "energy_eps_dissL2phid_bulk";
	hLabel.latexLabel = "{L\\bar{\\psi}_{\\mathrm{diss}}}/\\mathcal{E}_{\\mathrm{diss,lost}}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "strn";
	// strain bc - body
	for (int i = 0; i < DiM; ++i)
	{
		string tmps2, tmps = "eps_bc_minus_bulk" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps2 = "\\bar{\\epsilon}_{" + pINDcommas[i] + "\\mathrm{bc}}";
#if DiM1
		tmps = tmps2 + "-\\bar{\\epsilon}";
#else
		tmps = tmps2 + "-\\bar{\\epsilon}_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	hLabel.subgroup = "strs";
	// stress bc - body
	for (int i = 0; i < DiM; ++i)
	{
		string tmps2, tmps = "sigma_bc_minus_bulk" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps2 = "\\bar{\\sigma}_{" + pINDcommas[i] + "\\mathrm{bc}}";
#if DiM1
		tmps = tmps2 + "-\\bar{\\sigma}";
#else
		tmps = tmps2 + "-\\bar{\\sigma}_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// ratio relevant under certain loading and before fracture for 1D
#if DiM1
	hLabel.subgroup = "strs2strn";
	hLabel.textLabel = "stress2strain_bc";
	hLabel.latexLabel = "\\bar{\\sigma}_{\\mathrm{bc}}/\\bar{\\epsilon}_{\\mathrm{bc}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "stress2strain_bulk";
	hLabel.latexLabel = "\\bar{\\sigma}/\\bar{\\epsilon}";
	hLabels.push_back(hLabel);
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	// damage stats
	group = "damage";
	hLabel.group = group;

	hLabel.subgroup = "mean";
	hLabel.textLabel = "mean_interface_damage";
	hLabel.latexLabel = "\\bar{D}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "sdiv";
	hLabel.textLabel = "sdiv_interface_damage";
	hLabel.latexLabel = "\\varsigma_{D}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "max";
	hLabel.textLabel = "max_interface_damage";
	hLabel.latexLabel = "D_\\mathrm{max}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "min";
	hLabel.textLabel = "min_interface_damage";
	hLabel.latexLabel = "D_\\mathrm{min}";
	hLabels.push_back(hLabel);

	hLabel.subgroup = "max";
	hLabel.textLabel = "max_interface_damage_source";
	hLabel.latexLabel = "D_\\mathrm{src,max}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "max2EffDelU";
	hLabel.latexLabel = "\\Delta{u_\\mathrm{n,max}}";
	hLabels.push_back(hLabel);

	/////////////////////////////////////////////////////////////////
	// linear momenta
	group = "linmomentum";
	hLabel.group = group;
	hLabel.subgroup = "tT";

	/////////////
	// current linear momentum at current tme
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "linMomentum" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "P";
#else
		tmps = "P_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// impulses BC = BC_L + BC_R
	hLabel.subgroup = "BC";
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "impulse_BC" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "J_{" + pINDcommas[i] + "\\mathrm{bc}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "impulse_L" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "J_{" + pINDcommas[i] + "\\mathrm{L}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "impulse_R" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "J_{" + pINDcommas[i] + "\\mathrm{R}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	// source 
	hLabel.subgroup = "source";
#if HAVE_SOURCE
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "source_linMomentum" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "P_{" + pINDcommas[i] + "\\mathrm{src}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
#endif

	///////////// DISSIPATION part
#if HAVE_SOURCE_ORDER0_q
	hLabel.subgroup = "damping";
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "damping_linMomentum" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "P_{" + pINDcommas[i] + "\\mathrm{dmp}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
#endif

	///////////// Numerical error
	hLabel.subgroup = "numerical";
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "numerical_linmom_error" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "P_{" + pINDcommas[i] + "\\mathrm{num}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}

	/////////////////////////////////////////////////////////////////
	// powers
	group = "power";
	hLabel.group = group;

	hLabel.subgroup = "BC";
	hLabel.textLabel = "power_BC";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{BC}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "power_L";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{L}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "power_R";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{R}}";
	hLabels.push_back(hLabel);

#if HAVE_SOURCE
	hLabel.subgroup = "source";

	hLabel.textLabel = "power_source_e";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{src}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "power_source_e_v";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{src},v}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "power_source_e_sigma";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{src},\\sigma}";
	hLabels.push_back(hLabel);
#endif

	///////////// power: DISSIPATION
	hLabel.subgroup = "dissipation";
	hLabel.textLabel = "powerIDiss";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{diss,I}}";
	hLabels.push_back(hLabel);

#if !DiM1
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "powerIDiss" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\mathcal{P}_{\\mathrm{diss,I" + pINDs[i] + "}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
#endif

#if HAVE_SOURCE_ORDER0_q
	hLabel.subgroup = "damping";
	hLabel.textLabel = "damping_p";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{dmp}}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "damping_p_v";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{dmp},v}";
	hLabels.push_back(hLabel);

	hLabel.textLabel = "damping_p_sigma";
	hLabel.latexLabel = "\\mathcal{P}_{\\mathrm{dmp},\\sigma}";
	hLabels.push_back(hLabel);
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	// u and v
	group = "u_v";
	hLabel.group = group;
	hLabel.subgroup = "u";
	// u left and right
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "u_L" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\bar{u}_{" + pINDcommas[i] + "\\mathrm{L}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "u_R" + pINDs[i];
		hLabel.textLabel = tmps;
		tmps = "\\bar{u}_{" + pINDcommas[i] + "\\mathrm{R}}";
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "u_bulk" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "\\bar{u}";
#else
		tmps = "\\bar{u}_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
	hLabel.subgroup = "v";
	for (int i = 0; i < DiM; ++i)
	{
		string tmps = "v_bulk" + pINDs[i];
		hLabel.textLabel = tmps;
#if DiM1
		tmps = "\\bar{v}";
#else
		tmps = "\\bar{v}_{" + pINDs[i] + "}";
#endif
		hLabel.latexLabel = tmps;
		hLabels.push_back(hLabel);
	}
#if HAVE_SOURCE
#if RING_PROBLEM
	hLabel.textLabel = "v_r";
	hLabel.latexLabel = "\\bar{v}_r";
	hLabels.push_back(hLabel);
#endif
#endif
}

OneBulktwoSideInterfaceInfo::OneBulktwoSideInterfaceInfo()
{
	bulkPtr = NULL;
	interfaceLeftOfBulkPtr = NULL;
	interfaceRightOfBulkPtr = NULL;;

	bulk_cntr = 0;
	interfaceLeftOfBulk_cntr = 0, interfaceRightOfBulk_cntr = 0;
}

OneSubdomain_All_bulksConnectivityInfo::OneSubdomain_All_bulksConnectivityInfo()
{
	Zero1D_Averages();
}

void OneSubdomain_All_bulksConnectivityInfo::Zero1D_Averages()
{
	rhoAve = 0.0;
	EAve = 0.0;
	cHarmonicAve = 0.0;
	cAve = 0.0;
	ZAve = 0.0;
	c_fromAverages = 0.0;
	Z_fromAverages = 0.0;
	DvvAve = 0.0;
	EHarmonicAve = 0.0;
}

void OneSubdomain_All_bulksConnectivityInfo::PrintIndicesLengthsKeyRunParameters(ostream& out) const
{
	out << "lengths\t" << segment_lengths.size() << '\n';
	for (unsigned int i = 0; i < segment_lengths.size(); ++i)
		out << '\t' << segment_lengths[i];
	out << "\ndeltaCs\t" << subdomain_interfaces.size() << '\n';
	for (unsigned int i = 0; i < subdomain_interfaces.size(); ++i)
		out << '\t' << subdomain_interfaces[i]->getDeltaC();
	out << "\nsigmaCs\t" << subdomain_interfaces.size() << '\n';
	for (unsigned int i = 0; i < subdomain_interfaces.size(); ++i)
		out << '\t' << subdomain_interfaces[i]->getSigmaC();
	out << "\nlength\t" << length << '\n';
	out << "xm\t" << xm << '\n';
	out << "xM\t" << xM << '\n';
	out << "numSegments\t" << numSegments << '\n';
	out << "bulk_cn_st\t" << bulk_cn_st << '\n';
	out << "bulk_cn_en\t" << bulk_cn_en << '\n';
	out << "xs\t" << subdomain_interface_xs.size() << '\n';
	for (unsigned int i = 0; i < subdomain_interface_xs.size(); ++i)
		out << '\t' << subdomain_interface_xs[i];
	out << "\n";
}


void OneSubdomain_All_bulksConnectivityInfo::OneSubdomain_All_bulksConnectivityInfo_Initialize()
{
	xm = subdomain_interface_xs[0];
	if (g_domain->isPeriodic)
		xm = g_domain->x_min;
	segment_lengths.resize(numSegments);
	double xStart = xm;

	length = 0.0;

	double segLength;
	for (int segi = 0; segi < numSegments; ++segi)
	{
		segLength = subdomain_interface_xs[segi + 1] - subdomain_interface_xs[segi];
		if (segLength < 0.0)
			segLength += g_domain->L;
		segment_lengths[segi] = segLength;
		xStart += segLength;
	}
	xM = xStart;
	length = xM - xm;
	inv_length = 1.0 / length;
}

void Subdomain_oneTime_spatial_points::Subdomain_oneTime_spatial_points_Print_Header_Data(ostream& out, bool printHeader)
{
	bool print_segi_pi = true, print_x = true, print_time = false;
	if (printHeader)
		OnePoint_inBulk_Fields::OnePoint_inBulk_Fields_Print_Header(out, print_segi_pi, print_x, print_time);
	unsigned int num_seg = spatialPoints.size();
	if (num_seg == 0)
		return;
	unsigned int num_pt = spatialPoints[0].size();
	for (unsigned int segi = 0; segi < num_seg; ++segi)
		for (unsigned int pti = 0; pti < num_pt; ++pti)
		{
			out << segi << '\t' << pti << '\t';
			spatialPoints[segi][pti].OnePoint_inBulk_Fields_Print_Data(out, print_x, print_time);
		}
}

Subdomain_spacetime_pp_data::Subdomain_spacetime_pp_data()
{
	print_space_points = true;
	uniform_delt = true;
	out_sd_summary = NULL;
	maxIter = 10;
	relTol4Conv = 1e-4;

	allSpatialPointsRevOrder_size = 3;
	allSpatialPointsRevOrder.resize(allSpatialPointsRevOrder_size);
	for (unsigned int i = 0; i < allSpatialPointsRevOrder.size(); ++i)
		allSpatialPointsRevOrder[i] = NULL;

	for (int tir = 0; tir < NUM_AVE_SUM_TIMES; ++tir)
		spatial_ave_sums[tir] = NULL;

	hasLeft = false;
	hasRight = false;
}

Subdomain_spacetime_pp_data::~Subdomain_spacetime_pp_data()
{
	if (out_sd_summary != NULL)
		delete out_sd_summary;
	for (unsigned int i = 0; i < allSpatialPointsRevOrder.size(); ++i)
	{
		if (allSpatialPointsRevOrder[i] != NULL)
			delete allSpatialPointsRevOrder[i];
	}
	for (int tir = 0; tir < NUM_AVE_SUM_TIMES; ++tir)
	{
		if (spatial_ave_sums[tir] != NULL)
			delete spatial_ave_sums[tir];
	}
}

void Subdomain_spacetime_pp_data::Initialize_Subdomain_spacetime_pp_data(unsigned int numTimesIn, OneSubdomain_All_bulksConnectivityInfo* sdciPtrIn, bool uniform_deltIn, int numTimeStep_BulkInterfacePoints_PrintIn, int numTimeStep_Interface_DSU_Fragment_PrintIn,
	int numSpatialSubsegments_BulkInterfacePoints_PrintIn, bool useRepeatedSimpsonRuleForHigherOrders,
	unsigned int maxIterIn, double relTol4ConvIn)
{
	numTimes = numTimesIn;
	uniform_delt = uniform_deltIn;
	numTimeStep_BulkInterfacePoints_Print = numTimeStep_BulkInterfacePoints_PrintIn;
	numTimeStep_Interface_DSU_Fragment_Print = numTimeStep_Interface_DSU_Fragment_PrintIn;
	sdciPtr = sdciPtrIn;
	print_space_points = (numTimeStep_BulkInterfacePoints_Print > 0);
	print_fragmentation = (numTimeStep_Interface_DSU_Fragment_Print > 0);

	numSpatialSubsegments_BulkInterfacePoints_Print = numSpatialSubsegments_BulkInterfacePoints_PrintIn;
	maxIter = maxIterIn;
	relTol4Conv = relTol4ConvIn;


	string fileName;
	string specificName = "_Summary";
	GetSubdomainIndexed_TimeIndexed_FileName(fileName, sdciPtr->subdomain_number, -1, specificName);
	out_sd_summary = new fstream();
	out_sd_summary->open(fileName.c_str(), ios::out);

	numSpatialPointsPerSegment = SetNewtonCotes_Points_AndWeights(numSpatialSubsegments_BulkInterfacePoints_Print, 
		useRepeatedSimpsonRuleForHigherOrders, spatialIntegrationWeights, spatialIntegrationPoints);
	Set_xs();
}

void Subdomain_spacetime_pp_data::AddComputeTimeStep(int timeIndexIn, double timeValue)
{
	if (spatial_ave_sums[NUM_AVE_SUM_TIMES - 1] != NULL)
		delete spatial_ave_sums[NUM_AVE_SUM_TIMES - 1];
	for (int tir = NUM_AVE_SUM_TIMES - 1; tir >= 1; --tir)
		spatial_ave_sums[tir] = spatial_ave_sums[tir - 1];
	spatial_ave_sums[0] = new Subdomain_spatial_ave_sum();

	unsigned szInterface = sdciPtr->subdomain_interfaces.size();
	spatial_ave_sums[0]->Size_interface_sigma_delu_vec(szInterface);

	times.push_back(timeValue);

	int lastIndex = allSpatialPointsRevOrder.size() - 1;
	if (allSpatialPointsRevOrder[lastIndex] != NULL)
	{
		delete allSpatialPointsRevOrder[lastIndex];
		allSpatialPointsRevOrder[lastIndex] = NULL;
	}
	for (int i = lastIndex; i > 0; --i)
		allSpatialPointsRevOrder[i] = allSpatialPointsRevOrder[i - 1];

	timeIndex = timeIndexIn;
	times[timeIndex] = timeValue;

	delt = 0.0;
	timeIntWeights_rti_size = 0;

	if (timeIndex >= 1)
	{
		delt = times[timeIndex] - times[timeIndex - 1];
		if (!uniform_delt || (timeIndex % 2 == 1))
		{
			timeIntWeights_rti_size = 2;
			timeIntWeights_rti.resize(timeIntWeights_rti_size);
			timeIntWeights_rti[0] = 0.5 * delt;
			timeIntWeights_rti[1] = timeIntWeights_rti[0];
		}
		else // difference of forward Euler and Simpson rule
		{
			timeIntWeights_rti_size = 3;
			timeIntWeights_rti.resize(timeIntWeights_rti_size);
			double tmp = delt / 6.0;
			timeIntWeights_rti[0] =  2.0 * tmp;
			timeIntWeights_rti[1] =  5.0 * tmp;
			timeIntWeights_rti[2] = -1.0 * tmp;
		}
	}

	// setting spatial domain
	Subdomain_oneTime_spatial_points* siotsp = new Subdomain_oneTime_spatial_points();
	allSpatialPointsRevOrder[0] = siotsp;
	siotsp->spatialPoints.resize(sdciPtr->numSegments);
	for (int segi = 0; segi < sdciPtr->numSegments; ++segi)
	{
		siotsp->spatialPoints[segi].resize(numSpatialPointsPerSegment);
		for (int pti = 0; pti < numSpatialPointsPerSegment; ++pti)
		{
			siotsp->spatialPoints[segi][pti].x = xs[segi][pti];
			siotsp->spatialPoints[segi][pti].time = timeValue;
		}
	}
	
	// Above -> initialization
	////////////////////////////////////////////////////////////////////
	// Below -> computations
	// A: inside points
	for (int segi = 0; segi < sdciPtr->numSegments; ++segi)
	{
		for (int pti = 1; pti < numSpatialPointsPerSegment - 1; ++pti)
			Compute_Inside_Segment_pt(segi, pti, *sdciPtr->subdomain_bulk_segments[segi].bulkPtr);
	}
	// B: Interfaces

	SL_OneInterfaceAllTimes* subdomain_interface;
	unsigned int interface_starti = 0;
	if (g_domain->isPeriodic)
		interface_starti = 1;
	for (unsigned int interfacei = interface_starti; interfacei < szInterface; ++interfacei)
	{
		subdomain_interface = sdciPtr->subdomain_interfaces[interfacei];
		SL_interfacePPtData *ptSlnPtr = subdomain_interface->timeSeqData.GetCurrentPosition();
		Compute_End_Segment_pt(subdomain_interface, interfacei, *ptSlnPtr, subdomain_interface->ts_bulkProps->bulk_leftPtr, subdomain_interface->ts_bulkProps->bulk_rightPtr, timeValue);
	}

	// printing fragment data
	if (print_fragmentation)
		Print_Interface_DSU_Fragment_OneTimeStep();

	///////////////////////////////////////////////////////////////////////////////
	// Finalizing the spacetime stats and printing
	Finalize_Subdomain_spatial_ave_sum_One_TimeStep();
}

void Subdomain_spacetime_pp_data::Set_xs()
{
	double x, xStart = sdciPtr->xm;
	xs.resize(sdciPtr->numSegments);
	x_weights.resize(sdciPtr->numSegments);

	double segLength;
	for (int segi = 0; segi < sdciPtr->numSegments; ++segi)
	{
		segLength = sdciPtr->segment_lengths[segi];
		xs[segi].resize(numSpatialPointsPerSegment);
		x_weights[segi].resize(numSpatialPointsPerSegment);
		for (int pti = 0; pti < numSpatialPointsPerSegment; ++pti)
		{
			x = xStart + segLength * spatialIntegrationPoints[pti];
			xs[segi][pti] = x;
			x_weights[segi][pti] = segLength * spatialIntegrationWeights[pti];
		}
		xStart += segLength;
	}
}

void Subdomain_spacetime_pp_data::Integrate_v4u_nonIC(int segi, int pti)
{
	OnePoint_inBulk_Fields *point_fieldsNTPtr = &allSpatialPointsRevOrder[0]->spatialPoints[segi][pti];
	OnePoint_inBulk_Fields* point_fieldsPT_or_NTPIPtr = &allSpatialPointsRevOrder[1]->spatialPoints[segi][pti];
	for (int i = 0; i < DiM; ++i)
	{
		point_fieldsNTPtr->u[i] = point_fieldsPT_or_NTPIPtr->u[i];
		for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
			point_fieldsNTPtr->u[i] += timeIntWeights_rti[rti] * allSpatialPointsRevOrder[rti]->spatialPoints[segi][pti].v[i];
	}
}

void Subdomain_spacetime_pp_data::Compute_Inside_Segment_pt(int segi, int pti, SL_Bulk_Properties & segment)
{
	OnePoint_inBulk_Fields *point_fieldsNTPtr = &allSpatialPointsRevOrder[0]->spatialPoints[segi][pti];
	OnePoint_inBulk_Fields* point_fieldsPT_or_NTPIPtr = NULL;
	bool isIC = (timeIndex == 0);
	if (!isIC)
		point_fieldsPT_or_NTPIPtr = &allSpatialPointsRevOrder[1]->spatialPoints[segi][pti];
	double x = xs[segi][pti];

	segment.Compute_Bulk_Values(x, sdciPtr->subdomain_bulk_segments[segi].interfaceLeftOfBulkPtr, sdciPtr->subdomain_bulk_segments[segi].interfaceRightOfBulkPtr,
		*point_fieldsNTPtr, isIC, point_fieldsPT_or_NTPIPtr,
			maxIter, relTol4Conv);

	// compute u from v
	if (!isIC)
		Integrate_v4u_nonIC(segi, pti);

	// update spacetime stats
	Update_Domain_space_spacetime_Integrals_fromPoint_NoInterfaceParts(segi, pti, segment);
}

subdomainInterfaceType Subdomain_spacetime_pp_data::Compute_End_Segment_pt(SL_OneInterfaceAllTimes* subdomain_interface, int interfacei, SL_interfacePPtData& ptSln, SL_Bulk_Properties* segmentPtrLeft, SL_Bulk_Properties* segmentPtrRight, double timeValue)
{
	subdomainInterfaceType sdit = sdit_interior;
	vector<int> segis, ptis;
//	bool hasLRsides[NUM_SIDES];
	vector<int> sides;
	vector<SL_Bulk_Properties*> segmentPtrs;
	bool ring_open_up_bndry = false;
	if (interfacei == 0)
	{
		sdit = sdit_left;
		sides.push_back(SDR);
//		hasLRsides[SDL] = false;
//		hasLRsides[SDR] = true;
		segis.push_back(0);
		ptis.push_back(0);
		segmentPtrs.push_back(segmentPtrRight);
	}
	else
	{
		sides.push_back(SDL);
//		hasLRsides[SDL] = true;
		segis.push_back(interfacei - 1);
		ptis.push_back(numSpatialSubsegments_BulkInterfacePoints_Print);
		segmentPtrs.push_back(segmentPtrLeft);

		if (interfacei == sdciPtr->numSegments)
		{
			if (!g_domain->isPeriodic)
			{
				sdit = sdit_right;
				//			hasLRsides[SDR] = false;
			}
			else
			{
				//			hasLRsides[SDR] = true;
				ring_open_up_bndry = true;
				sides.push_back(SDR);
				segis.push_back(0);
				ptis.push_back(0);
				segmentPtrs.push_back(segmentPtrRight);
			}
		}
		else
		{
//			hasLRsides[SDR] = true;
			sides.push_back(SDR);
			segis.push_back(interfacei);
			ptis.push_back(0);
			segmentPtrs.push_back(segmentPtrRight);
		}
	}
	unsigned side_sz = ptis.size();
	double delu0 = 0.0;
	if (side_sz == 2)
		delu0 = ptSln.sl_side_ptData[sides[1]].u_downstream_final[0] - ptSln.sl_side_ptData[sides[0]].u_downstream_final[0];
	if (ring_open_up_bndry)
		delu0 += g_domain->ring_opened1D_al * timeValue;
	for (unsigned int sdi = 0; sdi < side_sz; ++sdi)
	{
		SL_Bulk_Properties* segmentPtr = segmentPtrs[sdi];
		int segi =segis[sdi], pti = ptis[sdi];
		OnePoint_inBulk_Fields *point_fieldsNTPtr = &allSpatialPointsRevOrder[0]->spatialPoints[segi][pti];
		SL_Interface_PtData_OneSide* sl_side_ptDataPtr = &ptSln.sl_side_ptData[sides[sdi]];
		point_fieldsNTPtr->damage = subdomain_interface->getEffectiveDamage_4_InterfacialDamage_TSRs(ptSln);
		point_fieldsNTPtr->delu0 = delu0;
		CopyVec(sl_side_ptDataPtr->sigma_downstream_final, point_fieldsNTPtr->sigma);
		CopyVec(sl_side_ptDataPtr->v_downstream_final, point_fieldsNTPtr->v);
		CopyVec(sl_side_ptDataPtr->u_downstream_final, point_fieldsNTPtr->u);
		segmentPtr->Compute_Strain_from_Stress(point_fieldsNTPtr->sigma, point_fieldsNTPtr->eps);
		point_fieldsNTPtr->x = xs[segi][pti];
		point_fieldsNTPtr->time = ptSln.interface_time;

#if HAVE_SOURCE
		g_SL_desc_data.GetNonRingNonZeroTerm_SourceTerm(point_fieldsNTPtr->x, point_fieldsNTPtr->time, segmentPtr->flag, point_fieldsNTPtr->source_v, point_fieldsNTPtr->source_sigma);
#if RING_PROBLEM
		point_fieldsNTPtr->v_r = ptSln.v_r_final;
		double sigma_theta_source_final = segmentPtr->E_iso * point_fieldsNTPtr->v_r / g_domain->ring_R; // E v_r / R (same on both sides)
		point_fieldsNTPtr->source_sigma[0] += sigma_theta_source_final;
#endif
#endif
		// updating spacetime integrals
		Update_Domain_space_spacetime_Integrals_fromPoint_NoInterfaceParts(segi, pti, *segmentPtr);
	}

	// updating spacetime interface statistics
	Subdomain_spatial_ave_sum* spatial_ave_sum = spatial_ave_sums[0];
	double D = subdomain_interface->getEffectiveDamage_4_InterfacialDamage_TSRs(ptSln);
	if (side_sz == 2)
	{
		spatial_ave_sum->mean_interface_damage += D;
		spatial_ave_sum->sdiv_interface_damage += D * D;
		spatial_ave_sum->max_interface_damage = MAX(spatial_ave_sum->max_interface_damage, D);
		spatial_ave_sum->min_interface_damage = MIN(spatial_ave_sum->min_interface_damage, D);
		spatial_ave_sum->max_interface_damage_source = MAX(spatial_ave_sum->max_interface_damage_source, ptSln.interface_damage_source_final);
		spatial_ave_sum->max2EffDelU = MAX(spatial_ave_sum->max2EffDelU, ptSln.maxEffDelU);
	}
	if (sdit == sdit_interior)
	{
		VEC delv, delu, sigma;
		CopyVec(ptSln.sl_side_ptData[SDL].sigma_downstream_final, sigma);
		SubtractVec(ptSln.sl_side_ptData[SDR].v_downstream_final, ptSln.sl_side_ptData[SDL].v_downstream_final, delv);
		SubtractVec(ptSln.sl_side_ptData[SDR].u_downstream_final, ptSln.sl_side_ptData[SDL].u_downstream_final, delu);
		if (ring_open_up_bndry)
		{
			delu[0] += g_domain->ring_opened1D_al * timeValue;
			delv[0] += g_domain->ring_opened1D_al;
		}
		// computing powerIDissipation
		for (int i = 0; i < DiM; ++i)
			spatial_ave_sum->powerIDiss_vec[i] += delv[i] * sigma[i];

		CopyVec(sigma, spatial_ave_sum->interface_sigma_delu_vec[interfacei].first);
		CopyVec(delu, spatial_ave_sum->interface_sigma_delu_vec[interfacei].second);

		for (int i = 0; i < DiM; ++i)
		{
			spatial_ave_sum->eps_ave_bulk_intfc[i] += delu[i]; // division by length is done later
			spatial_ave_sum->energyIDiss_Recoverable += delu[i] * sigma[i]; // multiplication by 0.5 is done later
		}
	}
	else if (sdit == sdit_left)
	{
		hasLeft = true;
		CopyVec(ptSln.sl_side_ptData[SDR].u_downstream_final, spatial_ave_sum->u_L);
		spatial_ave_sum->power_L = 0.0;
		for (int i = 0; i < DiM; ++i)
		{
			spatial_ave_sum->sigman_L[i] = ptSln.sl_side_ptData[SDR].sigma_downstream_final[i];
			spatial_ave_sum->power_L -= ptSln.sl_side_ptData[SDR].v_downstream_final[i] * spatial_ave_sum->sigman_L[i];
		}
	}
	else
	{
		hasRight = true;
		CopyVec(ptSln.sl_side_ptData[SDL].u_downstream_final, spatial_ave_sum->u_R);
		CopyVec(ptSln.sl_side_ptData[SDL].sigma_downstream_final, spatial_ave_sum->sigman_R);
		spatial_ave_sum->power_R = 0.0;
		for (int i = 0; i < DiM; ++i)
			spatial_ave_sum->power_R += ptSln.sl_side_ptData[SDL].v_downstream_final[i] * spatial_ave_sum->sigman_R[i];
	}
	return sdit;
}

void Subdomain_spacetime_pp_data::Finalize_Subdomain_spatial_ave_sum_One_TimeStep()
{
	Subdomain_spatial_ave_sum* spatial_ave_sum = spatial_ave_sums[0];
	spatial_ave_sum->timeIndex = timeIndex;
	spatial_ave_sum->timeVal = times[timeIndex];

	// averages
	FactorVec(spatial_ave_sum->u_ave, sdciPtr->inv_length);
	FactorVec(spatial_ave_sum->v_ave, sdciPtr->inv_length);
	FactorVec(spatial_ave_sum->eps_ave_bulk_blk , sdciPtr->inv_length);
	FactorVec(spatial_ave_sum->eps_ave_bulk_intfc, sdciPtr->inv_length);
	AddVec(spatial_ave_sum->eps_ave_bulk_blk, spatial_ave_sum->eps_ave_bulk_intfc, spatial_ave_sum->eps_ave_bulk);

	FactorVec(spatial_ave_sum->sigma_ave_bulk, sdciPtr->inv_length);

	if (!g_domain->isPeriodic)
	{
		for (int i = 0; i < DiM; ++i)
		{
			spatial_ave_sum->sigma_ave_bc[i] = (spatial_ave_sum->sigman_R[i] * sdciPtr->xM - spatial_ave_sum->sigman_L[i] * sdciPtr->xm) * sdciPtr->inv_length;
			spatial_ave_sum->eps_ave_bc[i] = (spatial_ave_sum->u_R[i] - spatial_ave_sum->u_L[i]) * sdciPtr->inv_length;
		}
		SubtractVec(spatial_ave_sum->sigma_ave_bc, spatial_ave_sum->sigma_ave_bulk, spatial_ave_sum->sigma_ave_bc_minus_bulk);
		SubtractVec(spatial_ave_sum->eps_ave_bc, spatial_ave_sum->eps_ave_bulk, spatial_ave_sum->eps_ave_bc_minus_bulk);
	}
	else
	{
		CopyVec(spatial_ave_sum->eps_ave_bulk, spatial_ave_sum->eps_ave_bc);
		CopyVec(spatial_ave_sum->sigma_ave_bulk, spatial_ave_sum->sigma_ave_bc);
		setValue(spatial_ave_sum->eps_ave_bc_minus_bulk, 0.0);
		setValue(spatial_ave_sum->sigma_ave_bc_minus_bulk, 0.0);
	}
	spatial_ave_sum->energy_eps_recoverable_bc = 0.0;
	spatial_ave_sum->energy_eps_recoverable_bulk = 0.0;
	for (int i = 0; i < DiM; ++i)
	{
		spatial_ave_sum->energy_eps_recoverable_bc += spatial_ave_sum->eps_ave_bc[i] * spatial_ave_sum->sigma_ave_bc[i];
		spatial_ave_sum->energy_eps_recoverable_bulk += spatial_ave_sum->eps_ave_bulk[i] * spatial_ave_sum->sigma_ave_bulk[i];
	}
	spatial_ave_sum->energy_eps_recoverable_bc *= 0.5;
	spatial_ave_sum->energy_eps_recoverable_bulk *= 0.5;

	spatial_ave_sum->phi_dsum = spatial_ave_sum->K_dsum + spatial_ave_sum->U_dsum;

#if HAVE_SOURCE
#if RING_PROBLEM
	spatial_ave_sum->v_r_ave *= sdciPtr->inv_length;
#endif
	spatial_ave_sum->source_e_dsum = spatial_ave_sum->source_e_v_dsum + spatial_ave_sum->source_e_sigma_dsum;
#if HAVE_SOURCE_ORDER0_q
	spatial_ave_sum->damping_e_dsum = spatial_ave_sum->damping_e_v_dsum + spatial_ave_sum->damping_e_sigma_dsum;
#endif
#endif

	// time integrations:
	bool isIC = (timeIndex == 0);
	Subdomain_spatial_ave_sum* spatial_ave_sumPrev = NULL;
	if (!isIC)
	{
		spatial_ave_sumPrev = spatial_ave_sums[1];
#if HAVE_SOURCE
		spatial_ave_sum->source_e_v_dtsum = spatial_ave_sumPrev->source_e_v_dtsum;
#if HAVE_SOURCE_ORDER0_q
		spatial_ave_sum->damping_e_v_dtsum = spatial_ave_sumPrev->damping_e_v_dtsum;
#endif
#endif
		spatial_ave_sum->energy_eps_total_bc = spatial_ave_sumPrev->energy_eps_total_bc;
		spatial_ave_sum->energy_eps_total_bulk = spatial_ave_sumPrev->energy_eps_total_bulk;
		spatial_ave_sum->powerIDiss = 0.0;
		spatial_ave_sum->eneIDiss = 0.0;
		for (int i = 0; i < DiM; ++i)
		{
			spatial_ave_sum->energy_eps_total_bc += 0.5 * (spatial_ave_sum->eps_ave_bc[i] - spatial_ave_sumPrev->eps_ave_bc[i]) *
				(spatial_ave_sum->sigma_ave_bc[i] + spatial_ave_sumPrev->sigma_ave_bc[i]);
			spatial_ave_sum->energy_eps_total_bulk += 0.5 * (spatial_ave_sum->eps_ave_bulk[i] - spatial_ave_sumPrev->eps_ave_bulk[i]) *
				(spatial_ave_sum->sigma_ave_bulk[i] + spatial_ave_sumPrev->sigma_ave_bulk[i]);

			if (!g_domain->isPeriodic)
			{
				spatial_ave_sum->impulse_L[i] = spatial_ave_sumPrev->impulse_L[i];
				for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
					spatial_ave_sum->impulse_L[i] -= timeIntWeights_rti[rti] * spatial_ave_sums[rti]->sigman_L[i];

				spatial_ave_sum->impulse_R[i] = spatial_ave_sumPrev->impulse_R[i];
				for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
					spatial_ave_sum->impulse_R[i] += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->sigman_R[i];
			}
			else
			{
				setValue(spatial_ave_sum->impulse_L, 0.0);
				setValue(spatial_ave_sum->impulse_R, 0.0);
			}
			spatial_ave_sum->eneIDiss_vec[i] = spatial_ave_sumPrev->eneIDiss_vec[i];
			if (timeIndex >= 1) // (timeIndex > 1)
			{
				unsigned int numInterface = spatial_ave_sum->interface_sigma_delu_vec.size();
				for (unsigned int interfacei = 0; interfacei < numInterface; ++interfacei)
				{
					VEC *sig_prev = &spatial_ave_sumPrev->interface_sigma_delu_vec[interfacei].first;
					VEC *delu_prev = &spatial_ave_sumPrev->interface_sigma_delu_vec[interfacei].second;
					VEC *sig_new = &spatial_ave_sum->interface_sigma_delu_vec[interfacei].first;
					VEC *delu_new = &spatial_ave_sum->interface_sigma_delu_vec[interfacei].second;
					spatial_ave_sum->eneIDiss_vec[i] += 0.5 * ((*sig_prev)[i] + (*sig_new)[i]) * ((*delu_new)[i] - (*delu_prev)[i]);
				}
#if HAVE_SOURCE
				for (int segi = 0; segi < sdciPtr->numSegments; ++segi)
				{
					SL_Bulk_Properties *segmentPtr = sdciPtr->subdomain_bulk_segments[segi].bulkPtr;
					double rho = segmentPtr->rho;
					for (int pti = 0; pti < numSpatialPointsPerSegment; ++pti)
					{
						double spatialWeight = x_weights[segi][pti], factor = 0.5 * rho * spatialWeight;
						OnePoint_inBulk_Fields *point_fieldsNTPtr = &allSpatialPointsRevOrder[0]->spatialPoints[segi][pti];
						OnePoint_inBulk_Fields *point_fieldsPT_or_NTPIPtr = &allSpatialPointsRevOrder[1]->spatialPoints[segi][pti];

						double delu, tmp_e_lin_m = 0.0, tmp_e_damping_v = 0.0, damping_source_v;
						for (int i = 0; i < DiM; ++i)
						{
							delu = point_fieldsNTPtr->u[i] - point_fieldsPT_or_NTPIPtr->u[i];
							tmp_e_lin_m += (point_fieldsNTPtr->source_v[i] + point_fieldsPT_or_NTPIPtr->source_v[i]) * delu;

#if HAVE_SOURCE_ORDER0_q
							damping_source_v =
								(segmentPtr->D_vv * point_fieldsNTPtr->v[i] + segmentPtr->D_vsigma * point_fieldsNTPtr->sigma[i]) +
								(segmentPtr->D_vv * point_fieldsPT_or_NTPIPtr->v[i] + segmentPtr->D_vsigma * point_fieldsPT_or_NTPIPtr->sigma[i]);

							tmp_e_damping_v += damping_source_v * delu;
#endif
						}
						spatial_ave_sum->source_e_v_dtsum += factor * tmp_e_lin_m;
#if HAVE_SOURCE_ORDER0_q
						spatial_ave_sum->damping_e_v_dtsum += factor * tmp_e_damping_v;
#endif
					}
				}
#endif
			}
			else
			{
				THROW("This is just a check. We shouldn't be getting here\n");
				for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
					spatial_ave_sum->eneIDiss_vec[i] += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->powerIDiss_vec[i];
#if HAVE_SOURCE
				for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
					spatial_ave_sum->source_e_v_dtsum += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->source_e_v_dsum;
#if HAVE_SOURCE_ORDER0_q
				for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
					spatial_ave_sum->damping_e_v_dtsum += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->damping_e_v_dsum;
#endif
#endif
			}
			spatial_ave_sum->powerIDiss += spatial_ave_sum->powerIDiss_vec[i];
			spatial_ave_sum->eneIDiss += spatial_ave_sum->eneIDiss_vec[i];

#if HAVE_SOURCE
			spatial_ave_sum->source_linMomentum_dtsum[i] = spatial_ave_sumPrev->source_linMomentum_dtsum[i];
			for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
				spatial_ave_sum->source_linMomentum_dtsum[i] += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->source_linMomentum_dsum[i];

#if HAVE_SOURCE_ORDER0_q
			spatial_ave_sum->damping_linMomentum_dtsum[i] = spatial_ave_sumPrev->damping_linMomentum_dtsum[i];
			for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
				spatial_ave_sum->damping_linMomentum_dtsum[i] += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->damping_linMomentum_dsum[i];
#endif
#endif
		}
		if (!g_domain->isPeriodic)
		{
			AddVec(spatial_ave_sum->impulse_L, spatial_ave_sum->impulse_R, spatial_ave_sum->impulse_BC);
			spatial_ave_sum->energy_L = spatial_ave_sumPrev->energy_L;
			spatial_ave_sum->energy_R = spatial_ave_sumPrev->energy_R;
#if 0
			for (unsigned int i = 0; i < DiM; ++i)
			{
				spatial_ave_sum->energy_L -= 0.5 * (spatial_ave_sumPrev->sigman_L[i] + spatial_ave_sum->sigman_L[i]) *
					(spatial_ave_sum->u_L[i] - spatial_ave_sumPrev->u_L[i]);
				spatial_ave_sum->energy_R += 0.5 * (spatial_ave_sumPrev->sigman_R[i] + spatial_ave_sum->sigman_R[i]) *
					(spatial_ave_sum->u_R[i] - spatial_ave_sumPrev->u_R[i]);
			}
#else
			for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
				spatial_ave_sum->energy_L += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->power_L;
			for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
				spatial_ave_sum->energy_R += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->power_R;
#endif
			spatial_ave_sum->energy_BC = spatial_ave_sum->energy_L + spatial_ave_sum->energy_R;
		}
		else
		{
			setValue(spatial_ave_sum->impulse_BC, 0.0);
			spatial_ave_sum->energy_L = 0.0;
			spatial_ave_sum->energy_R = 0.0;
			spatial_ave_sum->energy_BC = 0.0;
		}
#if HAVE_SOURCE
		spatial_ave_sum->source_e_sigma_dtsum = spatial_ave_sumPrev->source_e_sigma_dtsum;
		for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
			spatial_ave_sum->source_e_sigma_dtsum += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->source_e_sigma_dsum;
		spatial_ave_sum->source_e_dtsum = spatial_ave_sum->source_e_v_dtsum + spatial_ave_sum->source_e_sigma_dtsum;
#if HAVE_SOURCE_ORDER0_q
		spatial_ave_sum->damping_e_sigma_dtsum = spatial_ave_sumPrev->damping_e_sigma_dtsum;
		for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
			spatial_ave_sum->damping_e_sigma_dtsum += timeIntWeights_rti[rti] * spatial_ave_sums[rti]->damping_e_sigma_dsum;
		spatial_ave_sum->damping_e_dtsum = spatial_ave_sum->damping_e_v_dtsum + spatial_ave_sum->damping_e_sigma_dtsum;
#endif
#endif
		spatial_ave_sum->energy_eps_diss_bc = spatial_ave_sum->energy_eps_total_bc - spatial_ave_sum->energy_eps_recoverable_bc;
		spatial_ave_sum->energy_eps_diss_bulk = spatial_ave_sum->energy_eps_total_bulk - spatial_ave_sum->energy_eps_recoverable_bulk;

		if (g_domain->b_ring_opened1D_kinetic_energy_on_full_vTheta)
		{
			spatial_ave_sum->source_e_dtsum = spatial_ave_sumPrev->source_e_dtsum;
			double aL = g_domain->ring_opened1D_al;
			for (unsigned int rti = 0; rti < timeIntWeights_rti_size; ++rti)
				spatial_ave_sum->source_e_dtsum += aL * timeIntWeights_rti[rti] * spatial_ave_sums[rti]->sigma_ave_bulk[0];
		}
	}
	else
	{
		spatial_ave_sum->energy_eps_total_bc = spatial_ave_sum->energy_eps_recoverable_bc;
		spatial_ave_sum->energy_eps_total_bulk = spatial_ave_sum->energy_eps_recoverable_bulk;
		spatial_ave_sum->energy_eps_diss_bc = 0.0;
		spatial_ave_sum->energy_eps_diss_bulk = 0.0;
	}
	// calculating numerical errors
	if (timeIndex == 0)
		phi0 = spatial_ave_sum->phi_dsum;
	spatial_ave_sum->input_energy = phi0 + spatial_ave_sum->energy_BC;
	spatial_ave_sum->phys_diss_tot = spatial_ave_sum->eneIDiss;
	spatial_ave_sum->input_energy += spatial_ave_sum->source_e_dtsum;
#if HAVE_SOURCE
//	spatial_ave_sum->input_energy += spatial_ave_sum->source_e_dtsum;
#if HAVE_SOURCE_ORDER0_q
	spatial_ave_sum->phys_diss_tot += spatial_ave_sum->damping_e_dtsum;
#endif
#endif
	spatial_ave_sum->numerial_energy_diss = spatial_ave_sum->input_energy - spatial_ave_sum->phi_dsum - spatial_ave_sum->phys_diss_tot;

	if (timeIndex == 0)
	{
		linMomentumZero = 0.0;
		for (int i = 0; i < DiM; ++i)
			linMomentumZero = spatial_ave_sum->linMomentum_dsum[i];
	}
	for (int i = 0; i < DiM; ++i)
		spatial_ave_sum->numerical_linmom_error[i] = linMomentumZero - spatial_ave_sum->linMomentum_dsum[i] + spatial_ave_sum->impulse_BC[i];
#if HAVE_SOURCE
	for (int i = 0; i < DiM; ++i)
		spatial_ave_sum->numerical_linmom_error[i] += spatial_ave_sum->source_linMomentum_dtsum[i];
#if HAVE_SOURCE_ORDER0_q
	for (int i = 0; i < DiM; ++i)
		spatial_ave_sum->numerical_linmom_error[i] -= spatial_ave_sum->damping_linMomentum_dtsum[i];
#endif
#endif

	if (sdciPtr->numSegments > 0)
	{
		double fact = 1.0 / (sdciPtr->numSegments - 1);
		spatial_ave_sum->mean_interface_damage *= fact;
		spatial_ave_sum->sdiv_interface_damage = (spatial_ave_sum->sdiv_interface_damage - spatial_ave_sum->mean_interface_damage * spatial_ave_sum->mean_interface_damage) * fact;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////// stuff that can be done at the end of time step
	spatial_ave_sum->energyIDiss_Recoverable *= 0.5;
	spatial_ave_sum->phys_diss_lost = spatial_ave_sum->phys_diss_tot - spatial_ave_sum->energyIDiss_Recoverable;
	spatial_ave_sum->phys_diss_interface_lost = spatial_ave_sum->eneIDiss - spatial_ave_sum->energyIDiss_Recoverable;

	/////////////////////////////////////////////////////////////////////////////////////////////
	// printing the stat:
	if (timeIndex == 0)
		Subdomain_spatial_ave_sum::Subdomain_spatial_ave_sum_Data_Write_Header(*out_sd_summary);
	spatial_ave_sum->Subdomain_spatial_ave_sum_Write_Data(*out_sd_summary, sdciPtr->inv_length);

	bool b_print = false;
	if (print_space_points)
		b_print = (b_print || (timeIndex % numTimeStep_BulkInterfacePoints_Print == 0) || (timeIndex == numTimes));
	if (b_print)
	{
		string fileName;
		string specificName = "BulkInterface";
		GetSubdomainIndexed_TimeIndexed_FileName(fileName, sdciPtr->subdomain_number, timeIndex, specificName);
		fstream out(fileName.c_str(), ios::out);
		if (!out.is_open())
		{
			cout << "fileName\t" << fileName << '\n';
			THROW("Cannot open file\n");
		}
		allSpatialPointsRevOrder[0]->Subdomain_oneTime_spatial_points_Print_Header_Data(out);
	}
}

void Subdomain_spacetime_pp_data::Print_Interface_DSU_Fragment_OneTimeStep()
{
	if ((timeIndex % numTimeStep_Interface_DSU_Fragment_Print != 0) && (timeIndex != numTimes))
		return;
	SL_OneInterfaceAllTimes* subdomain_interface;
	unsigned szInterface = sdciPtr->subdomain_interfaces.size();
	string time_str;

	string fileName;
	string specificName = "Interface_DSU_Fragment";
	GetSubdomainIndexed_TimeIndexed_FileName(fileName, sdciPtr->subdomain_number, timeIndex, specificName);
	fstream outfrag(fileName.c_str(), ios::out);
	outfrag << "Damage" << "\t" << "maxEffDelU";
//	for (unsigned int di = 0; di < DiM; ++di)
//		outfrag << "\tdelu" << di;
	for (unsigned int di = 0; di < DiM; ++di)
		outfrag << "\tuL" << di;
	for (unsigned int di = 0; di < DiM; ++di)
		outfrag << "\tuR" << di;
	for (unsigned int di = 0; di < DiM; ++di)
		outfrag << "\tsigma" << di;
	outfrag << '\n';
	unsigned int lastInterfaceNo = szInterface - 1;
	unsigned int interface_starti = 0;
	if (g_domain->isPeriodic)
		interface_starti = 1;
	for (unsigned int interfacei = interface_starti; interfacei < szInterface; ++interfacei)
	{
		double x = sdciPtr->subdomain_interface_xs[interfacei];
		subdomain_interface = sdciPtr->subdomain_interfaces[interfacei];
		SL_interfacePPtData *ptSlnPtr = subdomain_interface->timeSeqData.GetCurrentPosition();
		outfrag << ptSlnPtr->interface_damage_final << '\t';
		outfrag << ptSlnPtr->maxEffDelU << '\t';
		if (!g_domain->b_ring_opened1D)
		{
			for (unsigned int di = 0; di < DiM; ++di)
			{
				//			double delu = ptSlnPtr->sl_side_ptData[SDR].u_downstream_final[di] - ptSlnPtr->sl_side_ptData[SDL].u_downstream_final[di];
				//			outfrag << delu << '\t';
				outfrag << ptSlnPtr->sl_side_ptData[SDL].u_downstream_final[di] << '\t';
			}
			for (unsigned int di = 0; di < DiM; ++di)
				outfrag << ptSlnPtr->sl_side_ptData[SDR].u_downstream_final[di] << '\t';
		}
		else
		{
			bool lastInterface = (interfacei == lastInterfaceNo);
			unsigned int di = 0;
			double uL = ptSlnPtr->sl_side_ptData[SDL].u_downstream_final[di];
			double uR = ptSlnPtr->sl_side_ptData[SDR].u_downstream_final[di];
			double ax = g_SL_desc_data.load_parameters[0] * x, axt = ax * times[timeIndex];
			uL -= axt;
			if (!lastInterface)
				uR -= axt;
			else
				uR += axt;
			outfrag << uL << '\t' << uR << '\t';
		}
		if (interfacei > 0)
		{
			for (unsigned int di = 0; di < DiM; ++di)
				outfrag << ptSlnPtr->sl_side_ptData[SDL].sigma_downstream_final[di] << '\t';
		}
		else
		{
			for (unsigned int di = 0; di < DiM; ++di)
				outfrag << ptSlnPtr->sl_side_ptData[SDR].sigma_downstream_final[di] << '\t';
		}
		outfrag << '\n';
	}
}


void Subdomain_spacetime_pp_data::Update_Domain_space_spacetime_Integrals_fromPoint_NoInterfaceParts(int segi, int pti, SL_Bulk_Properties & segment)
{
	Subdomain_spatial_ave_sum* spatial_ave_sum = spatial_ave_sums[0];
	OnePoint_inBulk_Fields *point_fieldsNTPtr = &allSpatialPointsRevOrder[0]->spatialPoints[segi][pti];
	double spatialWeight = x_weights[segi][pti];
	double K = 0.0, U = 0.0;// , phi;
	VEC linMomentum;
	VEC vel, u;
	double ax = 0.0;
	CopyVec(point_fieldsNTPtr->v, vel);
	CopyVec(point_fieldsNTPtr->u, u);
	if (g_domain->b_ring_opened1D_kinetic_energy_on_full_vTheta)
	{
		ax = g_SL_desc_data.load_parameters[0] * point_fieldsNTPtr->x;
		vel[0] -= ax;
		double time = times[timeIndex];
		u[0] -= ax * time;
	}
	// linear momentum and energy:
	for (int i = 0; i < DiM; ++i)
	{
		linMomentum[i] = segment.rho * vel[i];
		K += vel[i] * vel[i];
		U += point_fieldsNTPtr->sigma[i] * point_fieldsNTPtr->eps[i];
	}
	K *= (0.5 * segment.rho);
	K += g_domain->ring_opened1D_kinetic_energy_vr;
	U *= 0.5;
//	phi = K + U;

	//////////////////////////////////////

	// updating K, U, phi of domain
	spatial_ave_sum->K_dsum += K * spatialWeight;
	spatial_ave_sum->U_dsum += U * spatialWeight;
	// later simply added
//	spatial_ave_sum->phi_dsum += phi * spatialWeight;

	VEC *sigmaPtr = &point_fieldsNTPtr->sigma, *epsPtr = &point_fieldsNTPtr->eps;
	for (int i = 0; i < DiM; ++i)
	{
		spatial_ave_sum->linMomentum_dsum[i] += linMomentum[i] * spatialWeight;
		spatial_ave_sum->eps_ave_bulk_blk[i] += (*epsPtr)[i] * spatialWeight;
		spatial_ave_sum->sigma_ave_bulk[i] += (*sigmaPtr)[i] * spatialWeight;
		spatial_ave_sum->u_ave[i] += u[i] * spatialWeight;
		spatial_ave_sum->v_ave[i] += vel[i] * spatialWeight;
	}

	// spacetime quantities are added in space and delt part is incorporated later
#if HAVE_SOURCE
#if RING_PROBLEM
	spatial_ave_sum->v_r_ave += point_fieldsNTPtr->v_r * spatialWeight;
#endif
	double source_e_v = 0.0;
	double source_e_sigma = 0.0;
	VEC source_linMomentum;

	for (int i = 0; i < DiM; ++i)
	{
		source_linMomentum[i] = segment.rho * point_fieldsNTPtr->source_v[i];
		source_e_v += vel[i] * source_linMomentum[i];
		source_e_sigma += point_fieldsNTPtr->eps[i] * point_fieldsNTPtr->source_sigma[i];
	}

	for (int i = 0; i < DiM; ++i)
		spatial_ave_sum->source_linMomentum_dsum[i] += source_linMomentum[i] * spatialWeight;

	spatial_ave_sum->source_e_v_dsum += source_e_v * spatialWeight;
	spatial_ave_sum->source_e_sigma_dsum += source_e_sigma * spatialWeight;
	//	source_e = source_e_v + source_e_sigma;
	// spatial_ave_sum->source_e_dsum += source_e * spatialWeight;

	// damping contribution to linear momentum and energy
#if HAVE_SOURCE_ORDER0_q
	/// damping terms
	VEC damping_v, damping_sigma;
	setValue(damping_v, 0.0);
	setValue(damping_sigma, 0.0);
	for (int i = 0; i < DiM; ++i)
	{
		double v = point_fieldsNTPtr->v[i];
		if (g_domain->b_ring_opened1D_damping_on_full_vTheta)
			v = vel[i];
		damping_v[i] += (segment.D_vv * v + segment.D_vsigma * point_fieldsNTPtr->sigma[i]);
		damping_sigma[i] += (segment.D_sigmav * v + segment.D_sigmasigma * point_fieldsNTPtr->sigma[i]);
	}

	double damping_e_v = 0.0;
	double damping_e_sigma = 0.0;
	VEC damping_linMomentum;

	for (int i = 0; i < DiM; ++i)
	{
		double v = point_fieldsNTPtr->v[i];
		if (g_domain->b_ring_opened1D_damping_on_full_vTheta)
			v = vel[i];
		damping_linMomentum[i] = segment.rho * damping_v[i];
		damping_e_v += v * damping_linMomentum[i];
		damping_e_sigma += point_fieldsNTPtr->eps[i] * damping_sigma[i];
	}

	for (int i = 0; i < DiM; ++i)
		spatial_ave_sum->damping_linMomentum_dsum[i] += spatialWeight * damping_linMomentum[i];
	spatial_ave_sum->damping_e_v_dsum += damping_e_v * spatialWeight;
	spatial_ave_sum->damping_e_sigma_dsum += damping_e_sigma * spatialWeight;
#endif
#endif
}

