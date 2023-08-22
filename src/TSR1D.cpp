#include "TSR1D.h"
#include "commonMacros.h"
#include "globalFunctions.h"
#include "commonFunctions.h"

/////////////////
string getName(TSR_loadingStages dat)
{
	if (dat == TSR_loadingStages_none)
		return "none";
	if (dat == tsrl_dam0_sigI_neg)
		return "dam0_sigI_neg";
	if (dat == tsrl_dam0_sigI_less_sigC)
		return "dam0_sigI_less_sigC";
	if (dat == tsrl_dam0_sigI_over_sigC)
		return "dam0_sigI_over_sigC";
	if (dat == tsrl_damMid_loadingBranch)
		return "damMid_loadingBranch";
	if (dat == tsrl_damMid_unloadingBranchPosDelv)
		return "damMid_unloadingBranchPosDelv";
	if (dat == tsrl_damMid_unloadingBranchNegDelv)
		return "damMid_unloadingBranchNegDelv";
	if (dat == tsrl_damMid_delun_neg_sigI_neg)
		return "damMid_delun_neg_sigI_neg";
	if (dat == tsrl_damMid_delun_neg_sigI_pos)
		return "damMid_delun_neg_sigI_pos";
	if (dat == tsr1_dam1_delun_pos)
		return "dam1_delun_pos";
	if (dat == tsr1_dam1_delun_neg)
		return "dam1_delun_neg";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, TSR_loadingStages& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= TSR_loadingStages_SIZE)
			THROW("too large of a number\n");
		typeVal = (TSR_loadingStages)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < TSR_loadingStages_SIZE; ++i)
	{
		typeVal = (TSR_loadingStages)i; // casting integer to TSR_loadingStages, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading TSR_loadingStages\n");
}

//operator for output
ostream& operator<<(ostream& out, TSR_loadingStages dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, TSR_loadingStages& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

TSR1d_fields_1side::TSR1d_fields_1side()
{
	s = 0.0;
	u = 0.0;
	v = 0.0;
}

TSR1d_2Sides::TSR1d_2Sides()
{
	delu = 0.0;
	delv = 0.0;
}

TSR1d_stepn_to_np1::TSR1d_stepn_to_np1()
{
	tol_TSR = g_tol_TSR;
	Zeff = -1;
	step_n_max_delu = 0.0;
	iterN = 0;
	max_iterN = 4;
}

void TSR1d_stepn_to_np1::Zero_TSR1d_stepn_to_np1()
{
	step_n_tsrStage = TSR_loadingStages_none;
	validShortStage = TSR_OrtizStages_slnFound;

	step_np1_tsrStage = TSR_loadingStages_none;
	for (int i = 0; i < TSR_OrtizStages_SIZE; ++i)
		calcType[i] = -1;
	message = "";
}

void TSR1d_stepn_to_np1::Initialize_TSR1d_stepn_to_np1(double tol_TSRFactor)
{
	Zero_TSR1d_stepn_to_np1();

	double Zl = Zs[SDL], Zr = Zs[SDR];
	double ZsInv = 1.0 / (Zl + Zr);
	Zeff = 2.0 * Zl * Zr * ZsInv;
	step_n.delu = step_n.sides[SDR].u - step_n.sides[SDL].u;
	step_n.delv = step_n.sides[SDR].v - step_n.sides[SDL].v;
	sI = (Zl * ws[SDR] + Zr * ws[SDL]) * ZsInv;
	tol_TSR = g_tol_TSR * tol_TSRFactor;
	toldc = tol_TSR * deltaC;
	tolsc = tol_TSR * sigmaC;
	// these short terms can change later, but are good starting points
	maxDelUEqZero = step_n_max_delu < 1e-2 * toldc;
	if (maxDelUEqZero)
	{
		step_n_max_sigma = sigmaC;
		validShortStage = tsro_A;
	}
	else if (step_n_max_delu >= deltaC)
	{ 
		step_n_max_sigma = 0.0;
		validShortStage = tsro_D;
	}
	else
	{
		step_n_max_sigma = sigmaC * (1 - step_n_max_delu / deltaC);
		validShortStage = tsro_B;
	}
}

bool TSR1d_stepn_to_np1::Compute_Stage_Aux()
{
	if (validShortStage == tsro_D) // max del u > deltaC -> damage = 1
	{
		// test if we stay on D (delu > 0, separation) or E (delu = 0 - contact)
		TSR_OrtizStages stageShort2TestNow = Test_Stage_D();
		if (stageShort2TestNow == TSR_OrtizStages_invalid)
			return false;
		if (stageShort2TestNow == tsro_E)
		{
			TSR_OrtizStages stageShort2TestNow = Test_Stages_DelU_eqZero_A_or_E(tsro_E);
			if (stageShort2TestNow == TSR_OrtizStages_invalid)
				return false;
			if (stageShort2TestNow != TSR_OrtizStages_slnFound)
			{
				if (iterN == max_iterN)
					message = "started from branch D(full damage, delu > 0), shifted to branch E (contact) but this case failed. So, there is no valid solution";
				return false;
			}
		}
		else if (stageShort2TestNow != TSR_OrtizStages_slnFound)
		{
			if (iterN == max_iterN)
				message = "At full damage only separation (case D, delu > 0) or contact (case E, delu = 0) are valid.";
			return false;
		}
		return true;
	}
	else if (validShortStage == tsro_A) // max delu = 0 -> damage = 0
	{
		// can end up in case A (0 < sigma < sigmaC), case E (sigma < 0) or case B (loading branch)
		TSR_OrtizStages stageShort2TestNow = Test_Stages_DelU_eqZero_A_or_E(tsro_A);
		if (stageShort2TestNow == TSR_OrtizStages_invalid)
			return false;
		if (stageShort2TestNow == tsro_B)
		{
			stageShort2TestNow = Test_Stage_B();
			if (stageShort2TestNow == TSR_OrtizStages_invalid)
				return false;
			if (stageShort2TestNow != TSR_OrtizStages_slnFound)
			{
				if (iterN == max_iterN)
					message = "started from branch A(zero damage), shifted to branch B (loading) but this case failed. So, there is no valid solution";
				return false;
			}
		}
		else if (stageShort2TestNow != TSR_OrtizStages_slnFound)
		{
			if (iterN == max_iterN)
				message = "From branch A (zero damage), can only stay there (delu = 0, branches A,E) or get to loading branch B. This option is not valid";
			return false;
		}
		return true;
	}
	else // intermediate damage
	{
		// can end up in loading (caseB), unloading (caseC), full damage (caseD), or contact (caseE) modes now
		// this will start as one of the cases B, C, or E
		if (step_n.delu < 1e-2 * toldc) // test if contact is the correct solution
		{
			TSR_OrtizStages stageShort2TestNow = Test_Stages_DelU_eqZero_A_or_E(tsro_E);
			if (stageShort2TestNow == TSR_OrtizStages_invalid)
				return false;
			// 3 levels are added here as it may jump from E to C to B!
			if (stageShort2TestNow != TSR_OrtizStages_slnFound) // means case E is not valid, other cases to be tested
			{
				stageShort2TestNow = Test_newcase(stageShort2TestNow);
				if (stageShort2TestNow != TSR_OrtizStages_slnFound) // means case E is not valid, other cases to be tested
				{
					stageShort2TestNow = Test_newcase(stageShort2TestNow);
					if (stageShort2TestNow == TSR_OrtizStages_invalid)
						return false;
					if (stageShort2TestNow != TSR_OrtizStages_slnFound)
					{
						if (iterN == max_iterN)
							message = "started from branch E (contact mode), didn't stay there and the new suggested branch didn't work either";
						return false;
					}
				}
			}
			return true;
		}
		if (step_n_max_delu - step_n.delu <= 0) // toldc) // start with loading branch and see if it stays in loading branch
		{
			TSR_OrtizStages stageShort2TestNow = Test_Stage_B();
			if (stageShort2TestNow == TSR_OrtizStages_invalid)
				return false;
			if (stageShort2TestNow != TSR_OrtizStages_slnFound) // means case B is not valid, other cases to be tested
			{
				stageShort2TestNow = Test_newcase(stageShort2TestNow);
				if (stageShort2TestNow == TSR_OrtizStages_invalid)
					return false;
				if (stageShort2TestNow != TSR_OrtizStages_slnFound)
				{
					stageShort2TestNow = Test_newcase(stageShort2TestNow);
					if (stageShort2TestNow == TSR_OrtizStages_invalid)
						return false;
					if (stageShort2TestNow != TSR_OrtizStages_slnFound)
					{
						if (iterN == max_iterN)
							message = "started from branch B (contact mode), didn't stay there and the new suggested branch didn't work either\n";
						return false;
					}
				}
			}
			return true;
		}
		// the only valid branch is C (A and E are impossible to start from)
		TSR_OrtizStages stageShort2TestNow = Test_Stage_C();
		if (stageShort2TestNow == TSR_OrtizStages_invalid)
			return false;
		if (stageShort2TestNow != TSR_OrtizStages_slnFound) // means case C is not valid, other cases to be tested
		{
			stageShort2TestNow = Test_newcase(stageShort2TestNow);
			if (stageShort2TestNow == TSR_OrtizStages_invalid)
				return false;
			if (stageShort2TestNow != TSR_OrtizStages_slnFound)
			{
				stageShort2TestNow = Test_newcase(stageShort2TestNow);
				if (stageShort2TestNow == TSR_OrtizStages_invalid)
					return false;
				if (stageShort2TestNow != TSR_OrtizStages_slnFound)
				{
					if (iterN == max_iterN)
						message = "started from branch C (contact mode), didn't stay there and the new suggested branch didn't work either";
					return false;
				}
			}
		}
	}
	return true;
}

bool TSR1d_stepn_to_np1::Read_MainInputs(istream& in)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (in.eof())
		return false;
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "delT")
		{
			READ_NDOUBLE(in, buf, delT);
		}
		else if (buf == "sigmaC")
		{
			READ_NDOUBLE(in, buf, sigmaC);
		}
		else if (buf == "deltaC")
		{
			READ_NDOUBLE(in, buf, deltaC);
		}
		else if (buf == "step_n_max_delu")
		{
			READ_NDOUBLE(in, buf, step_n_max_delu);
		}
		else if (buf == "step_n_tsrStage")
		{
			in >> step_n_tsrStage;
		}
		else if (buf == "ws")
		{
			READ_NDOUBLE(in, buf, ws[SDL]);
			READ_NDOUBLE(in, buf, ws[SDR]);
		}
		else if (buf == "Zs")
		{
			READ_NDOUBLE(in, buf, Zs[SDL]);
			READ_NDOUBLE(in, buf, Zs[SDR]);
		}
		else if (buf == "us")
		{
			READ_NDOUBLE(in, buf, step_n.sides[SDL].u);
			READ_NDOUBLE(in, buf, step_n.sides[SDR].u);
		}
		else if (buf == "vs")
		{
			READ_NDOUBLE(in, buf, step_n.sides[SDL].v);
			READ_NDOUBLE(in, buf, step_n.sides[SDR].v);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	return true;
}

void TSR1d_stepn_to_np1::Write_MainInputs(ostream & out)
{
	out << "{\n";
	out << setprecision(22);
	out << "delT\t" << delT << '\n';
	out << "sigmaC\t" << sigmaC << '\n';
	out << "deltaC\t" << deltaC << '\n';
	out << "Zs\t" << Zs[SDL] << '\t' << Zs[SDR] << '\n';
	out << "step_n_tsrStage\t" << step_n_tsrStage << '\n';
	out << "step_n_max_delu\t" << step_n_max_delu << '\n';
	out << "ws\t" << ws[SDL] << '\t' << ws[SDR] << '\n';
	out << "us\t" << step_n.sides[SDL].u << '\t' << step_n.sides[SDR].u << '\n';
	out << "vs\t" << step_n.sides[SDL].v << '\t' << step_n.sides[SDR].v << '\n';
	out << "}\n";
}

void TSR1d_stepn_to_np1::Write_MainInputs_NonConvergentCase()
{
	fstream out("TestFiles/TestOrtiz1D_debug.txt", ios::out);
	Write_MainInputs(out);
}

double TSR1d_stepn_to_np1::Compute_step_np1_Solution(double& sStar, double& uLStar, double& vLStar, double& uRStar, double& vRStar, double& maxDelu, TSR_loadingStages& solutionStage, TSR_OrtizStages& stageShort)
{
	double tol_TSRFactor = 1;
	bool validSln;
	iterN = 0;
	for (iterN = 0; iterN <= max_iterN; ++iterN)
	{
		Initialize_TSR1d_stepn_to_np1(tol_TSRFactor);
		validSln = Compute_Stage_Aux();
		if (validSln)
			break;
		tol_TSRFactor *= 10.0;
	}
	if (!validSln)
	{
		cout << "message\n" << message << '\n';
		THROW_TSR("In the loop above, maybe the tolerance still needs to be loosened further\n");
	}
	stageShort = validShortStage;
	solutionStage = step_np1_tsrStage;
	TSR1d_2Sides* step_np1Ptr = &step_np1[(int)validShortStage];
	step_np1Ptr->delu = step_np1Ptr->sides[SDR].u - step_np1Ptr->sides[SDL].u;
	step_np1Ptr->delv = step_np1Ptr->sides[SDR].v - step_np1Ptr->sides[SDL].v;
	step_np1_max_delu = MAX(step_np1Ptr->delu, step_n_max_delu);
	maxDelu = step_np1_max_delu;

	sStar = step_np1Ptr->sides[SDL].s;
	uLStar = step_np1Ptr->sides[SDL].u;
	uRStar = step_np1Ptr->sides[SDR].u;
	vLStar = step_np1Ptr->sides[SDL].v;
	vRStar = step_np1Ptr->sides[SDR].v;
	return sI;
}

TSR_OrtizStages TSR1d_stepn_to_np1::Test_Stages_DelU_eqZero_A_or_E(TSR_OrtizStages stageShort)
{
	if (calcType[(int)stageShort] != -1)
	{
		if (iterN == max_iterN)
			message = "This option is tested before:Test_Stages_DelU_eqZero_A_or_E";
		return TSR_OrtizStages_invalid;
	}
	double sStar;
	TSR_OrtizStages retVal = Test_Stages_DelU_eqZero_A_or_E_Aux(sStar);
	if (retVal == TSR_OrtizStages_slnFound)
	{
		Set_StarValues_From_sStar(sStar, 0.0);
		calcType[(int)validShortStage] = 1;
		return retVal;
	}
	calcType[(int)stageShort] = 0;
	return retVal;
}

TSR_OrtizStages TSR1d_stepn_to_np1::Test_Stages_DelU_eqZero_A_or_E_Aux(double& sStar)
{
//	sStar = sI + (1-alpha)/2/alpha * Zeff * step_n.delv +  + 1/2/alpha/delt * Zeff * step_n.delu
#if TSR_BEuler
	sStar = sI;
#else
	sStar = sI + 0.5 * Zeff * step_n.delv;
#endif
	if (!maxDelUEqZero)
	{
#if TSR_BEuler
		sStar += 0.5 * step_n.delu * Zeff / delT;
#else
		sStar += step_n.delu * Zeff / delT;
#endif
		bool negStressCase = (sStar < 0.0);
#if 0
		if (!negStressCase)
		{
			double delv = 2.0 * (sI - sStar) / Zeff;
			if ((sStar < tolsc) && (delv < 0.0))
			{
				if (((step_n_max_delu >= deltaC) && (calcType[tsro_E] == 0)) || (calcType[tsro_C] == 0))
					negStressCase = true;
			}
		}
#endif
		if (negStressCase) // compressive stres is correct
		{
			step_np1_tsrStage = tsrl_damMid_delun_neg_sigI_neg;
			validShortStage = tsro_E;
			if (step_n_max_delu >= deltaC) // fully damaged
				step_np1_tsrStage = tsr1_dam1_delun_neg;
			else if (step_n_max_delu <= toldc)
				step_np1_tsrStage = tsrl_dam0_sigI_neg;
			return TSR_OrtizStages_slnFound;
		}
		if (step_n_max_delu >= deltaC) // fully damaged
		{
			step_np1_tsrStage = tsr1_dam1_delun_pos;
			validShortStage = tsro_D;
			return validShortStage;
		}
//		if (step_n_max_delu <= toldc)
//		{
//			THROW_TSR("Should not get here\n");
//		}
		// mid stage
		step_np1_tsrStage = tsrl_damMid_unloadingBranchPosDelv;
		validShortStage = tsro_C;
		return validShortStage;
	}
	else
	{
		// this is initial stage, tress is either too large (above sigmaC or not)
		if (sStar > sigmaC)
		{
//			sStar = sigmaC;
			step_n_tsrStage = tsrl_dam0_sigI_over_sigC;
			validShortStage = tsro_B;
			return validShortStage;
		}
		else if (sStar > 0.0) // tensile
		{
			step_n_tsrStage = tsrl_dam0_sigI_less_sigC;
			validShortStage = tsro_A;
//			sStar is good
			return TSR_OrtizStages_slnFound; // means the current choice is good
		}
		else // negative
		{
			step_n_tsrStage = tsrl_dam0_sigI_neg;
			validShortStage = tsro_E;
//			sStar is good
			return TSR_OrtizStages_slnFound; // means the current choice is good
		}
	}
}

TSR_OrtizStages TSR1d_stepn_to_np1::Test_Stage_B()
{
	TSR_OrtizStages stageShort = tsro_B;
	if (calcType[(int)stageShort] != -1)
	{
		if (iterN == max_iterN)
			message = "This option is tested before:Test_Stage_B";
		return TSR_OrtizStages_invalid;
	}
	// computin delu
	// numerator: Zeff * (step_n.delu + (1 - alpha) * delT * step_n.delv) + 2 * alpha * delT * (sI - sigmaC)
	// denominator: Zeff * deltaC - 2 * alpha * sigmaC * delT
#if TSR_BEuler
	double delu_numerator = Zeff * step_n.delu  + 2.0 * delT * (sI - sigmaC);
	double delu_denominator = Zeff * deltaC - 2.0 * sigmaC * delT;
#else
	double delu_numerator = Zeff * (step_n.delu + 0.5 * delT * step_n.delv) + delT * (sI - sigmaC);
	double delu_denominator = Zeff * deltaC - sigmaC * delT;
#endif
	if (delu_denominator < 0)
	{
		cout << "tauCohesive\t" << Zeff * deltaC / sigmaC << '\n';
		cout << "delT\t" << delT << '\n';
		THROW_TSR("delT > tauCohesve, reduce time step (half tauC for backward euler, equal to tauC for central)\n");
	}
	double delu_factor = computeRatio(delu_numerator, delu_denominator);
	double delu = delu_factor * deltaC;
	if (deltaC - delu < -toldc)
	{
		calcType[(int)stageShort] = 0;
		return tsro_D;
	}
	if (delu_factor < DBL_MIN)
	{
		if ((calcType[tsro_A] == 0) && (iterN == max_iterN))
			delu_factor = 0;
		else
		{
			calcType[(int)stageShort] = 0;
			return tsro_E;
		}
	}
	double sStar, delv;
	ComputeStarValues_From_delu(delu, sStar, delv);
	double tolv = 0; // tol_TSR * sigmaC / Zeff;
	if (delv < -tolv) // unloading -> should use he unloading branch
	{
		calcType[(int)stageShort] = 0;
		return tsro_C;
	}
	// valid solution from here on
	validShortStage = tsro_B;
	calcType[(int)stageShort] = 1;
	step_np1_tsrStage = tsrl_damMid_loadingBranch;
	Set_StarValues_From_sStar(sStar, delu);
	return TSR_OrtizStages_slnFound;
}

TSR_OrtizStages TSR1d_stepn_to_np1::Test_Stage_C()
{
	TSR_OrtizStages stageShort = tsro_C;
	if (calcType[(int)stageShort] != -1)
	{
		if (iterN == max_iterN)
			message = "This option is tested before:Test_Stage_C";
		return TSR_OrtizStages_invalid;
	}

	// computin delu
	//delu_numerator = Zeff * (step_n.delu + (1 - alpha) * delT * step_n.delv) + 2 * alpha * delT * sI;
	// delu_denominator = Zeff * step_n_max_delu + 2 * alpha * step_n_max_sigma * delT;

#if TSR_BEuler
	double delu_numerator = Zeff * step_n.delu + 2.0 * delT * sI;
	double delu_denominator = Zeff * step_n_max_delu + 2.0 * step_n_max_sigma * delT;
#else
	double delu_numerator = Zeff * (step_n.delu + 0.5 * delT * step_n.delv) + delT * sI;
	double delu_denominator = Zeff * step_n_max_delu + step_n_max_sigma * delT;
#endif
	double delu_factor = computeRatio(delu_numerator, delu_denominator);
	double delu = delu_factor * step_n_max_delu;
	if (1 - delu_factor < -tol_TSR) // it's back to loading branch
	{
		step_np1_tsrStage = tsrl_damMid_unloadingBranchPosDelv;
		calcType[(int)stageShort] = 0;
		if (deltaC - delu > -toldc)
		{
			if (calcType[(int)tsro_B] == -1)
				return tsro_B;
			else if (1 - delu_factor < 0.1 * tol_TSR)
			{
				cout << "delu_factor\t" << delu_factor << '\n';
				cout << "step_n_max_delu\t" << step_n_max_delu << '\n';
				cout << "delu\t" << delu << '\n';
				cout << "deltaC\t" << deltaC << '\n';
				{
					if (iterN == max_iterN)
						message = "This may be a tolerance issue and may adjust the tolerances";
					return TSR_OrtizStages_invalid;
				}
			}
		}
		else
			return tsro_D;
	}
	double sStar, delv;
	ComputeStarValues_From_delu(delu, sStar, delv);
	if (delv < 0.0) // unloading -> should use he unloading branch
		step_np1_tsrStage = tsrl_damMid_unloadingBranchNegDelv;
	else
		step_np1_tsrStage = tsrl_damMid_unloadingBranchPosDelv;

	if ((delu_factor < 0) && (calcType[tsro_E] == -1))// going to contact branch
	{
		TSR_OrtizStages retVal = Test_Stages_DelU_eqZero_A_or_E(tsro_E);
		if ((retVal == TSR_OrtizStages_slnFound) && (validShortStage == tsro_E))
			return retVal;
		else
		{
//			delu_factor = 0.0;
//			delu = 0.0;
		}
	}
	// valid solution from here on
	validShortStage = tsro_C;
	calcType[(int)stageShort] = 1;
	Set_StarValues_From_sStar(sStar, delu);
	return TSR_OrtizStages_slnFound;
}

TSR_OrtizStages TSR1d_stepn_to_np1::Test_Stage_D()
{
	TSR_OrtizStages stageShort = tsro_D;
	if (calcType[(int)stageShort] != -1)
	{
		if (iterN == max_iterN)
			message = "This option is tested before:Test_Stage_D";
		return TSR_OrtizStages_invalid;
	}


	// computin delu
	// double delu = step_n.delu + (1-alpha) * delT * step_n.delv + 2 * alpha * delT * sI / Zeff;
#if TSR_BEuler
	double delu = step_n.delu + 2.0 * delT * sI / Zeff;
#else
	double delu = step_n.delu + 0.5 * delT * step_n.delv + delT * sI / Zeff;
#endif
	if (delu < 0.0)
	{
		step_np1_tsrStage = tsr1_dam1_delun_neg;
		calcType[(int)stageShort] = 0;
		return tsro_E;
	}
	step_np1_tsrStage = tsr1_dam1_delun_pos;
	calcType[(int)stageShort] = 1;

	double sStar = 0.0;
	// valid solution from here on
	validShortStage = tsro_D;
	calcType[(int)stageShort] = 1;
	Set_StarValues_From_sStar(sStar, delu);
	return TSR_OrtizStages_slnFound;
}

void TSR1d_stepn_to_np1::Set_StarValues_From_sStar(double sStar, double delu)
{
	TSR1d_2Sides* step_np1Ptr = &step_np1[(int)validShortStage];
	step_np1Ptr->sides[SDL].s = sStar;
	step_np1Ptr->sides[SDR].s = sStar;

	step_np1Ptr->sides[SDR].v = (ws[SDR] - step_np1Ptr->sides[SDR].s) / Zs[SDR];
	step_np1Ptr->sides[SDL].v = (step_np1Ptr->sides[SDL].s - ws[SDL]) / Zs[SDL];
	// the velocities should be equal

#if TSR_BEuler
	step_np1Ptr->sides[SDL].u = step_n.sides[SDL].u + delT * step_np1Ptr->sides[SDL].v;
	step_np1Ptr->sides[SDR].u = step_n.sides[SDR].u + delT * step_np1Ptr->sides[SDR].v;
#else // central
	double fact = 0.5 * delT;
	step_np1Ptr->sides[SDL].u = step_n.sides[SDL].u + fact * (step_n.sides[SDL].v + step_np1Ptr->sides[SDL].v);
	step_np1Ptr->sides[SDR].u = step_n.sides[SDR].u + fact * (step_n.sides[SDR].v + step_np1Ptr->sides[SDR].v);
#endif
	double deluComputed = step_np1Ptr->sides[SDR].u - step_np1Ptr->sides[SDL].u;
	if (DoublesAreEqual(delu, deluComputed) == false)
	{
		cout << "delu\t" << delu << '\n';
		cout << "deluComputed\t" << deluComputed << '\n';
		THROW_TSR("delu's not equal\n");
	}
}

void TSR1d_stepn_to_np1::ComputeStarValues_From_delu(double delu, double& sStar, double& delv)
{
//	delv = 1 / alpha / delT * (delu - step_n.delu) - step_n.delv * (1 - alpha) / alpha;
#if TSR_BEuler
	delv = 1.0 / delT * (delu - step_n.delu);
#else
	delv = 2.0 / delT * (delu - step_n.delu) - step_n.delv;
#endif
	sStar = sI - 0.5 * Zeff * delv;
}

TSR_OrtizStages TSR1d_stepn_to_np1::Test_newcase(TSR_OrtizStages stageShort2Test)
{
	if (stageShort2Test == tsro_B)
		return Test_Stage_B();
	if (stageShort2Test == tsro_C)
		return Test_Stage_C();
	if (stageShort2Test == tsro_D)
		return Test_Stage_D();
	if (stageShort2Test == tsro_E)
		return Test_Stages_DelU_eqZero_A_or_E(tsro_E);
	if (stageShort2Test == tsro_A)
		return Test_Stages_DelU_eqZero_A_or_E(tsro_A);
	THROW_TSR("invalid branch\n");
}


int Test_TSR_Ortiz_1D(string nameWOExtIn)
{
	string input_file_name = nameWOExtIn + ".txt";
	fstream in(input_file_name.c_str(), ios::in);
	if (!in.is_open())
	{
		nameWOExtIn = "TestFiles/TestOrtiz1D";
		string input_file_name = nameWOExtIn + ".txt";
		in.open(input_file_name.c_str(), ios::in);
		if (!in.is_open())
		{
			THROW("Cannot open the input file name\n");
		}
	}
	string output_file_name = nameWOExtIn + ".out";
	fstream out(output_file_name.c_str(), ios::out);

	int cntr = 0;
	while (true)
	{
		TSR1d_stepn_to_np1 tsr1D;
		if (tsr1D.Read_MainInputs(in) == false)
			break;
		double sI, sStar, uLStar, vLStar, uRStar, vRStar, maxDelu;
		TSR_loadingStages solutionStage;
		TSR_OrtizStages stageShort;
		sI = tsr1D.Compute_step_np1_Solution(sStar, uLStar, vLStar, uRStar, vRStar, maxDelu, solutionStage, stageShort);
		out << setprecision(22);
		out << "cntr\t" << cntr++;
		out << "\tsI\t" << sI;
		out << "\tsStar\t" << sStar;
		out << "\tuLStar\t" << uLStar;
		out << "\tuRStar\t" << uRStar;
		out << "\tdelu\t" << (uRStar - uLStar);
		out << "\tvLStar\t" << vLStar;
		out << "\tvRStar\t" << vRStar;
		out << "\tdelv\t" << (vRStar - vLStar);
		out << "\tmaxDelu\t" << maxDelu;
		out << "\tstageShort\t" << stageShort;
		out << "\tsolutionStage\t" << solutionStage;
		out << endl;
		out.flush();
	}
	cout << "cntr\t" << cntr << '\n';
	getchar();
	return cntr;
}

TSR_XuNeedleman::TSR_XuNeedleman()
{
	a = 10.0;
	Z = 1.0;
	deltaC = 1.0;
	sigmaC = 1.0;
	delTFactor = 0.0001;
	sigmaCFactor4Zero = 0.01;
}

void TSR_XuNeedleman::Compute(ostream* outPtr)
{
	tauC = deltaC * Z / sigmaC;
	delT = tauC * delTFactor;
	double sigTol = sigmaCFactor4Zero * sigmaC;
	// sigma(delta) + Z/2 * deltaDot = at -> 
	// deltaDot = At - B delta exp(1.0 - delta/deltaC)
	double ZInv2 = 2.0 / Z;
	double deltaCInv = 1.0 / deltaC, deln;
	double A = ZInv2 * a;
	double B = ZInv2 * sigmaC * deltaCInv;
	double time = 0.0;
	double delta = 0.0;
	double sigma = 0.0;
	double deltaDot = 0.0;
	bool b_print = (outPtr != NULL);
	if (b_print)
	{
		*outPtr << "time\tdelta\tsigma\tdeltaDot\n";
		*outPtr << time << "\t" << delta << "\t" << sigma << "\t" << deltaDot << "\n";
	}
	double h = delT, hd2 = 0.5 * delT, h6 = h / 6.0;
	double k1, k2, k3, k4, del, t;
	bool deltaCPassed = false;
	bool cont = true;
	while (cont)
	{
		t = time;
		del = delta;
		k1 = A * t - B * del * exp(1.0 - del * deltaCInv);
		t = time + hd2;
		del = delta + 0.5 * k1;
		k2 = A * t - B * del * exp(1.0 - del * deltaCInv);
		del = delta + 0.5 * k2;
		k3 = A * t - B * del * exp(1.0 - del * deltaCInv);
		t = time + h;
		del = delta + k3;
		k4 = A * t - B * del * exp(1.0 - del * deltaCInv);
		time += h;
		delta += h6 * (k1 + 2.0 * (k2 + k3) + k4);
		deln = delta * deltaCInv;
		sigma = sigmaC * deln * exp(1.0 - deln);
		deltaDot = ZInv2 * (a * time - sigma);
		if ((!deltaCPassed) && (deln >= 1.000))
		{
			deltaCPassed = true;
			tSigmaMax = time;
		}
		if (b_print)
			*outPtr << time << "\t" << delta << "\t" << sigma << "\t" << deltaDot << "\n";
		if ((deltaCPassed) && (sigma < sigTol))
		{
			delnSigmaZero = deln;
			tSigmaZero = time;
			cont = false;
		}
	}
	if (b_print)
		outPtr->flush();
}

void TSR_XuNeedleman::OutMain(ostream & out)
{
	out << "tSigmaMax\t" << tSigmaMax;
	out << "\tdelnSigmaZero\t" << delnSigmaZero;
	out << "\ttSigmaZero\t" << tSigmaZero;
	out << '\n';
}


