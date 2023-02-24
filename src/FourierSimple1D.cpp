#include "FourierSimple1D.h"
#include "commonMacros.h"
#include "globalFunctions.h"
#include "SimpleConfigMaker.h"

string getName(LoadModeType dat)
{
	if (dat == LoadModeType_none)
		return "none";
	if (dat == lmt_Incident)
		return "Incident";
	if (dat == lmt_Impact)
		return "Impact";
	if (dat == lmt_Neumann)
		return "Neumann";
	if (dat == lmt_Dirichlet)
		return "Dirichlet";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, LoadModeType& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= LoadModeType_SIZE)
			THROW("too large of a number\n");
		typeVal = (LoadModeType)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < LoadModeType_SIZE; ++i)
	{
		typeVal = (LoadModeType)i; // casting integer to LoadModeType, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading LoadModeType\n");
}

//operator for output
ostream& operator<<(ostream& out, LoadModeType dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, LoadModeType& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

istream& operator>>(istream& in, TD_FD& dat)
{
	string buf, buf2;
	READ_NSTRING(in, buf2, buf);
	if (buf != "{")
		THROW("SpatialAllDirsVariationFunction starts with '{'!\n");

	while (buf != "}")
	{
		READ_NSTRING(in, buf2, buf);

		if (buf == "signalType")
		{
			in >> dat.signalType;
		}
		else if (buf == "modelParas")
		{
			ReadVectorDouble(in, dat.modelParas);
			dat.b_modelParas = true;
		}
		else if (buf == "modelFlags")
		{
			ReadVectorInteger(in, dat.modelFlags);
		}
		else if (buf == "relw4Dirac_w")
		{
			READ_NDOUBLE(in, buf2, dat.relw4Dirac_w);
		}
		else if (buf == "tmin")
		{
			READ_NDOUBLE(in, buf2, dat.tmin);
			dat.b_tmin = true;
			++dat.n_t_lims;
		}
		else if (buf == "tmax")
		{
			READ_NDOUBLE(in, buf2, dat.tmax);
			dat.b_tmax = true;
			++dat.n_t_lims;
		}
		else if (buf == "delt")
		{
			READ_NDOUBLE(in, buf2, dat.delt);
			dat.b_delt = true;
			++dat.n_t_lims;
		}
		else if (buf == "wmin")
		{
			READ_NDOUBLE(in, buf2, dat.wmin);
			dat.b_wmin = true;
			++dat.n_w_lims;
		}
		else if (buf == "wmax")
		{
			READ_NDOUBLE(in, buf2, dat.wmax);
			dat.b_wmax = true;
			++dat.n_w_lims;
		}
		else if (buf == "delw")
		{
			READ_NDOUBLE(in, buf2, dat.delw);
			dat.b_delw = true;
			++dat.n_w_lims;
		}
		else if (buf == "uniform_ts")
		{
			READ_NBOOL(in, buf2, dat.uniform_ts);
		}
	}
	string key;	map<string, string>* mpPtr;
	double value = -1;
	for (unsigned int i = 0; i < dat.modelParas.size(); ++i)
	{
		string str_i;
		toString(i, str_i);
		key = "tdPara" + str_i;
		if (Find_Version_Value(key, value, mpPtr) == true)
			 dat.modelParas [i] = value;
	}
	dat.Initialize_AfterRead();
	return in;
}

TD_FD::TD_FD()
{
	relw4Dirac_w = 0.5;
	uniform_ts = false;
	hasComputedFT = false;
	tmin = 0.0;
	tmax = 1.0;
	delt = 0.01;

	wmin = 0.0;
	wmax = 1.0;
	delw = 0.01;

	modelParas.resize(10);
	fill(modelParas.begin(), modelParas.end(), 0.0);
	modelFlags.resize(3);
	fill(modelFlags.begin(), modelFlags.end(), 0);

	b_modelParas = false;
	b_tmin = false;
	b_tmax = false;
	b_delt = false;
	n_t_lims = 0;
	b_wmin = false;
	b_wmax = false;
	b_delw = false;
	n_w_lims = 0;

	timeScale = 0.0;
}

void TD_FD::Set_ts_ws_external(vector<double>& modelParasIn, vector<double>& tsIn, int& uniform_ts_i, double& wminIn, double& wmaxIn, double& delwIn, bool initialize_w)
{
	if (tsIn.size() > 0)
	{
		numt = tsIn.size();
		tmin = tsIn[0];
		tmax = tsIn[numt - 1];
		if (uniform_ts_i >= 0)
			uniform_ts = (bool)uniform_ts_i;
	}
	else if (n_t_lims == 3)
	{
		numt = (int)floor((tmax - tmin) / delt + 1e-6);
		ts.resize(numt);
		for (int i = 0; i < numt; ++i)
			ts[i] = tmin + delt * i;
		uniform_ts = true;
		tsIn = ts;
		uniform_ts_i = (int)uniform_ts;
	}
	if (initialize_w)
	{
		if (wminIn > wmaxIn)
		{
			wmin = wminIn;
			wmax = wmaxIn;
			delw = delwIn;
		}
		else
		{
			if (b_wmin)
				wminIn = wmin;
			if (b_wmax)
				wmaxIn = wmax;
			if (b_delw)
				delwIn = delw;
		}
	}
	if (b_modelParas)
		modelParasIn = modelParas;
	else
		modelParas = modelParasIn;
}

void TD_FD::Initialize_AfterRead()
{
	if (signalType == td_sine_TD_Gaussian)
	{
		double n = modelParas[ptg_n];
		double wminV = wmin, wmaxV = wmax; // modelParas[ptg_wmin], wmaxV = modelParas[ptg_wmax];
		//double wmM = fabs(wminV) + fabs(wmaxV);
		bool wmMProveded = (b_wmin && b_wmax);
		if (wmMProveded)
		{
			modelParas[ptg_w0] = 0.5 * (wmaxV + wminV);
			modelParas[ptg_sigma] = 2.0 / (wmaxV - wminV);
		}
		else
		{
			double w0 = modelParas[ptg_w0];
			double sigma = modelParas[ptg_sigma];
			wmin = w0 - 1.0 / sigma;
			wmax = w0 + 1.0 / sigma;
			//			modelParas[ptg_wmin] = w0 - 1.0 / sigma;
			//			modelParas[ptg_wmax] = w0 + 1.0 / sigma;
		}
		if (n > 1e-7)
			modelParas[ptg_t0] = 2.0 * n * modelParas[ptg_sigma];
		else
			modelParas[ptg_n] = 0.5 * modelParas[ptg_t0] / modelParas[ptg_sigma];
		timeScale = 2 * PI / modelParas[ptg_w0];
	}
	else if (signalType == td_Friedlander)
	{
		hasComputedFT = false;
		if ((signalType == td_sine_TD_Gaussian) || (signalType == td_UnitSquare))
			hasComputedFT = true;
		else if (signalType == td_Friedlander)
		{
			hasComputedFT = true;
			int sz = modelParas.size();
			if (sz > frd_tramp)
			{
				double tramp = modelParas[frd_tramp];
				if (tramp > 1e-12)
					hasComputedFT = false;
				else if (sz > frd_tflat)
				{
					double tflat = modelParas[frd_tflat];
					if (tflat > 1e-12)
						hasComputedFT = false;
				}
			}
			timeScale = modelParas[frd_td];
		}
		else if (signalType == td_UnitSquare)
		{
			timeScale = modelParas[frd_td];
		}
		else
		{
			cout << "signalType\t" << signalType << '\n';
			THROW("option not implemented\n");
		}
		if (modelFlags.size() == 0)
		{
			modelFlags[0] = fldt_Friedlander;
		}
		if (modelFlags[0] != fldt_Friedlander)
			hasComputedFT = false;
	}
}

void TD_FD::Initialize_TD_FD(bool initialize_FD)
{
	Compute_TD_Values();
	if (initialize_FD)
		Compute_FD_Values();
}

void TD_FD::Compute_TD_Values()
{
	numt = ts.size();
	if (numt == 0)
	{
		numt = (int)floor((tmax - tmin) / delt + 1e-6);
		ts.resize(numt);
		for (int i = 0; i < numt; ++i)
			ts[i] = tmin + delt * i;
		uniform_ts = true;
	}
	if (td_vals.size() == 0)
	{
		if (signalType == td_gen)
		{
			THROW("TD vlaues must have been specified beforehand for general loading case\n");
		}
		td_vals.resize(numt);
		td_DTs.resize(numt);
		fill(td_DTs.begin(), td_DTs.end(), 0.0);
		for (int i = 0; i < numt; ++i)
			td_vals[i] = ComputeSpecific_TD_Normalized_Value(ts[i], td_DTs[i], true);
	}
}

void TD_FD::Compute_FD_Values()
{
	numw = ws.size();
	if (numw == 0)
	{
		numw = (int)floor((wmax - wmin) / delw + 1e-6);
		ws.resize(numw);
		for (int i = 0; i < numw; ++i)
			ws[i] = wmin + delw * i;
	}
	if (fd_vals.size() != numw)
	{
		fd_vals.resize(numw);
		if (hasComputedFT == false)
			ComputeFourierTransformValues_UniformTime(ws, ts, td_vals, fd_vals, uniform_ts);
		else
			for (int i = 0; i < numw; ++i)
				fd_vals[i] = ComputeSpecific_FD_Value(ws[i]);
	}
}

double TD_FD::ComputeSpecific_TD_Normalized_Value(double t, double& dValdT, bool compute_dValdT)
{
	if (signalType == td_sine_TD_Gaussian)
	{
		double t0 = modelParas[ptg_t0];
		double w0 = modelParas[ptg_w0];
		double sigma = modelParas[ptg_sigma];

		//		double delt = t - t0;
		//		double tx = x / c[lcr_l][wmode];
		//		double tp = delt - tx;
		t = 25;
		double tp = t - t0;
		double tcs = tp;
		double pwr = tp / sigma;
		pwr = -0.5 * pwr * pwr;
		double expTerm = exp(pwr);
		double val = sin(w0 * tcs) * expTerm;
		if (compute_dValdT)
			dValdT = (w0 * cos(w0 * tcs) - tp / sigma / sigma * sin(w0 * tcs)) * expTerm;
		return val;
	}
	else if (signalType == td_Friedlander)
	{
		double td = modelParas[frd_td];
		double b = modelParas[frd_b];
		double P = 1.0;
		double p0 = modelParas[frd_p02P];
		double tr = modelParas[frd_tr];
		double tramp = 0.0;
		double tflat = 0.0;
		int sz = modelParas.size();
		if (sz > frd_tramp)
		{
			tramp = modelParas[frd_tramp];
			if (sz > frd_tflat)
				tflat = modelParas[frd_tflat];
		}
		double tau = (t - tr);
		if (tau < 0)
		{
			dValdT = 0.0;
			return p0;
		}
		if (tau < tramp)
		{
			dValdT = (P - p0) / tramp;
			return p0 + dValdT * tau;
		}
		tau -= tramp;
		if (tau < tflat)
		{
			dValdT = 0;
			return P;
		}
		tau -= tflat;
		tau /= td;
		double val;
		if (modelFlags[0] == fldt_Friedlander)
		{
			double expPart = exp(-b * tau);
			val = p0 + P * (1 - tau) * expPart;
			if (compute_dValdT)
			{
				dValdT = -P / td * ((1 - tau) * b + 1) * expPart;
			}
		}
		else if (modelFlags[0] == fldt_linear)
		{
			if (tau < 1)
			{
				val = P * (1 - tau);
				if (compute_dValdT)
				{
					dValdT = -P / td;
				}
			}
			else
			{
				val = 0.0;
				dValdT = 0.0;
			}
		}
		else if (modelFlags[0] == lfdt_exponential)
		{
			val = P * exp(-tau);
			if (compute_dValdT)
			{
				dValdT = -val / td;
			}
		}
		else
		{
			THROW("Invalid option\n");
		}
		return val;
	}
	else if (signalType == td_UnitSquare)
	{
		double td = modelParas[frd_td];
		double P = 1.0;
		dValdT = 0.0;
		if (t <= td)
			return P;
		return 0;
	}
	THROW("Invalid option\n");
}

double TD_FD::ComputeSpecific_TD_Value(double t)
{
	double dValdT;
	return ComputeSpecific_TD_Normalized_Value(t, dValdT, false);
}

Dcomplex TD_FD::ComputeSpecific_FD_Value(double w)
{
	if (signalType == td_sine_TD_Gaussian)
	{
		double t0 = modelParas[ptg_t0];
		double w0 = modelParas[ptg_w0];
		double sigma = modelParas[ptg_sigma];

		//	if (xAtInterface)
		//		x = 0.0;
		double tt = t0; //x / c[lcr_l][wmode] + t0;
		double sigmawpw0_2 = sigma * (w + w0);
		sigmawpw0_2 = sigmawpw0_2 * sigmawpw0_2;

		double sigmawmw0_2 = sigma * (w - w0);
		sigmawmw0_2 = sigmawmw0_2 * sigmawmw0_2;

		double termp = exp(-0.5 * sigmawpw0_2), termm = exp(-0.5 * sigmawmw0_2);
		double factor = 0.5 * sigma * (termp - termm);

		double s = sin(w * tt), c = cos(w * tt);
		Dcomplex tmp(s * factor, c * factor);
		return tmp;
	}
	else if (signalType == td_Friedlander)
	{
		double td = modelParas[frd_td];
		double b = modelParas[frd_b];
		double P = 1.0;
		double p0 = modelParas[frd_p02P];
		double tr = modelParas[frd_tr];

		double wdirac_range = relw4Dirac_w * delw;
		Dcomplex I(0.0, 1.0);
		Dcomplex wTilde = td * w;

		double dydx;
		double dirac_w = getInfinitlySmoothDeltaDirac(w, dydx, wdirac_range, false);
		Dcomplex factor = P * td * exp(-I * w * tr);
		Dcomplex tmp = 1.0 / (b + I * wTilde);
		Dcomplex val = factor * tmp * (1.0 - tmp) / sqrt(2.0 * PI) + p0 * dirac_w;
		return val;
	}
	else if (signalType == td_UnitSquare)
	{
		double P = 1.0;
		double td = modelParas[frd_td];
		Dcomplex I(0.0, 1.0);
		Dcomplex wTilde = td * w;

		double factor = computeRatio(P / sqrt(2.0 * PI), w);
		Dcomplex val = factor * (sin(wTilde) - I * (1.0 - cos(wTilde)));
		return val;
	}
	THROW("Invalid option\n");
}

void TD_FD::Size_TD()
{
	numt = (int)floor((tmax - tmin) / delt + 1e-6);
	ts.resize(numt);
	td_vals.resize(numt);
	td_DTs.resize(numt);
	for (int i = 0; i < numt; ++i)
		ts[i] = tmin + delt * i;
	fill(td_vals.begin(), td_vals.end(), 0.0);
	fill(td_DTs.begin(), td_DTs.end(), 0.0);
}

void TD_FD::Compute_Inegral_load_load2(double timeMax, double delT, double & integral_loadOut, double & integral_load2Out)
{
	int n = (int)ceil((timeMax / delT) / 2) * 2;
	delT = timeMax / n;
	double t = 0.0;
	integral_loadOut = 0.0;
	integral_load2Out = 0.0;
	double val = ComputeSpecific_TD_Value(t);
	integral_loadOut += val;
	integral_load2Out += val * val;
	val = ComputeSpecific_TD_Value(timeMax);
	integral_loadOut += val;
	integral_load2Out += val * val;

	t = delT;
	int nm1 = n - 1;
	for (int i = 1; i < n; i += 2)
	{
		val = ComputeSpecific_TD_Value(t);
		integral_loadOut += 4.0 * val;
		integral_load2Out += 4.0 * val * val;
		if (i == nm1)
			break;
		t += delT;
		val = ComputeSpecific_TD_Value(t);
		integral_loadOut += 2.0 * val;
		integral_load2Out += 2.0 * val * val;
		t += delT;
	}
	double fact = delT / 3.0;
	integral_loadOut *= fact;
	integral_load2Out *= fact;
}

ostream & operator<<(ostream & out, const TD_FD & dat)
{
	out << "signalType\t" << dat.signalType << '\n';
	out << "modelParas\n" << dat.modelParas << '\n';
	out << "relw4Dirac_w\n" << dat.relw4Dirac_w << '\t';
	out << "tmin\n" << dat.tmin << '\t';
	out << "tmax\n" << dat.tmax << '\t';
	out << "delt\n" << dat.delt << '\t';
	out << "wmin\n" << dat.wmin << '\t';
	out << "wmax\n" << dat.wmax << '\t';
	out << "delw\n" << dat.delw << '\n';

	out << "numt\t" << dat.numt << '\n';
	out << "numw\t" << dat.numw << '\n';
	out << "ts\n" << dat.ts << '\n';
	out << "ws\n" << dat.ws << '\n';
	out << "td_vals\n" << dat.td_vals << '\n';
	out << "td_DTs\n" << dat.td_DTs << '\n';
	out << "fd_vals\n" << dat.fd_vals.size() << '\n';
	for (unsigned int i = 0; i < dat.fd_vals.size(); ++i)
		out << dat.fd_vals[i] << '\n';
	return out;
}

Dcomplex ComputeFourierTransform_UniformTime(double w, const vector<double>& ts, vector<double>& ft, bool uniform_ts)
{
	int sz = ts.size();
	vector<double> re(sz), im(sz);
	for (int i = 0; i < sz; ++i)
	{
		double wt = w * ts[i];
		re[i] =  ft[i] * cos(wt);
		im[i] = -ft[i] * sin(wt);
	}
	vector<double> re_Int(sz), im_Int(sz);
	Integrate(ts, re, re_Int, sz, uniform_ts);
	Integrate(ts, im, im_Int, sz, uniform_ts);
	Dcomplex fhat_w = 1.0 / sqrt(2.0 * PI) * Dcomplex(re_Int[sz - 1], im_Int[sz - 1]);
	return fhat_w;
}

void ComputeFourierTransformValues_UniformTime(const vector<double>& ws, const vector<double>& ts, vector<double>& ft, vector<Dcomplex>& fd_vals, bool uniform_ts)
{
	int numw = ws.size();
	fd_vals.resize(numw);
	for (int i = 0; i < numw; ++i)
		fd_vals[i] = ComputeFourierTransform_UniformTime(ws[i], ts, ft, uniform_ts);
}

string getName(TD_signalT dat)
{
	if (dat == td_gen)
		return "gen";
	if (dat == td_sine_TD_Gaussian)
		return "sine_TD_Gaussian";
	if (dat == td_Friedlander)
		return "Friedlander";
	if (dat == td_UnitSquare)
		return "UnitSquare";
	cout << "dat\t" << (int)dat << '\n';
	THROW("invalid type\n");
}

bool fromName(const string& name, TD_signalT& dat)
{
	string nameTest;
	for (int i = 0; i < td_SIZE; ++i)
	{
		nameTest = getName((TD_signalT)i);
		if (nameTest == name)
		{
			dat = (TD_signalT)i;
			return true;
		}
	}
	THROW("invalid type\n");
	//	return false;
}

ostream& operator<<(ostream& out, TD_signalT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

istream& operator>>(istream& in, TD_signalT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, buf);
	fromName(buf, dat);
	return in;
}

void Test_TD_FD_Class(TD_signalT signalType)
{
	TD_FD tdfd;
	tdfd.signalType = signalType;
	tdfd.tmin = 0.0;
	tdfd.tmax = 50.0;
	tdfd.delt = 0.01;
	tdfd.wmin = -2.0;
	tdfd.wmax = 15.0;
	tdfd.delw = 0.05;
	tdfd.relw4Dirac_w = 0.001;

	if (signalType == td_sine_TD_Gaussian)
	{
		double t0 = 2.0;
		double w0 = 4.0;
		double sigma = 1.0;

		tdfd.modelParas[ptg_t0] = t0;
		tdfd.modelParas[ptg_w0] = w0;
		tdfd.modelParas[ptg_sigma] = sigma;
	}
	else if ((signalType == td_Friedlander) || (signalType == td_UnitSquare))
	{
		double td = 20.0;
		double b = 3.0;
		double P = 4.0;
		double p0 = 0.0;
		double tr = 0.0; // 2.0;
		tdfd.modelParas[frd_td] = td;
		tdfd.modelParas[frd_b] = b;
		tdfd.modelParas[frd_p02P] = p0 / P;
		tdfd.modelParas[frd_tr] = tr;
		tdfd.modelParas[frd_tramp] = 5.0;
		tdfd.modelParas[frd_tflat] = 15.0;

		if (tdfd.modelFlags.size() == 0)
			tdfd.modelFlags.resize(1);
		tdfd.modelFlags[0] = fldt_linear;
	}
	else
	{
		tdfd.Size_TD();

		for (int i = 0; i < tdfd.numt; ++i)
		{
			double t = tdfd.ts[i];
			// square pulse
			if ((t >= -1e-6) && ((t - 1.0) < 1e-6))
				tdfd.td_vals[i] = 1.0;
		}
		tdfd.uniform_ts = true;
	}
	tdfd.Initialize_TD_FD(true);
	vector<Dcomplex> fd_vals_integrated;
	ComputeFourierTransformValues_UniformTime(tdfd.ws, tdfd.ts, tdfd.td_vals, fd_vals_integrated, tdfd.uniform_ts);

	fstream out("TestFiles/test_TDFD_members.txt", ios::out);
	out << tdfd << '\n';
	out.close();
	out.open("TestFiles/test_TDFD_FDTest.txt", ios::out);
	for (unsigned int i = 0; i < tdfd.ws.size(); ++i)
	{
		double w = tdfd.ws[i];
		out << w << '\t' << tdfd.fd_vals[i].real() << '\t' << tdfd.fd_vals[i].imag() << '\t' << fd_vals_integrated[i].real() << '\t' << fd_vals_integrated[i].imag() << '\n';
	}
	out.close();
	out.open("TestFiles/test_TDFD_TDTest.txt", ios::out);
	for (unsigned int i = 0; i < tdfd.ts.size(); ++i)
	{
		double t = tdfd.ts[i];
		out << t << '\t' << tdfd.td_vals[i] << '\n';
	}
}
