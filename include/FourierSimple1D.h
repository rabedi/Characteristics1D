#ifndef FOURIER_SIMPLE_1D__H
#define FOURIER_SIMPLE_1D__H

#include <iostream>
#include <fstream>

#include <vector>
#include <cstdlib>
using namespace std;

#include "Dims.h"
#include "commonFunctions.h"

// td_gen is a general time domain signal (for our applications nonzero for t >= 0)
// td_sine_TD_Gaussian:
// f(t) = sine(omega_0 t) exp(-[(t - t0)/sigma]^2) -> fhat(omega) = 1/2 sigma exp(-omega t0) [exp(0.5 sigma^2(omega + omega_0)^2) - exp(0.5 sigma^2(omega - omega_0)^2) ]

// td_Friedlander
// equation Appendix (2), Mitchell_2016_Pandolfi_Ortiz_Effect of Brittle Fracture in a Metaconcrete Slab under Shock Loading.pdf 
// f(t) = p0 + P(1 - tau) exp(-b tau), tau = (t - tr)/Td for tau > 0, f(t) = p0 for tau < 0 
// fhat(omega) = p0 delta(omega) + P T_d exp(-i omega tr)[1/omegap/sqrt(2 pi)(1/omegap - i) + sqrt(pi/2)(delta(omegap) - i delta'(omegap)]
// td_UnitSquare -> 1 from 0 to Td, zero later
// omegap = T_d * w - i b
typedef enum {td_gen, td_sine_TD_Gaussian, td_Friedlander, td_UnitSquare, td_SIZE} TD_signalT;

// originally had Friedlander decay, but can consider other types too
typedef enum {fldt_Friedlander, fldt_linear, lfdt_exponential, FriedlanderGrpDecayT_SIZE} FriedlanderGrpDecayT;

// parameters for td_sine_TD_Gaussian
typedef enum {ptg_n, ptg_t0, ptg_w0, ptg_sigma, ptg_SIZE} sine_TD_GaussianPara; // ptg_wmin, ptg_wmax, 
// parameters for Friedlander, and UnitSquare 
typedef enum {frd_td, frd_b, frd_p02P, frd_tr, frd_tramp, frd_tflat, frd_SIZE} FriedlanderPara;

typedef enum { LoadModeType_none = -1, lmt_Incident, lmt_Impact, lmt_Neumann, lmt_Dirichlet, LoadModeType_SIZE} LoadModeType;

string getName(LoadModeType dat);
void name2Type(string& name, LoadModeType& typeVal);
ostream& operator<<(ostream& out, LoadModeType dat);
istream& operator>>(istream& in, LoadModeType& dat);

// Fourier transform defined by:
// fhat(omega) = 1/sqrt(2 pi) \int_{-infty}^{infty} f(t) exp(-i w t) dt
class TD_FD
{
	friend istream& operator>>(istream& in, TD_FD& dat);
	friend ostream& operator<<(ostream& out, const TD_FD& dat);
public:
	// inputs
	double timeScale;

	bool hasComputedFT;
	TD_signalT signalType;
	vector<double> modelParas;
	vector<int> modelFlags;
	// some cases involve dirac function in frequency domain. For that the range of regulalized diract function is [-wdirac, wdirac]
	// wdirac = relw4Dirac_w * delw
	double relw4Dirac_w;
	double tmin;
	double tmax;
	double delt;
	double wmin;
	double wmax;
	double delw;
	bool uniform_ts;

	TD_FD();
	// values above are provided. 
	// if the type is general, 	vector<double> td_vals must have been set
	void Set_ts_ws_external(vector<double>& modelParasIn, vector<double>& tsIn, int& uniform_ts_i, double& wminIn, double& wmaxIn, double& delwIn, bool initialize_w);
	void Initialize_AfterRead();
	void Initialize_TD_FD(bool initialize_FD);

	void Compute_TD_Values();
	void Compute_FD_Values();

	// computed
	// times
	int numt;
	int numw;
	vector<double> ts;
	// frequencies
	vector<double> ws;
	vector<double> td_vals;
	vector<double> td_DTs;

	vector<Dcomplex> fd_vals;

	////////////////////////////////////////

	bool b_modelParas;
	bool b_tmin;
	bool b_tmax;
	bool b_delt;
	int n_t_lims;
	bool b_wmin;
	bool b_wmax;
	bool b_delw;
	int n_w_lims;

	// return is val
	double ComputeSpecific_TD_Normalized_Value(double t, double& dValdT, bool compute_dValdT);
	double ComputeSpecific_TD_Value(double t);

	Dcomplex ComputeSpecific_FD_Value(double w);
	// TD values, ... are sized based on mint, maxt, delt
	void Size_TD();

	//// used to compute input energy through incident wave, ... and average stress
	void Compute_Inegral_load_load2(double timeMax, double delT, double& integral_loadOut, double& integral_load2Out);
};

// return value fhat_w
Dcomplex ComputeFourierTransform_UniformTime(double w, const vector<double>& ts, vector<double>& ft, bool uniform_ts);
void ComputeFourierTransformValues_UniformTime(const vector<double>& ws, const vector<double>& ts, vector<double>& ft, vector<Dcomplex>& fd_vals, bool uniform_ts);

string getName(TD_signalT dat);
bool fromName(const string& name, TD_signalT& dat);
ostream& operator<<(ostream& out, TD_signalT dat);
istream& operator>>(istream& in, TD_signalT& dat);

void Test_TD_FD_Class(TD_signalT signalType = td_sine_TD_Gaussian);

#endif