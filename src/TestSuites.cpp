#include "TestSuites.h"
#include "SLDescriptorData.h"

void Test_OnePointRiemannSolutions(string inputFileWOExt)
{
	string fileNameInput = inputFileWOExt + ".txt";
	fstream in(fileNameInput.c_str(), ios::in);
	SLInterfaceCalculator sifc;
	sifc.Read(in);
	string fileName = inputFileWOExt + "_scalars.txt";
	fstream outScalars(fileName, ios::out);
	fileName = inputFileWOExt + "_finalSolution.txt";
	fstream outFinalSln(fileName, ios::out);
	fileName = inputFileWOExt + "_adaptivity.txt";
	fstream outAdaptivity(fileName, ios::out);
	fileName = inputFileWOExt + "_iterationConv.txt";
	fstream outIterationConv(fileName, ios::out);
	bool accept_point;
	IOF_type iof_finalVals = iof_ascii, iof_scalarVals = iof_ascii;
	double space_or_time = 0.0;
	AdaptivityS as = sifc.Main_Compute_OnePoint(accept_point, iof_finalVals, iof_scalarVals, space_or_time, &outScalars, &outFinalSln, &outAdaptivity, &outIterationConv);
	cout << "accept_point\t" << accept_point << '\n';
	cout << "as\t" << as << '\n';
}

void Test_print_CyclicLoading()
{
	vector<double> paras(2);
	paras[0] = 5.0; // sigma
	paras[1] = 2.0; // period
	CyclicLoading cl;
	cl.Initialize_FromParamaters(paras);
	
	int numPer = 4;
	int numTPerStep = 100;
	double tMax = numPer * paras[1];
	double delT = paras[1] / numTPerStep;
	int sz = numPer * numTPerStep + 1;
//	vector<double> times(sz);
	fstream log("Test_print_CyclicLoading.txt", ios::out);
	for (int i = 0; i < sz; ++i)
	{
		double time = i * delT;
		double val = cl.getValue(time);
		log << time << '\t' << val << '\n';
	}
}
