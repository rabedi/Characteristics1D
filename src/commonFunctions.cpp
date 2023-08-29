#include "commonFunctions.h"
#include <float.h>
#include "commonMacros.h"
#include "globalFunctions.h"
#include <algorithm>
using namespace std;

int RUN_GVN = 0;
int RUN_GSN = 0;
int g_i_serial = -1;
int g_existenceMode = 2;

// less than zero means it's inactive
void gChangeSerialNum(string& filename)
{
	AddChangeShortVSerialNumber2FileName(filename, g_i_serial, g_existenceMode);
}

char* ITOA( int value, char* result, int base )
{
	// check that the base if valid
	
	if (base < 2 || base > 16) { *result = 0; return result; }
	
    
	
	char* out = result;
	
	int quotient = value;
	
    
	
	do {
        
		*out = "0123456789abcdef"[ std::abs( quotient % base ) ];
        
		++out;
        
		quotient /= base;
        
	} while ( quotient );
	
    
	
	// Only apply negative sign for base 10
	
	if ( value < 0 && base == 10) *out++ = '-';
	
    
	
	std::reverse( result, out );
	
	*out = 0;
	
	return result;
	
}

int getNewComponentPosition(vector<double>& vec, double newComp, double tol, bool& foundInVector)
{
	int i, sz = vec.size();

	for (i = 0; i < sz; ++i)
		if (DoublesAreEqual(vec[i], newComp, tol) == true)
			break;
	if (i < sz)
	{
		foundInVector = true;
		return i;
	}
	for (i = 0; i < sz; ++i)
		if (vec[i] > newComp)
			break;
	foundInVector = false;
	return i;
}


void addNewComponentPosition(vector<double>& vec, double newComp, double tol, int& pos, bool& foundInVector)
{
	if (pos == -1) // it means	
		pos = getNewComponentPosition(vec, newComp, tol, foundInVector);
	if (foundInVector == false)
		addNewComponentPosition(vec, newComp, pos);
}

void SplitFileName(const string& filename, string& path, string& name, string& ext, string& folderSpliter)
{
    size_t lastdot = filename.find_last_of(".");
	string remainder;
    if (lastdot == std::string::npos)
	{
		ext = "";
		remainder = filename;
	}
 	else
	{
		remainder = filename.substr(0, lastdot);
		ext = filename.substr(lastdot + 1);
	}
	folderSpliter = "/";
    size_t lastdota = remainder.find_last_of("/");
    size_t lastdotb = remainder.find_last_of("\\");
	lastdot = MIN(lastdota, lastdotb);

    if (lastdot == std::string::npos)
	{
		path = "";
		name = remainder;
	}
 	else
	{
		if (lastdotb < lastdota)
			folderSpliter = "\\";
		path = remainder.substr(0, lastdot);
		name = remainder.substr(lastdot + 1);
	}
}

void MergeFileName(string& filename, const string& path, const string& name, const string& ext, string folderSpliter)
{
	size_t sz = path.length();
	if (sz > 0)
	{
		filename = path;
		if (path[sz - 1] != folderSpliter[0]) // path does not have folder spliter
			filename += folderSpliter;
	}
	else
		filename = string("");

	filename += name;

	sz = ext.length();
	if (sz > 0)
	{
		if (ext[0] != '.') 
			filename += string(".");
		filename += ext;
	}
}

// existenceMode = 0, does not check existence of file
//				 = 1, checks if file exists, if not, exits the code
//				 = 2, checks if file exists, if not, uses original file name
void AddChangeSerialNumber2FileName(const string& filename, string& filenameNew, int serialNum, int existenceMode)
{
	string path, name, ext, folderSpliter;
	SplitFileName(filename, path, name, ext, folderSpliter);
	unsigned int sz = name.size();
	int indexEn = sz - 1, indexSt;
	for (indexSt = indexEn; indexSt >= 0; --indexSt)
	{
		char c = name[indexSt];
		bool isNum = ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') ||
			(c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
		if (!isNum)
			break;
	}
	//	char buf[255];
	//	name.copy(buf, indexSt, 0);
	//	string filenamePartNew = buf;
	string filenamePartNew = "";
	for (int i = 0; i <= indexSt; ++i)
		filenamePartNew += name[i];
	int szNum = indexEn - indexSt;
	string serial;
	toString(serialNum, serial);
	int szSerial = serial.size();
	int zeroPadSize = MAX(szNum - szSerial, 0);
	for (int i = 0; i < zeroPadSize; ++i)
		serial = "0" + serial;
	filenamePartNew += serial;
	MergeFileName(filenameNew, path, filenamePartNew, ext, folderSpliter);
	if (existenceMode > 0)
	{
		fstream in(filenameNew.c_str(), ios::in);
		if (!in.is_open())
		{
			if (existenceMode == 1)
			{
				cout << "filename\t" << filename << '\n';
				cout << "filenameNew\t" << filenameNew << '\n';
				THROW("Cannot open filenameNew made from filename\n");
			}
			else
				filenameNew = filename;
		}
	}
}

void AddChangeShortVSerialNumber2FileName(string& filename, int serialNum, int existenceMode)
{
	string filenameBK = filename;
	AddChangeSerialNumber2FileName(filenameBK, filename, serialNum, existenceMode);
}

void ChangeFileNameExtensionPath(string& newfilename, const string& filename, string newExt, string newPath, string newName, bool newPathAddedAtTheEndOfCurrenctPath, 
		string subtractEndOfFileName, string addEndOfFileName, string subtractBeginningOfFileName, string addBeginningFileName)
{
	string path, name, ext, folderSpliter;
	SplitFileName(filename, path, name, ext, folderSpliter);
	if (newExt.length() == 0)
		newExt = ext;
	else if ((newExt == "none") || (newExt == "None"))
		newExt = "";
	if (newPath.length() == 0)
		newPath = path;
	else if ((newPath == "none") || (newPath == "None"))
		newPath = "";
	else if (newPathAddedAtTheEndOfCurrenctPath == true)
		newPath = path + folderSpliter + newPath;
	if (newName.length() == 0)
		newName = name;

	// modifying the base file name
	// end part:
	size_t sz = subtractEndOfFileName.length();
	if ((sz > 0) && (newName.compare(newName.size()-sz, sz, subtractEndOfFileName) == 0))
		newName = newName.substr(0, newName.size()-sz);

	newName = newName + addEndOfFileName;

	// beginning of the file
	sz = subtractBeginningOfFileName.length();
	if ((sz > 0) && (newName.compare(0, sz, subtractBeginningOfFileName) == 0))
		newName = newName.substr(sz, newName.size());

	newName = addBeginningFileName + newName;
	MergeFileName(newfilename, newPath, newName, newExt, folderSpliter);
}

bool AddSerialNumber2FileName(const string & fileNameWOSerial, string & fileNameOut, int serialNo, string txtBeforeSerial, bool checkIfFileCanBeOpened)
{
	fileNameOut = fileNameWOSerial;
	if (serialNo < 0)
	{
		if (checkIfFileCanBeOpened)
		{
			fstream in(fileNameOut.c_str(), ios::in);
			if (!in.is_open())
			{
				cout << "fileNameOut\t" << fileNameOut << " cannot be opened\n";
				THROW("File cannot be opened\n");
				//				return false;
			}
		}
		return true;
	}
	string path, name, ext, spliter;
	string meshNameWOExt = removeExtension(fileNameWOSerial);
	SplitFileName(fileNameWOSerial, path, name, ext, spliter);
	string dotext = "";
	if (ext.size() > 0)
		dotext = "." + ext;

	string ser;
	toString(serialNo, ser);
	meshNameWOExt = meshNameWOExt + txtBeforeSerial + ser;
	fileNameOut = meshNameWOExt + dotext;
	fstream in(fileNameOut.c_str(), ios::in);
	if (!in.is_open())
	{
		if (serialNo == 0)
		{
			fileNameOut = fileNameWOSerial;
			in.close();
			in.open(fileNameOut.c_str(), ios::in);
		}
		if (checkIfFileCanBeOpened && (!in.is_open()))
		{
			cout << "fileNameOut\t" << fileNameOut << " cannot be opened\n";
			THROW("File cannot be opened\n");
			//				return false;
		}
	}
	return true;
}

void print_state (const std::ios& stream)
{
  std::cout << " good()=" << stream.good();
  std::cout << " eof()=" << stream.eof();
  std::cout << " fail()=" << stream.fail();
  std::cout << " bad()=" << stream.bad() << '\n';
}

bool AreSame(const double a, const double b, double tolerance)
{
	double TOL;
	TOL = (tolerance < 0.0) ? DBL_EPSILON : tolerance;

	return fabs(a - b) < TOL;
}

bool LessThan(const double a, const double b, double tolerance)
{
	double TOL;
	TOL = (tolerance < 0.0) ? DBL_EPSILON : tolerance;

	return (b - a) >= TOL;
}

bool LessThanOrEqual(const double a, const double b, double tolerance)
{
	return LessThan(a, b, tolerance) || AreSame(a, b, tolerance);
}

bool GreaterThan(const double a, const double b, double tolerance)
{
	double TOL;
	TOL = (tolerance < 0.0) ? DBL_EPSILON : tolerance;

	return (a - b) >= TOL;
}

bool GreaterThanOrEqual(const double a, const double b, double tolerance)
{
	return GreaterThan(a, b, tolerance) || AreSame(a, b, tolerance);
}

bool AreSame(const vector<double>& a, const vector<double>& b, double tolerance)
{
	double TOL;
	TOL = (tolerance < 0.0) ? DBL_EPSILON : tolerance;

	bool same = true;

	if ((a.size() != b.size()) || a.empty())
	{
		return false;
	}

	for (int i = 0; i < (int)a.size(); ++i)
		same = same && AreSame(a[i], b[i], TOL);
	//same = same && (fabs(a[i] - b[i]) < TOL);

	return same;
}

double nFactorial(int n)
{
	return (n == 1 || n == 0) ? 1.0 : nFactorial(n - 1) * n;
}

void multivariatePolynomialTerms(vector<double>& crd, int level, vector<int>& I, vector<int>& Imax, vector<double>& terms, int& arbitrary_counter, bool computeDerivative, int wrtJ, int nthDiff)
{
	int dim = I.size();

	if (level > 0)
	{
		for (I[level - 1] = 0; I[level - 1] <= Imax[level - 1]; ++I[level - 1])
		{
			multivariatePolynomialTerms(crd, level - 1, I, Imax, terms, arbitrary_counter, computeDerivative, wrtJ, nthDiff);
		}
	}
	else
	{

		double tmpTerm = 1.0;
		for (int d = 0; d < dim; ++d)
		{
			//int row = arbitrary_counter;
			//int col = d;
			if (!computeDerivative)
			{
				tmpTerm *= pow(crd[d], I[d]);
			}
			else
			{
				if (wrtJ != d)
				{
					tmpTerm *= pow(crd[d], I[d]);
				}
				else
				{
					int effPow = I[d] - nthDiff;
					double effCoeff = (effPow < 0) ? 0.0 : nFactorial(I[d]) / nFactorial(effPow);
					effPow = (effPow < 0) ? 1 : effPow;
					tmpTerm *= (effCoeff * pow(crd[d], effPow));
				}
			}
		}
		terms.push_back(tmpTerm);

		++arbitrary_counter;

		bool end_of_loop = true;
		for (int e = 0; e < dim; ++e)
		{
			end_of_loop = end_of_loop && (I[e] == Imax[e]);
		}

		if (end_of_loop)
		{
			int fullDof = 1;
			for (int j = 0; j < dim; ++j)
			{
				fullDof *= Imax[j] + 1;
			}

			if (fullDof != (arbitrary_counter))
				THROW("number of matrix rows must agree with number of total polynomial degrees of freedom.");
		}
	}
}

double PolyValue(double x, vector<double>& coef, int p)
{
	if (p < 0)
		p = coef.size() - 1;
	double val = 0.0;
	double xPow = 1.0;
	for (int i = 0; i <= p; ++i)
	{
		val += xPow * coef[i];
		xPow *= x;
	}
	return val;
}

double genPolynomialFunc(vector<double>& crd, vector<int> polyOrder, vector<double>& coefficients, bool computeDerivative, int wrtJ, int nthDiff)
{
	int numVar = crd.size();
	if (numVar != polyOrder.size())
		THROW("number of polynomial orders specified must be same length as coordinate components.");

	int numCoeffs = coefficients.size();

	int spaceDim = numVar;
	int counter = 0;

	vector<int>I(numVar);
	vector<double> Xterms;
	multivariatePolynomialTerms(crd, numVar, I, polyOrder, Xterms, counter, computeDerivative, wrtJ, nthDiff);

	int numTerms = Xterms.size();
	if (numTerms != numCoeffs)
	{
		stringstream ss;
		ss << "Insufficient number of coefficients (" << numCoeffs << ")--polynomial needs (" << numTerms << ") coefficients.";
	}

	double valOut = 0.0;
	for (int i = 0; i < numTerms; ++i)
	{
		valOut += (coefficients[i] * Xterms[i]);
	}
	return valOut;
}

double genPolynomialDerivative(vector<double>& crd, vector<int> polyOrder, vector<double>& coefficients, int j, int n)
{
	return genPolynomialFunc(crd, polyOrder, coefficients, true, j, n);
}

double genUnivariateExponentialFunc(vector<double>& crd, int numTerms, vector<double>& coefficients)
{
	int NUMCOEFFSPERTERM = 3;
	int numCoeff = coefficients.size();
	if (crd.size() > 1)
	{
		THROW("NO IMPLEMENTATION FOR MULTIVARIATE CASE.");
	}
	if ((NUMCOEFFSPERTERM * numTerms) != numCoeff)
		THROW("Insufficient number of coefficients.");

	double valOut = 0.0;
	for (int i = 0; i < numCoeff; i += NUMCOEFFSPERTERM)
	{
		double a = coefficients[i];
		double b = coefficients[i + 1];
		double c = coefficients[i + 2];
		valOut += (a * exp(b * crd[0] + c));
	}
	return valOut;
}

double genUnivariateExponentialDerivative(vector<double>& crd, int numTerms, vector<double>& coefficients, int wrtJ, int nthDiff)
{
	int NUMCOEFFSPERTERM = 3;
	int numCoeff = coefficients.size();
	if (crd.size() > 1)
	{
		THROW("NO IMPLEMENTATION FOR MULTIVARIATE CASE.");
	}
	if ((NUMCOEFFSPERTERM * numTerms) != numCoeff)
		THROW("Insufficient number of coefficients.");

	double valOut = 0.0;
	for (int i = 0; i < numCoeff; i += NUMCOEFFSPERTERM)
	{
		double a = coefficients[i];
		double b = coefficients[i + 1];
		double c = coefficients[i + 2];
		valOut += (a * pow(b, nthDiff) * exp(b * crd[0] + c));
	}
	return valOut;
}

double genUnivariateHarmonicFunc(vector<double>& crd, int numTerms, vector<double>& coefficients)
{
	int NUMCOEFFSPERTERM = 3;
	int numCoeff = coefficients.size();
	if (crd.size() > 1)
	{
		THROW("NO IMPLEMENTATION FOR MULTIVARIATE CASE.");
	}
	if ((NUMCOEFFSPERTERM * numTerms) != numCoeff)
		THROW("Insufficient number of coefficients.");

	double valOut = 0.0;
	for (int i = 0; i < numCoeff; i += NUMCOEFFSPERTERM)
	{
		double a = coefficients[i];
		double b = coefficients[i + 1];
		double c = coefficients[i + 2];
		valOut += (a * sin(b * crd[0] + c));
	}
	return valOut;
}

double genUnivariateHarmonicDerivative(vector<double>& crd, int numTerms, vector<double>& coefficients, int wrtJ, int nthDiff)
{
	int NUMCOEFFSPERTERM = 3;
	int numCoeff = coefficients.size();
	if (crd.size() > 1)
	{
		THROW("NO IMPLEMENTATION FOR MULTIVARIATE CASE.");
	}
	if ((NUMCOEFFSPERTERM * numTerms) != numCoeff)
		THROW("Insufficient number of coefficients.");

	double valOut = 0.0;
	for (int i = 0; i < numCoeff; i += NUMCOEFFSPERTERM)
	{
		double a = coefficients[i];
		double b = coefficients[i + 1];
		double c = coefficients[i + 2];
		if (((nthDiff - 0) % 4) == 0)
			valOut += (a * pow(b, nthDiff) * 1.0 * sin(b * crd[0] + c));
		else if (((nthDiff - 1) % 4) == 0)
			valOut += (a * pow(b, nthDiff) * 1.0 * cos(b * crd[0] + c));
		else if (((nthDiff - 2) % 4) == 0)
			valOut += (a * pow(b, nthDiff) * -1.0 * sin(b * crd[0] + c));
		else if (((nthDiff - 3) % 4) == 0)
			valOut += (a * pow(b, nthDiff) * -1.0 * cos(b * crd[0] + c));
	}
	return valOut;
}

double genUnivariateLogFunc(vector<double>& crd, int numTerms, vector<double>& coefficients)
{
	int NUMCOEFFSPERTERM = 3;
	int numCoeff = coefficients.size();
	if (crd.size() > 1)
	{
		THROW("NO IMPLEMENTATION FOR MULTIVARIATE CASE.");
	}
	if ((NUMCOEFFSPERTERM * numTerms) != numCoeff)
		THROW("Insufficient number of coefficients.");

	double valOut = 0.0;
	for (int i = 0; i < numCoeff; i += NUMCOEFFSPERTERM)
	{
		double a = coefficients[i];
		double b = coefficients[i + 1];
		double c = coefficients[i + 2];
		valOut += (a * log(b * crd[0] + c));
	}
	return valOut;
}

double genUnivariateLogDerivative(vector<double>& crd, int numTerms, vector<double>& coefficients, int wrtJ, int nthDiff)
{
	int NUMCOEFFSPERTERM = 3;
	int numCoeff = coefficients.size();
	if (crd.size() > 1)
	{
		THROW("NO IMPLEMENTATION FOR MULTIVARIATE CASE.");
	}
	if ((NUMCOEFFSPERTERM * numTerms) != numCoeff)
		THROW("Insufficient number of coefficients.");

	double valOut = 0.0;
	for (int i = 0; i < numCoeff; i += NUMCOEFFSPERTERM)
	{
		double a = coefficients[i];
		double b = coefficients[i + 1];
		double c = coefficients[i + 2];
		valOut += ((pow(-1.0, (nthDiff - 1))*nFactorial(nthDiff - 1)*pow(b, nthDiff)) / (pow((b*crd[0] + c), nthDiff)));
	}
	return valOut;
}

int Sign(double x)
{
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

bool IsRelativeErrorAccepable(double valA, double valB, double relTol, double& relError)
{
	double denom = MAX(abs(valA), abs(valB));
	double er = abs(valA - valB);
	double absDenom;
	double r = denom / relTol;
	if (r < 1)
		absDenom = (1.0 - r * (1.0 - denom));
	else
		absDenom = denom;
	relError = er / absDenom;
	return (relError < relTol);
}

void ReadSetInteger(istream& in, std::set<int>& dat)
{
	int tmp;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (fromString(buf, tmp) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("Value must be integer\n");
		}
		dat.insert(tmp);
		READ_NSTRING(in, buf, buf);
	}
}

void ReadVectorInteger(istream& in, vector<int>& dat)
{
	int tmp;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (fromString(buf, tmp) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("Value must be integer\n");
		}
		dat.push_back(tmp);
		READ_NSTRING(in, buf, buf);
	}
}

void ReadVectorString(istream& in, vector<string>& dat)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		dat.push_back(buf);
		READ_NSTRING(in, buf, buf);
	}
}

void ReadMapInteger2Integer(istream& in, std::map<int, int>& dat)
{
	dat.clear();
	int key;
	int val;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf != "(")
			THROW("Entry for map should be enclosed in (\n");
		READ_NSTRING(in, buf, buf);
		if (fromString(buf, key) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("key must be integer\n");
		}
		READ_NSTRING(in, buf, buf);
		if (buf != ",")
			THROW("Entries must be separated by , (\n");
		READ_NSTRING(in, buf, buf);
		if (fromString(buf, val) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("val must be integer\n");
		}
		READ_NSTRING(in, buf, buf);
		if (buf != ")")
			THROW("Entry for map should be enclosed in )\n");
		dat[key] = val;
		READ_NSTRING(in, buf, buf);
	}
}

void ReadMapString2String(istream& in, std::map<string, string>& dat)
{
	dat.clear();
	string key;
	string val;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf != "(")
		{
			cout << "buf\t" << buf << '\n';
			THROW("Entry for map should be enclosed in (\n");
		}
		READ_NSTRING(in, buf, buf);
		key = buf;
		READ_NSTRING(in, buf, buf);
		if (buf != ",")
		{
			cout << "buf\t" << buf << '\n';
			THROW("Entries must be separated by , (\n");
		}
		READ_NSTRING(in, buf, buf);
		val = buf;
		READ_NSTRING(in, buf, buf);
		if (buf != ")")
		{
			cout << "buf\t" << buf << '\n';
			THROW("Entry for map should be enclosed in )\n");
		}
		dat[key] = val;
		READ_NSTRING(in, buf, buf);
	}
}

void ReadMapInteger2Double(istream& in, std::map<int, double>& dat)
{
	dat.clear();
	int key;
	double val;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf != "(")
			THROW("Entry for map should be enclosed in (\n");
		READ_NSTRING(in, buf, buf);
		if (fromString(buf, key) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("key must be integer\n");
		}
		READ_NSTRING(in, buf, buf);
		if (buf != ",")
			THROW("Entries must be separated by , (\n");
		READ_NSTRING(in, buf, buf);
		if (fromString(buf, val) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("val must be double\n");
		}
		READ_NSTRING(in, buf, buf);
		if (buf != ")")
			THROW("Entry for map should be enclosed in )\n");
		dat[key] = val;
		READ_NSTRING(in, buf, buf);
	}
}

void ReadSetDouble(istream& in, std::set<double>& dat)
{
	double tmp;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (fromString(buf, tmp) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("Value must be double\n");
		}
		dat.insert(tmp);
		READ_NSTRING(in, buf, buf);
	}
}

void ReadVectorDouble(istream& in, vector<double>& dat)
{
	double tmp;
	string buf;
	READ_NSTRING(in, buf, buf);
	dat.clear();
	if (buf != "{")
		THROW("istream should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (fromString(buf, tmp) == false)
		{
			cout << "buf\t" << buf << '\n';
			THROW("Value must be double\n");
		}
		dat.push_back(tmp);
		READ_NSTRING(in, buf, buf);
	}
}

void WriteSetInteger(ostream& out, const std::set<int>& dat)
{
	out << "{";
	set<int>::const_iterator it;
	for (it = dat.begin(); it != dat.end(); ++it)
		out << '\t' << *it;
	out << "\t}";
}

void WriteVectorInteger(ostream& out, const vector<int>& dat)
{
	out << "{";
	for (unsigned int i = 0; i < dat.size(); ++i)
		out << '\t' << dat[i];
	out << "\t}";
}

void WriteVectorString(ostream& out, const vector<string>& dat)
{
	out << "{";
	for (unsigned int i = 0; i < dat.size(); ++i)
		out << '\t' << dat[i];
	out << "\t}";
}

void WriteMapInteger2Integer(ostream& out, const std::map<int, int>& dat)
{
	out << "{";
	map<int, int>::const_iterator it;
	for (it = dat.begin(); it != dat.end(); ++it)
		out << "\t( " << it->first << " , " << it->second << " )";
	out << "\t}";
}

void WriteSetDouble(ostream& out, const std::set<double>& dat)
{
	out << "{";
	set<double>::const_iterator it;
	for (it = dat.begin(); it != dat.end(); ++it)
		out << '\t' << *it;
	out << "\t}";
}

void WriteVectorDouble(ostream& out, const vector<double>& dat)
{
	out << "{";
	for (unsigned int i = 0; i < dat.size(); ++i)
		out << '\t' << dat[i];
	out << "\t}";
}

void WriteMapInteger2Double(ostream& out, const std::map<int, double>& dat)
{
	out << "{";
	map<int, double>::const_iterator it;
	for (it = dat.begin(); it != dat.end(); ++it)
		out << "\t( " << it->first << " , " << it->second << " )";
	out << "\t}";
}

////////////////////////////////////// Infinitely smooth functions ///////////////////////////////////////////////////////////////
double getInfinitelySmoothTransitionFunctionX_0_2_Infinity_Y_0_2_1(double x, double& dydx, bool compute_dydx)
{
	dydx = 0.0;
	if (x < DBL_MIN)
		return 0.0;
	x = 1.0 / x;
	if (x > DBL_MAX)
		return 0.0;
	double val = exp(-x);
	if (compute_dydx == true)
	{
		dydx = x * x * val;
	}
	return val;
}

double getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1(double x, double& dydx, bool compute_dydx)
{
	dydx = 0.0;
	if (x <= 0.0)
		return 0.0;
	if (x >= 1.0)
		return 1.0;

	double dfx, df1mx;
	double fx = getInfinitelySmoothTransitionFunctionX_0_2_Infinity_Y_0_2_1(x, dfx, compute_dydx);
	double f1mx = getInfinitelySmoothTransitionFunctionX_0_2_Infinity_Y_0_2_1(1.0 - x, df1mx, compute_dydx);
	double denom = (fx + f1mx);
	if (compute_dydx == true)
		dydx = (dfx * f1mx + df1mx * fx) / denom / denom;
	return fx / denom;
}

double getInfinitelySmoothBumpFunctionX_0_2_1_Y_0_2_1_2_0(double x, double& dydx, bool compute_dydx)
{
	dydx = 0.0;
	if (x <= 0.0)
		return 0.0;
	if (x >= 1.0)
		return 0.0;

	double dfx, df1mx;
	double fx = getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1(x, dfx, compute_dydx);
	double f1mx = getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1(1.0 - x, df1mx, compute_dydx);
	if (compute_dydx == true)
		dydx = 4.0 *  (dfx * f1mx - df1mx * fx);
	return 4.0 * fx * f1mx;
}

double getInfinitlySmoothDeltaDirac(double x, double & dydx, double diracPositivieWidth, bool compute_dydx)
{
	static double factor4_integral1 = 2.6512611503134806412163015920669;
	//	double x_between0to1 = (x + diracPositivieWidth) / 2 / diracPositivieWidth;
	double x_between0to1 = 0.5 * (x / diracPositivieWidth + 1);
	double yFactor = 0.5 * factor4_integral1 / diracPositivieWidth;
	double dy_dxTmp;
	double y = getInfinitelySmoothBumpFunctionX_0_2_1_Y_0_2_1_2_0(x_between0to1, dy_dxTmp, compute_dydx);
	y *= yFactor;
	if (compute_dydx)
	{
		dydx = yFactor * 0.5 / diracPositivieWidth * dy_dxTmp;
	}
	return y;
}

double getInfinitlySmoothRampFunction(double x, double& dydx, double xLimMax, double functionValMax, double xLimMin, double functionValMin, bool compute_dydx)
{
	dydx = 0.0;
	if (x <= xLimMin)
		return functionValMin;
	if (x >= xLimMax)
		return functionValMax;
	double dgdx;
	double dely = (functionValMax - functionValMin);
	double delx = (xLimMax - xLimMin);
	double val = functionValMin + dely * getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1((x - xLimMin) / delx, dgdx, compute_dydx);
	if (compute_dydx == true)
		dydx = dely / delx * dgdx;
	return val;
}

double getInfinitelySmoothBumpFunction(double x, double& dydx, double xLimMin, double xLimMax, double bumpMaxValue, double bumpMinValue, bool compute_dydx)
{
	dydx = 0.0;
	double delx = (xLimMax - xLimMin);
	double dely = bumpMaxValue - bumpMinValue;
	double xg = (x - xLimMin) / delx;
	double dydxg;
	double g = getInfinitelySmoothBumpFunctionX_0_2_1_Y_0_2_1_2_0(xg, dydxg, compute_dydx);
	if (compute_dydx == true)
		dydx = dydxg * dely / delx;
	return bumpMinValue + g * dely;
}

double getInfinitelySmoothBumpWithPlateauFunction(double x, double& dydx, double a, double b, double c, double d, double bumpMaxValue, double bumpMinValue, bool compute_dydx)
{
	dydx = 0.0;
	double dely1 = (b - a);
	double dely2 = (d - c);
	double x1 = (x - a) / dely1;
	double x2 = (d - x) / dely2;
	double dy1, dy2;
	double delBumpVal = bumpMaxValue - bumpMinValue;
	double y1 = getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1(x1, dy1, compute_dydx);
	double y2 = getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1(x2, dy2, compute_dydx);
	if (compute_dydx == true)
		dydx = (dy1 * y2 / dely1 - dy2 * y1 / dely2) * delBumpVal;
	return bumpMinValue + y1 * y2 * delBumpVal;
}
