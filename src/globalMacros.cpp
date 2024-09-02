#include "globalMacros.h"
#include "SLInterfaceFracturePF.h"
#include "DomainPostProcessingS2.h"

int g_interfaceFlags_IncImp[4];
int g_bulkFlags_IncImp[NUM_SIDES];

vector<string> pINDs;
vector<string> pINDcommas;
string g_prefileName = "";
string g_prefileNameWOSlash = "";
string g_prefileNamePPS2 = "";
double invalidNonnegativeNum = -1.0;
double invalidNum = 1e40;
int g_serialNumber = 0;
double g_time = 0.0;
int g_low_disk_space = 0;
int g_1_interface_low_disk_space = 0;
double gt0 = 0.0;

bool b_db_p = false;
fstream g_logout;

#if DB_ON
fstream dbout("_db.txt", ios::out);
#else
fstream dbout;
#endif

bool Is_Valid_Num(double val)
{
	return (val < 1e39);
}

void setGlobalMembers()
{
	if (g_slf_conf == NULL)
		g_slf_conf = new SLFractureGlobal_Configuration();

	//	printf("%d\n", _getmaxstdio());
//	_setmaxstdio(32768);
	//_setmaxstdio(8192);
//	printf("%d\n", _getmaxstdio());

	pINDs.resize(DiM);
	pINDcommas.resize(DiM);
#if DiM1
	pINDs[0] = "";
	pINDcommas[0] = "";
#else
	pINDs[0] = "1";
	pINDs[1] = "2";

	pINDcommas[0] = "1,";
	pINDcommas[1] = "2,";
#if DiM3
	pINDs[2] = "3";
	pINDcommas[2] = "3,";
#endif
#endif
	g_interfaceFlags_IncImp[LO_INT] = 2001;
	g_interfaceFlags_IncImp[LI_INT] = 2002;
	g_interfaceFlags_IncImp[RI_INT] = 2003;
	g_interfaceFlags_IncImp[RO_INT] = 2004;

	g_bulkFlags_IncImp[SDL] = 1001;
	g_bulkFlags_IncImp[SDR] = 1003;

	Contact_Damage_State_IO_Stat_1Field::SetStatics_Contact_Damage_State_IO_Stat_1Field();
}

void CopyFile2OutputDirectory(const string& baseFileName)
{
	string fileNameNew = getFileNameFromPath(baseFileName);
	fileNameNew = g_prefileName + "/" + fileNameNew;
	string source = baseFileName;
	fileOperation(copyF, source, fileNameNew);
}
 