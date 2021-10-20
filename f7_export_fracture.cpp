#include "stdifm.h"
#include "f7_export_fracture.h"
#include <fstream>  // to write in a text file
#include <iomanip>	// to format column width in text file

IfmModule g_pMod;  /* Global handle related to this plugin */

#pragma region IFM_Definitions
/* --- IFMREG_BEGIN --- */
/*  -- Do not edit! --  */

static IfmResult OnBeginDocument (IfmDocument);
static void OnEndDocument (IfmDocument);
static void Serialize (IfmDocument, IfmArchive);
static void OnEditDocument (IfmDocument, Widget);
static void OnLeaveSimulator (IfmDocument);
static void PostTimeStep (IfmDocument);

/*
 * Enter a short description between the quotation marks in the following lines:
 */
static const char szDesc[] = 
  "Please, insert a plug-in description here!";

#ifdef __cplusplus
extern "C"
#endif /* __cplusplus */

IfmResult RegisterModule(IfmModule pMod)
{
  if (IfmGetFeflowVersion (pMod) < IFM_REQUIRED_VERSION)
    return False;
  g_pMod = pMod;
  IfmRegisterModule (pMod, "SIMULATION", "F7_EXPORT_FRACTURE", "f7_export_fracture", 0x1000);
  IfmSetDescriptionString (pMod, szDesc);
  IfmSetCopyrightPath (pMod, "f7_export_fracture.txt");
  IfmSetHtmlPage (pMod, "f7_export_fracture.htm");
  IfmSetPrimarySource (pMod, "f7_export_fracture.cpp");
  IfmRegisterProc (pMod, "OnBeginDocument", 1, (IfmProc)OnBeginDocument);
  IfmRegisterProc (pMod, "OnEndDocument", 1, (IfmProc)OnEndDocument);
  IfmRegisterProc (pMod, "Serialize", 1, (IfmProc)Serialize);
  IfmRegisterProc (pMod, "OnEditDocument", 1, (IfmProc)OnEditDocument);
  IfmRegisterProc (pMod, "OnLeaveSimulator", 1, (IfmProc)OnLeaveSimulator);
  IfmRegisterProc (pMod, "PostTimeStep", 1, (IfmProc)PostTimeStep);
  return True;
}

static void Serialize (IfmDocument pDoc, IfmArchive pArc)
{
  CF7ExportFracture::FromHandle(pDoc)->Serialize (pDoc, pArc);
}
static void OnEditDocument (IfmDocument pDoc, Widget wParent)
{
  CF7ExportFracture::FromHandle(pDoc)->OnEditDocument (pDoc, wParent);
}
static void OnLeaveSimulator (IfmDocument pDoc)
{
  CF7ExportFracture::FromHandle(pDoc)->OnLeaveSimulator (pDoc);
}
static void PostTimeStep (IfmDocument pDoc)
{
  CF7ExportFracture::FromHandle(pDoc)->PostTimeStep (pDoc);
}

/* --- IFMREG_END --- */
#pragma endregion


static IfmResult OnBeginDocument (IfmDocument pDoc)
{
  if (IfmDocumentVersion (pDoc) < IFM_CURRENT_DOCUMENT_VERSION)
    return false;

  try {
    IfmDocumentSetUserData(pDoc, new CF7ExportFracture(pDoc));
  }
  catch (...) {
    return false;
  }

  return true;
}

static void OnEndDocument (IfmDocument pDoc)
{
  delete CF7ExportFracture::FromHandle(pDoc);
}

///////////////////////////////////////////////////////////////////////////
// Implementation of CF7ExportFracture

// Constructor
CF7ExportFracture::CF7ExportFracture (IfmDocument pDoc)
  : m_pDoc(pDoc)
{
  /*
   * TODO: Add your own code here ...
   */
}

// Destructor
CF7ExportFracture::~CF7ExportFracture ()
{
  /*
   * TODO: Add your own code here ...
   */
}

// Obtaining class instance from document handle
CF7ExportFracture* CF7ExportFracture::FromHandle (IfmDocument pDoc)
{
  return reinterpret_cast<CF7ExportFracture*>(IfmDocumentGetUserData(pDoc));
}

// Callbacks

char*	g1_filename;
double  g2_eff_timestep;			// effective time step, dissolution reaction will be extrapolated to this value in days
double  g3_stabilization;			// stabilization time after fracture is modified in days
double	g4_interval;				// interval at which results are exported (dissolution reaction time) in years
double  g5_accum_time;				// accumulates flow simulation time
int		g6_export_simple;			// number to make export files shorter

void CF7ExportFracture::Serialize (IfmDocument pDoc, IfmArchive pArc)
{
	switch (IfmioGetMode(pArc)) {				// decide what's going on: INIT, STORE, LOAD, FREE? 
	case IfmIO_STORE:							// activated when click saving
		IfmioString(pArc, &g1_filename, 100);
		IfmioDouble(pArc, &g2_eff_timestep);
		IfmioDouble(pArc, &g3_stabilization);
		IfmioDouble(pArc, &g4_interval);
		IfmioInt(pArc, &g6_export_simple);
		break;
	case IfmIO_LOAD:							// activated when file is open
		if (IfmioGetVersion(pArc) >= 0x1000) {
			// Serialize.
			IfmioString(pArc, &g1_filename, 100);
			IfmioDouble(pArc, &g2_eff_timestep);
			IfmioDouble(pArc, &g3_stabilization);
			IfmioDouble(pArc, &g4_interval);
			IfmioInt(pArc, &g6_export_simple);
		}
		break;
	case IfmIO_INIT:							// initialize variables
		break;
	case IfmIO_FREE:							// clean memory when file is closed 
		break;
	}
}

void CF7ExportFracture::OnEditDocument (IfmDocument pDoc, Widget wParent)
{
	char*	input1 = g1_filename;
	double	input2 = g2_eff_timestep;
	double	input3 = g3_stabilization;
	double	input4 = g4_interval;
	int	input5 = g6_export_simple;

	IfmProperty props[] = {
		{ "1 Export file", IfmPROP_STRING, &input1, &g1_filename, "export filename example C : \\ directory \\ file.txt" },
		{ "2 Effective dissolution time-step", IfmPROP_DOUBLE, &input2, &g2_eff_timestep, "Effective dissolution time-step in DAYS" },
		{ "3 Mass transport stabilization time", IfmPROP_DOUBLE, &input3, &g3_stabilization, "Mass transport stabilization time in DAYS" },
		{ "4 Export interval", IfmPROP_DOUBLE, &input4, &g4_interval, "Data export interval time in YEARS" },
		{ "5 Short export", IfmPROP_INT, &input5, &g6_export_simple, "Export every X fracture. Minimum value = 1" }
	};
	IfmEditProperties(pDoc, "export data", "My parameters: ", props, 5); // Last number is amount of input values.

	g1_filename = input1;
	g2_eff_timestep = input2;
	g3_stabilization = input3;
	g4_interval = input4;
	g6_export_simple = input5;
}

void CF7ExportFracture::OnLeaveSimulator (IfmDocument pDoc)
{
  /* 
   * TODO: Add your own code here ...
   */
}

void CF7ExportFracture::PostTimeStep (IfmDocument pDoc)
{
	double time_inc = IfmGetCurrentTimeIncrement(pDoc);
	g5_accum_time = g5_accum_time + time_inc;
	double export_interval_days = g4_interval * 365;
	double time_factor = g2_eff_timestep / g3_stabilization;

	if (g5_accum_time * time_factor > export_interval_days)	{

		using namespace std;


		double current_time = IfmGetAbsoluteSimulationTime(pDoc);
		double export_time = current_time * time_factor / 365;
		int fracture_count = IfmGetNumberOfTotalFractureElements(pDoc);

		const char* filepath = IfmGetProblemPath(pDoc);

		// current date/time based on current system
		time_t now = time(NULL);
		char dt[26] = {};
		ctime_s(dt, 26, &now);
		
		// Dynamic memory allocation (new command)) for pointer 'ptr' -> IfmFracNodes FeFlow C++ structure (int nbn;int nop[];) */
		IfmFracNodes *ptr = new IfmFracNodes;

		if (export_time < 1.5 * g4_interval)
		{
			ofstream ofs(g1_filename, ios::out);

			ofs << filepath << endl;
			ofs << "Simulation started on : " << dt << endl;
			
			ofs << std::setfill(' ') << std::setw(8) << "fracture" << " , ";
			ofs << std::setfill(' ') << std::setw(8) << "X1" << " , ";
			ofs << std::setfill(' ') << std::setw(8) << "Y1" << " , ";
			ofs << std::setfill(' ') << std::setw(8) << "X2" << " , ";
			ofs << std::setfill(' ') << std::setw(8) << "Y2" << endl;

			for (int n = 0; n < fracture_count; n = n + g6_export_simple)
			{
				IfmGetNodalArrayOfFractureElement(pDoc, n, ptr);	// Get number of nodes (nbn) and node mesh number (nop[]) and store it in structure IfmFracNodes pointed by 'ptr' */
				int nodes_per_frac = ptr->nbn;						// Get number of nodes per fracture from ptr.nbn and store it in 'nodes per frac' */
				double* c = new double[nodes_per_frac];				// Dynamic memory allocation for vector 'c' than can hold 'nodes_per_frac' number of int variables */

				double x1 = IfmGetX(pDoc, ptr->nop[0]);
				double y1 = IfmGetY(pDoc, ptr->nop[0]);
				double x2 = IfmGetX(pDoc, ptr->nop[1]);
				double y2 = IfmGetY(pDoc, ptr->nop[1]);

				ofs << std::setfill(' ') << std::setw(8) << n << " , ";
				ofs << std::setfill(' ') << std::setw(8) << x1 << " , ";
				ofs << std::setfill(' ') << std::setw(8) << y1 << " , ";
				ofs << std::setfill(' ') << std::setw(8) << x2 << " , ";
				ofs << std::setfill(' ') << std::setw(8) << y2 << endl;				
			}
			ofs << endl;							//separation line in text file
			ofs.close();
		}

		ofstream ofs(g1_filename, ios::app);

		ofs << "Data export : " << dt ;
		ofs << "time = " << export_time << " years " << endl;
		//		ofs << "fracture ,    X    , Y , feflow_cond ,  true_aper  ,   head  ,    area    ,    c_ceq    , flow_ls , Reynolds" << endl;
		
		ofs << std::setfill(' ') << std::setw(7)  << "n"		<< " , ";
		ofs << std::setfill(' ') << std::setw(10) << "true_aper"<< " , ";
		ofs << std::setfill(' ') << std::setw(7)  << "flow_eq"	<< " , ";
		ofs << std::setfill(' ') << std::setw(9)  << "frac_area"<< " , ";
		ofs << std::setfill(' ') << std::setw(10) << "flow_rate"<< " , ";
		ofs << std::setfill(' ') << std::setw(9)  << "conc"		<< " , ";
		ofs << std::setfill(' ') << std::setw(7)  << "head"		<< endl;

		for (int n = 0; n < fracture_count; n = n + g6_export_simple) {				//n +10 instead of n++ to reduce ouput file size

			IfmGetNodalArrayOfFractureElement(pDoc, n, ptr);	// Get number of nodes (nbn) and node mesh number (nop[]) and store it in structure IfmFracNodes pointed by 'ptr' */
			int nodes_per_frac = ptr->nbn;						// Get number of nodes per fracture from ptr.nbn and store it in 'nodes per frac' */
			double* c = new double[nodes_per_frac];				// Dynamic memory allocation for vector 'c' than can hold 'nodes_per_frac' number of int variables */

			/*MASS */

			for (int frac_node_index = 0; frac_node_index < nodes_per_frac; frac_node_index++)
			{
				c[frac_node_index] = IfmGetResultsTransportMassValue(pDoc, ptr->nop[frac_node_index]); // Get solute data from model using value from structure ptr.nop[] and store in vector c[]
			}

			double sum_solute = 0;

			for (int frac_node_index = 0; frac_node_index < nodes_per_frac; frac_node_index++) {	// START of subloop 2 level 1 trough array c[] of 'nodes_per_frac' size to accum solute values
				sum_solute = c[frac_node_index] + sum_solute;
			}	// END of subloop 2 level that accumulates solute values per fracture
			double ave_solute = sum_solute / nodes_per_frac;
			double conc = ave_solute;

			/* VELOCITY */

			double frac_flow_vel = IfmGetResultsVelocityNormValue(pDoc, ptr->nop[1]);				// fracture flow velocity in meters / day

			double feflow_frac_cross_A = IfmGetFracArea(pDoc, n, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);

			double flow_rate = feflow_frac_cross_A * frac_flow_vel;

			int frac_law = IfmGetFracLaw(pDoc, n, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES);	// Get fracture friction law, HP = 1, MS = 2

			double true_frac_aper = NULL;
			double true_frac_width = NULL;
			if (frac_law == 1) // HagenP
			{
			double feflow_frac_aper = IfmGetFracFlowConductivity(pDoc, n, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);

			double quadratic_IF = pow(feflow_frac_cross_A, 2) - 4 * feflow_frac_aper * feflow_frac_aper * feflow_frac_cross_A; // feflow_frac_aper = a , feflow_frac_cross_A = b
			
			if (quadratic_IF < 0)
			{
				true_frac_aper = feflow_frac_cross_A / (2 * feflow_frac_aper);
			}
			else
			{
				double up = feflow_frac_cross_A - pow(quadratic_IF, 0.5);
				double down = 2 * feflow_frac_aper;
				true_frac_aper = up / down;
			}
			
			true_frac_width = feflow_frac_cross_A / true_frac_aper;	// Estimate  fracture width = area / aperture in m
			}
			else if(frac_law == 2) // Manning-Strickler
			{
				double feflow_M_corr = IfmGetFracFlowConductivity(pDoc, n, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);
				double M_corr_ratio = feflow_M_corr / ( 50 *86400 ) ; //M_corr / M in seconds
				double C = pow(M_corr_ratio, 1.5) / 2;
				double square_root = pow(feflow_frac_cross_A - 4 * feflow_frac_cross_A*pow(C, 2), 0.5);
				true_frac_aper = (pow(feflow_frac_cross_A,0.5) - square_root) / (2 * C);
				true_frac_width = (pow(feflow_frac_cross_A,0.5) + square_root) / (2 * C);
			}

			double H_D = (2 * true_frac_aper * true_frac_width) / (true_frac_aper + true_frac_width);		// hydraulic diameter
			double Reynolds = (1000 * frac_flow_vel * H_D) / (86400 * 1.2e-3);	// reynolds number	

			double x = IfmGetX(pDoc, ptr->nop[1]);
			double y = IfmGetY(pDoc, ptr->nop[1]);


			double head(0.0);					// the misterious ifmbool *ine, see feflow message board
			IfmBool inDomain(false);			// http://forum.mikepoweredbydhi.com/index.php/topic,1094.msg2711.html#msg2711

			head = IfmGetResultsFlowHeadValueAtXYSlice(pDoc, x, y, 0, &inDomain);	// head at X,Y, slice number

			ofs << std::setfill(' ') << std::setw(7) << n << " , ";
//			ofs << std::setfill(' ') << std::setw(7) << x << " , ";
//			ofs << std::setfill(' ') << std::setw(1) << y << " , ";
//			ofs << std::setfill(' ') << std::setw(11) << feflow_cond << " , ";
			ofs << std::setfill(' ') << std::setw(10) << true_frac_aper << " , ";
			ofs << std::setfill(' ') << std::setw(7) << frac_law << " , ";
			ofs << std::setfill(' ') << std::setw(9) << feflow_frac_cross_A << " , ";
			ofs << std::setfill(' ') << std::setw(10) << flow_rate << " , ";
			ofs << std::setfill(' ') << std::setw(9) << conc << " , ";
			ofs << std::setfill(' ') << std::setw(7) << head << endl;
//			ofs << std::setfill(' ') << std::setw(7) << Reynolds << endl;

			delete[] c;

		}
		ofs << endl;							//separation line in text file

		ofs.close();

		g5_accum_time = 0 + (g5_accum_time * time_factor - g4_interval * 365);

		delete[] ptr;

	}
}

