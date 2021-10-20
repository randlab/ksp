#include "stdifm.h"
#include "2020_dreybrodt_F7.0.h"

#include <stdlib.h>						// library for atof = string to double 
#include <string.h>						// library for std::to_string = double to string
#include <sstream>
#include <math.h>						// Linux compilation needs this header

IfmModule g_pMod;  /* Global handle related to this plugin */

#pragma region IFM_Definitions
/* --- IFMREG_BEGIN --- */
/*  -- Do not edit! --  */

static IfmResult OnBeginDocument (IfmDocument);
static void OnEndDocument (IfmDocument);
static void Serialize (IfmDocument, IfmArchive);
static void OnEditDocument (IfmDocument, Widget);
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
  IfmRegisterModule (pMod, "SIMULATION", "2020_DREYBRODT_F7_0", "2020_dreybrodt_F7.0", 0x1000);
  IfmSetDescriptionString (pMod, szDesc);
  IfmSetCopyrightPath (pMod, "2020_dreybrodt_F7.0.txt");
  IfmSetHtmlPage (pMod, "2020_dreybrodt_F7.0.htm");
  IfmSetPrimarySource (pMod, "2020_dreybrodt_F7.0.cpp");
  IfmRegisterProc (pMod, "OnBeginDocument", 1, (IfmProc)OnBeginDocument);
  IfmRegisterProc (pMod, "OnEndDocument", 1, (IfmProc)OnEndDocument);
  IfmRegisterProc (pMod, "Serialize", 1, (IfmProc)Serialize);
  IfmRegisterProc (pMod, "OnEditDocument", 1, (IfmProc)OnEditDocument);
  IfmRegisterProc (pMod, "PostTimeStep", 1, (IfmProc)PostTimeStep);
  return True;
}

static void Serialize (IfmDocument pDoc, IfmArchive pArc)
{
  C2020DreybrodtF70::FromHandle(pDoc)->Serialize (pDoc, pArc);
}
static void OnEditDocument (IfmDocument pDoc, Widget wParent)
{
  C2020DreybrodtF70::FromHandle(pDoc)->OnEditDocument (pDoc, wParent);
}
static void PostTimeStep (IfmDocument pDoc)
{
  C2020DreybrodtF70::FromHandle(pDoc)->PostTimeStep (pDoc);
}

/* --- IFMREG_END --- */
#pragma endregion


static IfmResult OnBeginDocument (IfmDocument pDoc)
{
  if (IfmDocumentVersion (pDoc) < IFM_CURRENT_DOCUMENT_VERSION)
    return false;

  try {
    IfmDocumentSetUserData(pDoc, new C2020DreybrodtF70(pDoc));
  }
  catch (...) {
    return false;
  }

  return true;
}

static void OnEndDocument (IfmDocument pDoc)
{
  delete C2020DreybrodtF70::FromHandle(pDoc);
}

///////////////////////////////////////////////////////////////////////////
// Implementation of C2020DreybrodtF70

// Constructor
C2020DreybrodtF70::C2020DreybrodtF70 (IfmDocument pDoc)
  : m_pDoc(pDoc)
{
  /*
   * TODO: Add your own code here ...
   */
}

// Destructor
C2020DreybrodtF70::~C2020DreybrodtF70 ()
{
  /*
   * TODO: Add your own code here ...
   */
}

// Obtaining class instance from document handle
C2020DreybrodtF70* C2020DreybrodtF70::FromHandle (IfmDocument pDoc)
{
  return reinterpret_cast<C2020DreybrodtF70*>(IfmDocumentGetUserData(pDoc));
}

// Callbacks

double g1_CEQ_do;					// Calcite EQ concentration in Moles/liter in double for calculations
char*  g1_CEQ_st;					// Calcite EQ concentrarion in Moles/Liter in string for the OnEdit window

double g2_k_1_do;					// laminar flow: constant page 31 Dreybrodt book double value
char*  g2_k_1_st;					// laminar flow: constant page 31 Dreybrodt book string value

double g3_k_n_do;					// laminar flow: constant page 31 Dreybrodt book double value
char*  g3_k_n_st;					// laminar flow: constant page 31 Dreybrodt book string value

double g4_calcite_density;			// calcite density in moles per cubic meter
double g5_init_wait;				// initial time wait before the plug-in runs in the simulation 
int	   g6_messages;					// message level
double g7_eff_time_step;			// effective time step, dissolution reaction will be extrapolated to this value in days
double g8_pert_time;				// stabilization time after fracture is modified in days

double g9_dynamic_vis_do;			// dynamic viscosity in Pa*s = N *s /m^2
char*  g9_dynamic_vis_st;

double g10_diff_coeff_do;			// diffusion coeffcient dreybrodt page 31 
char*  g10_diff_coeff_st;

double g11_water_density_do;		// water density in kg/m3
char*  g11_water_density_st;

double g12_Manning_Coeff_do;		// manning coefficient (unitless)
char*  g12_Manning_Coeff_st;

int    g13_transition_int;			// switch to activate transition to turbulent

double g98_accum_update = 0;		// effective dissolution time in days
double g99_accum_pert = 0;			// perturbation time = 0

void C2020DreybrodtF70::Serialize (IfmDocument pDoc, IfmArchive pArc)
{
	switch (IfmioGetMode(pArc)) {				// decide what's going on: INIT, STORE, LOAD, FREE? 
	case IfmIO_STORE:							// activated when click saving
		IfmioString(pArc, &g1_CEQ_st, 10);
		IfmioString(pArc, &g2_k_1_st, 15);
		IfmioString(pArc, &g3_k_n_st, 15);
		IfmioDouble(pArc, &g4_calcite_density);
		IfmioDouble(pArc, &g5_init_wait);
		IfmioInt(pArc, &g6_messages);
		IfmioDouble(pArc, &g7_eff_time_step);
		IfmioDouble(pArc, &g8_pert_time);
		IfmioString(pArc, &g9_dynamic_vis_st, 15);
		IfmioString(pArc, &g10_diff_coeff_st, 15);
		IfmioString(pArc, &g11_water_density_st, 15);
		IfmioString(pArc, &g12_Manning_Coeff_st, 15);
		IfmioInt   (pArc, &g13_transition_int);
		break;
	case IfmIO_LOAD:							// activated when file is open
		if (IfmioGetVersion(pArc) >= 0x1000) {
			// Serialize.
			IfmioString(pArc, &g1_CEQ_st, 10);
			IfmioString(pArc, &g2_k_1_st, 15);
			IfmioString(pArc, &g3_k_n_st, 15);
			IfmioDouble(pArc, &g4_calcite_density);
			IfmioDouble(pArc, &g5_init_wait);
			IfmioInt(pArc, &g6_messages);
			IfmioDouble(pArc, &g7_eff_time_step);
			IfmioDouble(pArc, &g8_pert_time);
			IfmioString(pArc, &g9_dynamic_vis_st, 15);
			IfmioString(pArc, &g10_diff_coeff_st, 15);
			IfmioString(pArc, &g11_water_density_st, 15);
			IfmioString(pArc, &g12_Manning_Coeff_st, 15);
			IfmioInt(pArc, &g13_transition_int);
		}
		break;
	case IfmIO_INIT:							// initialize variables
		break;
	case IfmIO_FREE:							// clean memory when file is closed 
		break;
	}
}

void C2020DreybrodtF70::OnEditDocument (IfmDocument pDoc, Widget wParent)
{
	char *	input1 = g1_CEQ_st;
	char *	input2 = g2_k_1_st;
	char *	input3 = g3_k_n_st;
	double	input4 = g4_calcite_density;
	double	input5 = g5_init_wait;
	int		input6 = g6_messages;
	double	input7 = g7_eff_time_step;
	double	input8 = g8_pert_time;
	char*	input9 = g9_dynamic_vis_st;
	char*	input10 = g10_diff_coeff_st;
	char*   input11 = g11_water_density_st;
	char*   input12 = g12_Manning_Coeff_st;
	int		input13 = g13_transition_int;

	IfmProperty props[] = {
		{ "01 Calcite EQ", IfmPROP_STRING, &input1, NULL, "Input a double variable, the Calcite EQ concentration in Mols/L" },
		{ "02 k_1", IfmPROP_STRING, &input2, NULL, "k_1 constant in mol / (cm^2 * sec)" },
		{ "03 k_n", IfmPROP_STRING, &input3, NULL, "k_n constant in mol / (cm^2 * sec)" },
		{ "04 calcite density", IfmPROP_DOUBLE, &input4, NULL, "calcite density in Moles per cubic meter" },
		{ "05 Initial wait time", IfmPROP_DOUBLE, &input5, NULL, "initial wait time in days" },
		{ "06 Message level", IfmPROP_INT, &input6, NULL, "0 all messages, 1 just update, 2 no messages" },
		{ "07 Effective dissolution time step", IfmPROP_DOUBLE, &input7, NULL, "extrapolation time for dissolution reaction in days" },
		{ "08 Perturbation stabilization time", IfmPROP_DOUBLE, &input8, NULL, "days left for mass transport to stabilize after fracture changes" },
		{ "09 dynamic viscosity", IfmPROP_STRING, &input9, NULL, "dynamic viscosity in Pa*s = N *s /m^2, temperature dependant" },
		{ "10 diffusion coefficient", IfmPROP_STRING, &input10, NULL, "diffusion coefficient pag 31 dreybrodt in m^s / sec" },
		{ "11 water density", IfmPROP_STRING, &input11, NULL, "water density in kg/m3" },
		{ "12 Manning Strickler coefficient", IfmPROP_STRING, &input12, NULL, "Mannig-Strickler coefficient" },
		{ "13 Activate transition to turbulent?", IfmPROP_INT, &input13, NULL, "0 = NO, 1 = YES" },
	};

	IfmEditProperties(pDoc, "Constant inputs for Dreybrodt dissolution model", "My parameters: ", props, 13);

	g1_CEQ_st = input1;   // I save CEQ as string because double values only store 2 decimals. At the start of PostTimeStep I convert it to double.
	g2_k_1_st = input2;
	g3_k_n_st = input3;
	g4_calcite_density = input4;
	g5_init_wait = input5;
	g6_messages = input6;
	g7_eff_time_step = input7;
	g8_pert_time = input8;
	g9_dynamic_vis_st = input9;
	g10_diff_coeff_st = input10;
	g11_water_density_st = input11;
	g12_Manning_Coeff_st = input12;
	g13_transition_int = input13;
}

void C2020DreybrodtF70::PostTimeStep (IfmDocument pDoc)
{
	/* 0 Values with more than 2 decimals can't be stored on FEM file. Store them as strings and convert to double values to use them in calculations */

	g1_CEQ_do = atof(g1_CEQ_st);					// CEQ equilibrium concentration stored as string converted to double
	g2_k_1_do = atof(g2_k_1_st);					// k_1 string to double
	g3_k_n_do = atof(g3_k_n_st);					// k_n string to double	
	g9_dynamic_vis_do = atof(g9_dynamic_vis_st);	// dynamic viscosity string to double for Reynolds calculations
	g10_diff_coeff_do = atof(g10_diff_coeff_st);	// diffusion coefficient section 3.1.4.3
	g11_water_density_do = atof(g11_water_density_st);
	g12_Manning_Coeff_do = atof(g12_Manning_Coeff_st);

	
	/* 1 - Time controls and clock, PART 1 */
	double start_time = clock();

	double days = IfmGetAbsoluteSimulationTime(pDoc);			// This one makes the plug-in wait until initial wait time (g5) is accomplished
	if (days < g5_init_wait) { goto init_wait; }

	double delta_t;
	delta_t = IfmGetCurrentTimeIncrement(pDoc);		// This ones controls the interval in fracture size upgrades
	g99_accum_pert = g99_accum_pert + delta_t;		// could be done when transport has stabilized after fracture size is modified
	if (g99_accum_pert > g8_pert_time) { goto calculate_new_frac; }			// a if that reads fracture and changes source terms until new_source is close to last_source
	else { goto no_frac_update; }

	// ---------------------------------------

	calculate_new_frac:

	g98_accum_update = g98_accum_update + g7_eff_time_step;		// Accumulating effective dissolution time 

	/* 2 - Mesh variables */

	int ndim;
	ndim = IfmGetNumberOfDimensions(pDoc);
	int fractures;
	fractures = IfmGetNumberOfTotalFractureElements(pDoc);
	int warning_count = 0;		// warning count
	int W_C = 0;				//warning count

	/* 3 - Dinamic memory allocations for fracture update loop*/

	IfmFracNodes *ptr;
	ptr = new IfmFracNodes;						// Dynamic memory allocation (new command)) for pointer 'ptr' -> IfmFracNodes FeFlow C++ structure (int nbn;int nop[];) */
	double* coordi_1;
	coordi_1 = new double[ndim];				// Dynamic memory allocation 'coordi_1' to store coordinates, it holds 'ndim' number of doubles */
	double* coordi_2;
	coordi_2 = new double[ndim];				// Dynamic memory allocation 'coordi_2' to store coordinates, it holds 'ndim' number of doubles */


	/******************************************************************************************************************************/
	/* 4 START of big fracture update loop*/
	for (int frac_index = 0; frac_index < fractures; frac_index++)
	{
		IfmGetNodalArrayOfFractureElement(pDoc, frac_index, ptr);// Get number of nodes (nbn) and node mesh number (nop[]) and store it in structure IfmFracNodes pointed by 'ptr' */
		int nodes_per_frac = ptr->nbn;							 // Get number of nodes per fracture from ptr.nbn and store it in 'nodes per frac' */
		double* c = new double[nodes_per_frac];					 // Dynamic memory allocation for vector 'c' than can hold 'nodes_per_frac' number of int variables */

		/*4.0 Check for flow model : Hagen-Poiseuille or Manning-Strickler?*/

		int frac_law = IfmGetFracLaw(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES);	// Get fracture friction law, HP = 1, MS = 2

		/* 4.1 The distance between nodes, 2D & 3D, the 'x' in Dreybrodt 1990 */

		double frac_length = 0;

		if (ndim == 2)
		{
			coordi_1[0] = IfmGetX(pDoc, ptr->nop[0]);
			coordi_1[1] = IfmGetY(pDoc, ptr->nop[0]);
			coordi_2[0] = IfmGetX(pDoc, ptr->nop[1]);
			coordi_2[1] = IfmGetY(pDoc, ptr->nop[1]);
			frac_length = pow((pow(coordi_2[0] - coordi_1[0], 2) + pow(coordi_2[1] - coordi_1[1], 2)), 0.5);     // Fracure lenght 2D = Pythagoras ( X^2 + Y^2 ) ^ 0.5
		}
		else   // 3D problems
		{
			coordi_1[0] = IfmGetX(pDoc, ptr->nop[0]);
			coordi_1[1] = IfmGetY(pDoc, ptr->nop[0]);
			coordi_1[2] = IfmGetZ(pDoc, ptr->nop[0]);
			coordi_2[0] = IfmGetX(pDoc, ptr->nop[1]);
			coordi_2[1] = IfmGetY(pDoc, ptr->nop[1]);
			coordi_2[2] = IfmGetZ(pDoc, ptr->nop[1]);
			frac_length = pow((pow(coordi_2[0] - coordi_1[0], 2) + pow(coordi_2[1] - coordi_1[1], 2) + pow(coordi_2[2] - coordi_1[2], 2)), 0.5);     // Fracure lenght 2D = Pythagoras ( X^2 + Y^2 ) ^ 0.5
		}

		/* 4.2 Get fracture cross section dimensions */

		// Get fracture flow conductivity (aperture in M for Hagen Poiseuille, Roughness coefficient for Manning)
		// http://www.feflow.info/html/help71/feflow/13_Programming/IFM/API/IfmGetFracFlowConductivity.html


		double true_frac_aper = NULL;
		double true_frac_width = NULL;		
		double feflow_frac_cross_A = IfmGetFracArea(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);	// Get fracture cross section area in m^2

		if (frac_law == 1) // HagenP
		{
			double feflow_frac_aper = IfmGetFracFlowConductivity(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);

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
			double feflow_M_corr = IfmGetFracFlowConductivity(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);
			double M_corr_ratio = feflow_M_corr / ( g12_Manning_Coeff_do *86400 ) ; //M_corr / M in seconds
			double C = pow(M_corr_ratio, 1.5) / 2;
			double square_root = pow(feflow_frac_cross_A - 4 * feflow_frac_cross_A*pow(C, 2), 0.5);
			true_frac_aper = (pow(feflow_frac_cross_A,0.5) - square_root) / (2 * C);
			true_frac_width = (pow(feflow_frac_cross_A,0.5) + square_root) / (2 * C);
		}
		else
		{
			IfmInfo(pDoc, "Discrete features not Hagen-Poisuille or Manning-Strickler, simulation stopped");
			IfmSetSimulationControlFlag(pDoc, IfmCTL_BREAK);
		}
		
		/* 4.3.1 START of subloop 1 level 1. Get solute concentration values from fracture nodes, a fracture may have 2, 3 or 4 nodes */
		for (int frac_node_index = 0; frac_node_index < nodes_per_frac; frac_node_index++)
		{
			c[frac_node_index] = IfmGetResultsTransportMassValue(pDoc, ptr->nop[frac_node_index]); // Get solute data from model using value from structure ptr.nop[] and store in vector c[]
		} /* END of subloop 1 level 1 */

		/* 4.2.2 Estimate fracture average solute value from from fracture nodal values stored in array c[] of 'frac_node_index' size, see 4.1.1 */

		double sum_solute = 0;

		for (int frac_node_index = 0; frac_node_index < nodes_per_frac; frac_node_index++) {	// START of subloop 2 level 1 trough array c[] of 'nodes_per_frac' size to accum solute values
			sum_solute = c[frac_node_index] + sum_solute;
		}	// END of subloop 2 level that accumulates solute values per fracture

		double ave_solute = sum_solute / nodes_per_frac;									// Fracture averaged solute value 

		/***** condition of concentration*******/

		double head_0 = IfmGetResultsFlowHeadValue(pDoc, ptr->nop[0]);
		double head_1 = IfmGetResultsFlowHeadValue(pDoc, ptr->nop[1]);
		double ave_solute_2 = NULL;

		if      (c[0] <= c[1] && head_0 >= head_1) { ave_solute_2 = ave_solute; }
		else if (c[0] < c[1] && head_0 < head_1) { ave_solute_2 = c[1]; }
		else if (c[0] > c[1] && head_0 > head_1) { ave_solute_2 = c[0]; }
		else if (c[0] >= c[1] && head_0 <= head_1) { ave_solute_2 = ave_solute; } 
		else    ( IfmInfo(pDoc, "BIG ERROR lines 399-403 ave_solute_2, in fracture %i", frac_index) );
		
		//ave_solute_2 = ave_solute;  //////DEBUG!!!!!!!!!!!!!!!!

		// frac flow velocity for Hagen P
		//double cond_HP = (pow(true_frac_aper, 2) * 1000 * 9.81) / (12 * 0.0013); //
		//double gradient = (head_0 - head_1) / frac_length;
		//double frac_flow_vel = -1 * (cond_HP * (gradient) );  //true_frac_aper * true_frac_width *
			// frac flow velocity for Manning
			//double hyd_radius = (true_frac_aper * true_frac_width) / 2*(true_frac_aper + true_frac_width);		// hydraulic diameter (meters)
			//double cond_manning = (true_frac_aper * true_frac_width);
			//double flow_vel_manning = -1 * ( 50 * pow(hyd_radius, 0.6666666666)*pow(gradient,0.5) );

		//double frac_flow_vel_0 = IfmGetResultsVelocityNormValue(pDoc, ptr->nop[0]);
		//double frac_flow_vel_1 = IfmGetResultsVelocityNormValue(pDoc, ptr->nop[1]);

		//double ave_solute_2 = NULL;
		//
		//if (condition == 1)
		//{
		//	if (frac_flow_vel > 0) { ave_solute_2 = ave_solute;	}
		//	else { ave_solute_2 = c[1];	}
		//}
		//else
		//{
		//	if (frac_flow_vel > 0) { ave_solute_2 = c[0]; }
		//	else { ave_solute_2 = ave_solute; }
		//}
		/*****************************************************************/

		/* 4.1.5 Reynolds in fracture */
		//double frac_flow_vel_feflow = IfmGetResultsVelocityNormValue(pDoc, ptr->nop[1]);	// Ask feflow for velocity in node in meters per day (m/d)
		//double D_H = (2 * feflow_frac_cross_A) / (true_frac_aper + true_frac_width);		// hydraulic diameter (meters)
		//double Reynolds = (g11_water_density_do * frac_flow_vel_feflow * D_H) / (86400 * g9_dynamic_vis_do);// Re = (density = 1000 * vel * charac lenght = D_H ) / (dyn visc = 0.0012 N*s/m2 = kg*m*s/s2*m2 )

		/* 4.1.6 IF to choose rate depending on CEQ and laminar (1) or manning (2) */
		double dissolution_rate = NULL;

		if (ave_solute_2 < 0.9 * g1_CEQ_do)		/* 1 - 1st order dissolution kinetics when C < 0.9 CEQ  */
		{
												/* 2 - 4.1.6.1 1st oder dissolution rate when averaged solut  < 0.9 * CEQ under laminar flow conditions */
			if (frac_law == 1)					// Rate in mol/(m^2*s) for 1st order dissolution, Dreybrodt book pag 31 laminar flow
			{
				dissolution_rate = (g2_k_1_do / (1 + ((g2_k_1_do*true_frac_aper) / (3 * g10_diff_coeff_do*g1_CEQ_do)))*(1 - (ave_solute_2 / g1_CEQ_do))); 

			}
			else if (frac_law == 2)				/* 3 - 4.1.6.2  1st order dissolution rate when averaged solute < 0.9 CEQ under turbulent flow conditions*/
			{		 
				dissolution_rate = g2_k_1_do * (1 - (ave_solute_2 / g1_CEQ_do));
			}		
		}										//  end of ave_solute < 0.9 * CEQ
		else 
		{									
			if (ave_solute_2 > g1_CEQ_do)		// 7 - IN EQUILIBRIUM = no reaction because solute concentration is larger than equilibrium concentration
			{
				warning_count = warning_count + 1; W_C = warning_count;
			}		
			else								/* 8 - 4th order dissolution kinetics when 0.9 CEQ < averaged_solute < CEQ */
			{
				if (frac_law == 1 )				/* 9 - 4.1.6.3  4th order dissolution rate when 0.9 CEQ < average solute < CEQ under laminar flow conditions*/
				{								
					double fourth_order = g3_k_n_do * pow((1 - ave_solute_2 / g1_CEQ_do), 4);	// Rate in mol/(m^2*s) for 4th order dissolution kinetics, dreybrodt book page 31 laminar flow 
					
					// fix to discontinuity from 1st to 4th order rates
					double first_order = (g2_k_1_do / (1 + ((g2_k_1_do*true_frac_aper) / (3 * g10_diff_coeff_do*g1_CEQ_do)))*(1 - (ave_solute_2 / g1_CEQ_do))); 

					if (fourth_order >= first_order)
					{
						dissolution_rate = first_order;
						//IfmInfo(pDoc, "4th order rate > 1st order rate, in fracture %i, 10.03.2021 fix applied", frac_index);
					}
					else if (fourth_order < first_order) { dissolution_rate = fourth_order; }
					else { IfmInfo(pDoc, "Just here to catch IF ELSE fails in rate calculations, code line 475, fracture %i", frac_index); }					 
				}
				else if (frac_law == 2)			/* 10 - 4.1.6.4  4th order dissolution rate when 0.9 CEQ < averaged solute < CEQ under turbulent flow conditions*/
				{								
					dissolution_rate = g2_k_1_do * pow((1 - (ave_solute_2 / g1_CEQ_do)), 4);
				}

			}
		}

		/* 4.1.7 alternative retreat mixing corrosion*/
		//		double saline = (ave_solute - 174) / (34370 - 174);   // 1 simple mixing corrosion
		//		rate_copy = g2_k_1_do * (saline - pow(saline, 2));
		//		if (rate_copy < 0) { rate_copy = 0; }
		//		double wall_retreat = 2 * rate_copy * g7_eff_time_step * 86400 / g4_calcite_density;

		//double comp = IfmGetFracFlowCompressibility(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS) * 10000;  // NEW NEW NEW

		double wall_retreat = 2 * dissolution_rate * g7_eff_time_step * 86400  / g4_calcite_density;  // 86400 to go from mol/(m^2*s) to mol /(m^2*day) / 10'000 to go from mol/cm^2 to mol m^2
		//double wall_retreat = 2 * rate_copy * g7_eff_time_step * comp * 86400 / g4_calcite_density;  // 86400 to go from mol/(m^2*s) to mol /(m^2*day)

		/* 4.1.8 New fracture aperture and cross section area */
		double neu_true_frac_aper = true_frac_aper + wall_retreat;
		double neu_true_frac_width = true_frac_width + wall_retreat;       
		double neu_frac_cross_area = neu_true_frac_aper * neu_true_frac_width;

		/*DECIDE HERE NEW FRICTION LAW*/

		// Hagen P vel 
		double neu_cond_HP = (pow(neu_true_frac_aper, 2) * 1000 * 9.81) / (12 * 0.0013); //																  
		double abs_gradient = abs ( (head_0 - head_1)/frac_length );
		double abs_HP_vel = neu_cond_HP * abs_gradient;  //true_frac_aper * true_frac_width *	
		// frac flow velocity for Manning
		double neu_hyd_radius = (neu_true_frac_aper * neu_true_frac_width) / ( 2 * (neu_true_frac_aper + neu_true_frac_width) );		// hydraulic diameter (meters)
		double abs_manning_vel = (g12_Manning_Coeff_do * pow(neu_hyd_radius, 0.666666666666666 )*pow(abs_gradient, 0.5));
		if (g6_messages == 3)
		{
			IfmInfo(pDoc, "frac=%i coeff=%.2f hyd_rad=%.2e grad=%.2e manning=%.2e", frac_index,g12_Manning_Coeff_do,neu_hyd_radius,abs_gradient,abs_manning_vel);
		}

		// bottom part of laminar flow C < 0.9 CEQ equation pag 31 Dreybodt
		double neu_feflow_frac_aper = (neu_true_frac_aper * neu_true_frac_width) / (neu_true_frac_aper + neu_true_frac_width);

		double diff_limit = 1 + (g2_k_1_do * neu_true_frac_aper) / (3 * g10_diff_coeff_do * g1_CEQ_do);	 //	x100 m to cm, x1e.6 mmol/L to mol/cm^3 - 2021
		double av_ratio = (neu_true_frac_aper + neu_true_frac_width) / (500 * neu_frac_cross_area);  // area volume ratio
		double neu_decay_laminar = (av_ratio / diff_limit) * 86400;				// frac decay rate in 1/day to account for large fracture effect on disolution rate limited by diffusion */
		
		//new turbulent parmeters
		double neu_M_corr = pow((2 * neu_frac_cross_area) / ((neu_true_frac_aper + neu_true_frac_width) * pow(neu_frac_cross_area, 0.5)), 0.66666666666666666);//table 9.12 feflow white papers
		neu_M_corr = g12_Manning_Coeff_do * neu_M_corr * 86400;
		double new_decay_turbulent = av_ratio * 86400; //decay rate in turbulent = av_ratio
		
		/* 4.1.10 New fracture aperture and cross section area IF for laminar and turbulent*/

		double transition_trigger = 2;

		if (g13_transition_int == 0) // 0 == no transition, stay in laminar  
		{			
			IfmSetFracFlowConductivity(pDoc, frac_index, neu_feflow_frac_aper, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);	// Set new fracture hydraulic aperture 
			IfmSetFracArea(pDoc, frac_index, neu_frac_cross_area, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);				// Set new cross section area	
			IfmSetFracMassDecayRate(pDoc, frac_index, neu_decay_laminar, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);		// set new decay rate (account for Diff limit on big fractures with laminar flow)
		}
		else if (g13_transition_int == 1 && frac_law == 1 && transition_trigger * abs_manning_vel > abs_HP_vel) // IF manning_vel > HP_vel NO TRANSITION
		{																	
			IfmSetFracFlowConductivity(pDoc, frac_index, neu_feflow_frac_aper, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);	// Set new fracture hydraulic aperture 
			IfmSetFracArea(pDoc, frac_index, neu_frac_cross_area, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);				// Set new cross section area	
			IfmSetFracMassDecayRate(pDoc, frac_index, neu_decay_laminar, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmALL_FRAC_LAWS);		// set new decay rate (account for Diff limit on big fractures with laminar flow)
			if (g6_messages == 1){ 
				IfmInfo(pDoc, "ELSEIF #1 frac %i x=%.2f y=%.2f rate=%.2e M=%e, HP=%e NO transition", frac_index, coordi_1[0], coordi_1[1], dissolution_rate,abs_manning_vel, abs_HP_vel);}
		}

		/*1.2* neu_true_vel_manning, failsafe to avoid manning in 1st time step*/
		else if ( g13_transition_int == 1 && transition_trigger * abs_manning_vel <  abs_HP_vel ) // If Manning offers more resistance to flow than HagenPoisuille = TRANSITION
		{	
			IfmSetFracLaw(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);
			IfmSetFracArea(pDoc, frac_index, neu_frac_cross_area, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);	// Set new cross section area	
			IfmSetFracFlowConductivity(pDoc, frac_index, neu_M_corr, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);	// Set new fracture Manning coeff (area)
			IfmSetFracMassDecayRate(pDoc, frac_index, new_decay_turbulent, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);		// set new decay rate (for turbulent)
			if (g6_messages == 1) {
				IfmInfo(pDoc, "ELSEIF #2 frac %i x=%.2f y=%.2f rate=%.2e M=%e, HP=%e transition to TURB", frac_index, coordi_1[0], coordi_1[1], dissolution_rate, abs_manning_vel, abs_HP_vel); }
			
		}
		else if (frac_law == 2) // IF ALREADY MANNING, stay there
		{
			//IfmSetFracLaw(pDoc, frac_index, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);
			IfmSetFracArea(pDoc, frac_index, neu_frac_cross_area, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);	// Set new cross section area	
			IfmSetFracFlowConductivity(pDoc, frac_index, neu_M_corr, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);	// Set new fracture Manning coeff (area)
			IfmSetFracMassDecayRate(pDoc, frac_index, new_decay_turbulent, IfmALL_FRAC_TYPES, IfmALL_FRAC_MODES, IfmMANNING_LAW);		// set new decay rate (for turbulent)
			if (g6_messages == 1) {
				IfmInfo(pDoc, "ELSEIF #3 frac %i x=%.2f y=%.2f rate=%.2e M=%e, HP=%e already TURB", frac_index, coordi_1[0], coordi_1[1], dissolution_rate, abs_manning_vel, abs_HP_vel); }
		}
		else { if (g6_messages == 4) 
		{ 
			IfmInfo(pDoc, "NO growth frac %i at x=%.2f y=%.2f ave_sol=%.2f R=%e, M=%e, HP=%e", frac_index,coordi_1[0],coordi_1[1], ave_solute_2,dissolution_rate,abs_manning_vel,abs_HP_vel );}
		}
	
		/* 4.1.11 delete alloc   */
		delete[] c;							// dynamic memory allocation for vector'c' created at start of section 4 is liberated here
		
	} /* 4 ENd of big loop to update fracture */
	  /******************************************************************************************************************************/

	  /* 9 - End of dynamic memory allocations for fracture update loop */

	delete[] ptr;							// dynamic memory allocation for the pointer 'ptr' created in section 3 is liberated here
	delete[] coordi_1;						// dynamic memory allocation for vector'coordi_1' created in section 3 is liberated here
	delete[] coordi_2;						// dynamic memory allocation for vector'coordi_2' created in section 3 is liberated here

	/* 10 Time controls and clock, Part 2 */
	double end_time;
	end_time = clock();
	double clock_time;
	clock_time = end_time - start_time;
	double timeInSeconds;
	timeInSeconds = clock_time / CLOCKS_PER_SEC;
	//IfmInfo(pDoc, "waiting for reaction and transport stabilize = %g ", g98_accum_update); // FeFlow log screen message 
	if (g6_messages <= 2)
	{
		IfmInfo(pDoc, "In this time step %i fractures where not modified because C > CEQ", W_C);
		IfmInfo(pDoc, "Effective Sim time = %g days = %g years, calc time = %g clock() clicks, %.3f sec", g98_accum_update, g98_accum_update/365, clock_time, timeInSeconds); // FeFlow log screen message 
	}
	//	IfmInfo(pDoc, "Effective Sim time = %g d, calc time = %g clock() clicks, %f sec, flow = %f m3/d & Re = %f", g98_accum_update, clock_time, timeInSeconds, flow, reynolds);

	g99_accum_pert = 0 + (g99_accum_pert - g8_pert_time);					// resets g99_accum_time in section 1

	/* 11 - INIT WAIT GOTO */
	init_wait:

	if (g6_messages > 0);
	else {
		if (days < g5_init_wait) { IfmInfo(pDoc, "No frac update, simulation time %g days, waiting for initial stabilization time %g days", days, g5_init_wait); }
	}

	/* 12 - NO FRAC UPDATE GO TO*/
	no_frac_update: {; }
}

