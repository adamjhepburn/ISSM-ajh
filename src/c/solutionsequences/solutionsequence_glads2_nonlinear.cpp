/*!\file: solutionsequence_glads2_nonlinear.cpp
 * \brief: core of the GlaDS2 sheet thickness solve, using fixed-point iteration
 */

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_glads2_nonlinear(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;
	Vector<IssmDouble>* ug  = NULL;
	Vector<IssmDouble>* uf  = NULL;
	Vector<IssmDouble>* old_uf = NULL;
	Vector<IssmDouble>* pf  = NULL;
	Vector<IssmDouble>* df  = NULL;
	Vector<IssmDouble>* ys  = NULL;

	/*parameters:*/
	int max_nonlinear_iterations;
	IssmDouble eps_res,eps_rel,eps_abs;
	HydrologyGlaDS2Analysis* analysis = new HydrologyGlaDS2Analysis();

	/*Recover parameters (FIXME: from Stress balance for now :( )*/
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->UpdateConstraintsx();

	int  count=0;
	bool converged=false;

	/*Start non-linear iteration using input sheet thickness: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);

	while(!converged){

		/*save pointer to old solution*/
		delete old_uf;old_uf=uf;
		delete ug;

		/*Diagnostic chain: mean cavity height is held fixed at its old value here,
		 *it is only updated once the sheet thickness itself has converged (see hydrology_core.cpp)*/
		analysis->UpdateWaterPressure(femmodel);
		analysis->UpdateFlowingSheetHeight(femmodel);
		analysis->UpdateHydraulicPotential(femmodel);

		SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys); delete Kfs;
		femmodel->profiler->Start(SOLVER);
		Solverx(&uf, Kff, pf, old_uf, df, femmodel->parameters);
		femmodel->profiler->Stop(SOLVER);
		Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete ys;

		convergence(&converged,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs); delete Kff; delete pf; delete df;
		InputUpdateFromSolutionx(femmodel,ug);

		/*Increase count: */
		count++;
		if(count>=max_nonlinear_iterations && !converged){
			_printf0_("   maximum number of nonlinear iterations ("<<max_nonlinear_iterations<<") exceeded\n");
			converged=true;
		}
	}
	if(VerboseConvergence()) _printf0_(setw(50) << left << "   converged in "<<count<<" iterations\n");

	/*Bring pw/hw/phi in sync with the converged sheet thickness before returning*/
	analysis->UpdateWaterPressure(femmodel);
	analysis->UpdateFlowingSheetHeight(femmodel);
	analysis->UpdateHydraulicPotential(femmodel);

	/*clean-up*/
	delete uf;
	delete ug;
	delete old_uf;
	delete analysis;
}
