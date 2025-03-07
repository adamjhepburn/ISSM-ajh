/*!\file: solutionsequence_nonlinear.cpp
 * \brief: core of a non-linear solution, using fixed-point method 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#define AEPS        2.2204460492503131E-015
/*local routine to check qr/lh convergence:{{{*/
bool lakelhQrconvergence(Vector<IssmDouble>* lh_new,Vector<IssmDouble>* lh_old, Vector<IssmDouble>* qr_new,Vector<IssmDouble>* qr_old, IssmDouble eps_res);
/*}}}*/

/*main routine:{{{*/
void solutionsequence_glads_nonlinear(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;
	Vector<IssmDouble>* ug  = NULL;
	Vector<IssmDouble>* uf  = NULL;
	Vector<IssmDouble>* old_uf = NULL;
	Vector<IssmDouble>* pf  = NULL;
	Vector<IssmDouble>* df  = NULL;
	Vector<IssmDouble>* ys  = NULL;

	/*lh qr convergence criterion*/
	Vector<IssmDouble>* lh_new = NULL;
	/*Vector<IssmDouble>* lh = NULL;*/
	Vector<IssmDouble>* lh_old = NULL;
	Vector<IssmDouble>* qr_new = NULL;
	/*Vector<IssmDouble>* qr = NULL;*/
	Vector<IssmDouble>* qr_old = NULL;
	GetVectorFromInputsx(&lh_old,femmodel,HydrologyLakeHeightOldEnum,VertexSIdEnum);
	GetVectorFromInputsx(&qr_old,femmodel,HydrologyLakeOutletQrOldEnum,VertexSIdEnum);
	
	/*parameters:*/
	int max_nonlinear_iterations;
	IssmDouble eps_res,eps_rel,eps_abs;
	bool islakes;
	HydrologyGlaDSAnalysis* analysis = new HydrologyGlaDSAnalysis();

	/*Recover parameters (FIXME: from Stress balance for now :( )*/
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	femmodel->UpdateConstraintsx();

	int  count_out=0;
	bool converged_out=false;
	int  count_in=0;
	bool converged_in=false;

	/*Start non-linear iteration using input velocity: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);

	while(!converged_out){
		/*save pointers to old lh and qr*/
		/*delete lh;lh=lh_new;
		delete qr;qr=qr_new;*/

		count_in=0;
		converged_in=false;
		if(islakes){
			if(VerboseConvergence()) _printf0_("   updating phi at the lake outlet\n");
			analysis->UpdateLakeOutletPhi(femmodel);
		}
		
		/*begin inner loop*/
		while(!converged_in){
			/*save pointer to old solution*/
			delete old_uf;old_uf=uf;
			delete ug;

			SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
			CreateNodalConstraintsx(&ys,femmodel->nodes);
			Reduceloadx(pf, Kfs, ys); delete Kfs;
			femmodel->profiler->Start(SOLVER);
			Solverx(&uf, Kff, pf, old_uf, df, femmodel->parameters);
			femmodel->profiler->Stop(SOLVER);
			Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete ys;

			convergence(&converged_in,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs); delete Kff; delete pf; delete df;
			InputUpdateFromSolutionx(femmodel,ug);

			/*Increase count: */
			count_in++;
			if(count_in>=max_nonlinear_iterations && !converged_in){
				//_printf0_("   maximum number of nonlinear iterations of inner loop (" << max_nonlinear_iterations << ") exceeded\n"); 
				converged_in = true;
			}
		}

		if(VerboseConvergence()) _printf0_(setw(50) << left << "   Inner loop converged in "<<count_in<<" iterations\n");

		if(VerboseConvergence()) _printf0_("   updating sheet thickness\n");
		analysis->UpdateSheetThickness(femmodel);
		if(!islakes){
			if(VerboseConvergence()) _printf0_("   updating channel cross section\n");
			analysis->UpdateChannelCrossSection(femmodel);
			/*Converged if inner loop converged in one solution*/
			if(count_in==1) converged_out = true;
		}
		
		else if(islakes){
			converged_out = true; /*initially set to true, change to false if any lakes do not converge*/
			if(VerboseConvergence()) _printf0_("   updating channel cross section and lake outlet flux\n");
			analysis->UpdateChannelCrossSection(femmodel);
			analysis->UpdateLakeOutletDischarge(femmodel);
			if(VerboseConvergence()) _printf0_("   updating lake depth\n");
			analysis->UpdateLakeDepth(femmodel);
			/*check convergence of lake height*/
			GetVectorFromInputsx(&lh_new,femmodel,HydrologyLakeHeightEnum,VertexSIdEnum);
			GetVectorFromInputsx(&qr_new,femmodel,HydrologyLakeOutletQrEnum,VertexSIdEnum);
				if(!lakelhQrconvergence(lh_new, lh_old, qr_new, qr_old, eps_res)){
					converged_out = false;
				}
		}

		/*Increase count: */
		count_out++;
		if(count_out>=max_nonlinear_iterations && !converged_out){
			_printf0_("   maximum number of nonlinear iterations of outer loop (" << max_nonlinear_iterations << ") exceeded\n"); 
			converged_out = true;
		}
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count_out<<"x"<<count_in<<"\n");

	/*clean-up*/
	delete uf;
	delete ug;
	delete old_uf;
	delete analysis;
}/*}}}*/


bool lakelhQrconvergence(Vector<IssmDouble>* lh_new, Vector<IssmDouble>* lh_old, Vector<IssmDouble>* qr_new, Vector<IssmDouble>* qr_old, IssmDouble eps_res){/*}}}*/
    bool converged = true;
    bool lh_converged = true;
    bool qr_converged = true;

    if (!xIsNan<IssmDouble>(eps_res)) {
        // Calculate the norm of the difference between lh and lh_old
        IssmDouble nlh = lh_new->Norm(NORM_TWO);
		IssmDouble nlh_old = lh_old->Norm(NORM_TWO);
		IssmDouble lh_diff_norm = fabs((nlh-nlh_old)/nlh_old+AEPS);

        if (lh_diff_norm < eps_res) {
            if (VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion met: lh/lh_old " << lh_diff_norm * 100 << " < " << eps_res * 100 << " %\n");
        } else {
            if (VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion exceeded: lh/lh_old " << lh_diff_norm * 100 << " > " << eps_res * 100 << " %\n");
             lh_converged = false;
        }
		// Calculate the norm of the difference between qr and qr_old
		IssmDouble nqr = qr_new->Norm(NORM_TWO);
		IssmDouble nqr_old = qr_old->Norm(NORM_TWO);
		IssmDouble qr_diff_norm = fabs((nqr-nqr_old)/nqr_old+AEPS);

		if (qr_diff_norm <eps_res) {
			if (VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion met: qr/qr_old " << qr_diff_norm * 100 << " < " << eps_res * 100 << " %\n");
		} else {
			if (VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion exceeded: qr/qr_old " << qr_diff_norm * 100 << " > " << eps_res * 100 << " %\n");
			qr_converged = false;
		}
	}
	
	if (!lh_converged || !qr_converged) converged = false;

	/*assign output*/
	return converged;
} /*}}}*/
