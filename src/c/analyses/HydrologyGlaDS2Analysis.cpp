#include <float.h> /* defines DBL_EPSILON*/
#include "./HydrologyGlaDS2Analysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyGlaDS2Analysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologyGlaDS2Enum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spch",HydrologyGlaDS2AnalysisEnum,P1Enum);

}/*}}}*/
void HydrologyGlaDS2Analysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

    /*Fetch parameters: */
    int  hydrology_model;
    iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

    /*Now, do we really want GlaDS2?*/
    if(hydrology_model!=HydrologyGlaDS2Enum) return;

    /*Add channels? UNDER CONSTRUCTION*/

    /*Create discrete loads for Moulins UNDER CONSTRUCTION*/

	/*Deal with Neumann BC*/
	int M,N;
	int *segments = NULL;
	iomodel->FetchData(&segments,&M,&N,"md.mesh.segments");

	/*Check that the size seem right*/
	_assert_(N==3); _assert_(M>=3);
	for(int i=0;i<M;i++){
		if(iomodel->my_elements[segments[i*3+2]-1]){
			loads->AddObject(new Neumannflux(i+1,i,iomodel,segments));
		}
	}
	xDelete<int>(segments);

}/*}}}*/
void HydrologyGlaDS2Analysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

    /*Fetch parameters: */
    int  hydrology_model;
    iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

    /*Now, do we really want GlaDS2?*/
    if(hydrology_model!=HydrologyGlaDS2Enum) return;

    if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
    ::CreateNodes(nodes,iomodel,HydrologyGlaDS2AnalysisEnum,P1Enum);
    iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  HydrologyGlaDS2Analysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
    return 1;
}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

    /*Fetch parameters: */
    int  hydrology_model;
    iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
    int    meltflag;	
	iomodel->FindConstant(&meltflag,"md.hydrology.melt_flag");

    /*Now, do we really want GlaDS2?*/
    if(hydrology_model!=HydrologyGlaDS2Enum) return;

    /*Update elements: */
    int counter=0;
    for(int i=0;i<iomodel->numberofelements;i++){
        if(iomodel->my_elements[i]){
            Element* element=(Element*)elements->GetObjectByOffset(counter);
            element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
            counter++;
        }
    }
    
    iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
    if(meltflag==2){
        iomodel->FetchDataToInput(inputs,elements,"md.smb.runoff",SmbRunoffEnum);
    }
    if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sheet_conductivity",HydrologySheetConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);
    iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",HydrologySheetHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.hydraulic_potential",HydraulicPotentialEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.rheology_B_base",HydrologyRheologyBBaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum);
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);
}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want GlaDS2?*/
	if(hydrology_model!=HydrologyGlaDS2Enum) return;
    parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.pressure_melt_coefficient",HydrologyPressureMeltCoefficientEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.cavity_spacing",HydrologyCavitySpacingEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.melt_flag",HydrologyMeltFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.sheet_alpha",HydrologySheetAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.sheet_beta",HydrologySheetBetaEnum));
    parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.englacial_void_ratio",HydrologyEnglacialVoidRatioEnum));

    /*Friction*/
	FrictionUpdateParameters(parameters, iomodel);

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");
}/*}}}*/
/*Finite Element Analysis*/
void           HydrologyGlaDS2Analysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyGlaDS2Analysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyGlaDS2Analysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyGlaDS2Analysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyGlaDS2Analysis::CreateKMatrix(Element* element){/*{{{*/
    /*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
    IssmDouble  Jdet,dphi[3];
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
    element->GetVerticesCoordinates(&xyz_list);

    /*Retrieve inputs and parameters from the current Picard iterate*/
    IssmDouble alpha     = element->FindParam(HydrologySheetAlphaEnum);
    IssmDouble beta      = element->FindParam(HydrologySheetBetaEnum);
    IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
    IssmDouble evr       = element->FindParam(HydrologyEnglacialVoidRatioEnum);
    IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
    IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
    IssmDouble g         = element->FindParam(ConstantsGEnum);
    Input* h_input   = element->GetInput(HydrologySheetHeightEnum); _assert_(h_input);
    Input* hg_input  = element->GetInput(HydrologyMeanCavityHeightEnum); _assert_(hg_input);
    Input* hw_input  = element->GetInput(HydrologyFlowingSheetHeightEnum); _assert_(hw_input);
    Input* phi_input = element->GetInput(HydraulicPotentialEnum); _assert_(phi_input);
    Input* H_input   = element->GetInput(ThicknessEnum); _assert_(H_input);
    Input* k_input   = element->GetInput(HydrologySheetConductivityEnum); _assert_(k_input);

    /*Backward-Euler mass term and implicit Picard flux term*/
    Gauss* gauss=element->NewGauss(2);
    while(gauss->next()){
        IssmDouble h,hg,hw,H,k;
        element->JacobianDeterminant(&Jdet,xyz_list,gauss);
        element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
        element->NodalFunctions(basis,gauss);
        h_input->GetInputValue(&h,gauss);
        hg_input->GetInputValue(&hg,gauss);
        hw_input->GetInputValue(&hw,gauss);
        H_input->GetInputValue(&H,gauss);
        k_input->GetInputValue(&k,gauss);
        phi_input->GetInputDerivativeValue(&dphi[0],xyz_list,gauss);

        IssmDouble normgradphi=sqrt(dphi[0]*dphi[0]+dphi[1]*dphi[1]);
        if(normgradphi<DBL_EPSILON) normgradphi=DBL_EPSILON;
        IssmDouble conductivity=k*pow(max(hw,DBL_EPSILON),alpha)*pow(normgradphi,beta-2.);

        /*Derivative of the regularized pressure closure used in UpdateWaterPressure*/
        IssmDouble c=rho_water*g;
        IssmDouble storage_capacity=evr*(rho_ice/rho_water)*H;
        if(storage_capacity>DBL_EPSILON){
            IssmDouble x=(h-hg)/storage_capacity;
            IssmDouble delta=0.01;
            IssmDouble drx1=0.;
            IssmDouble drx2=0.;
            if(x>-delta && x<delta) drx1=(x+delta)/(2.*delta);
            else if(x>=delta) drx1=1.;
            if(x-1.>-delta && x-1.<delta) drx2=(x-1.+delta)/(2.*delta);
            else if(x-1.>=delta) drx2=1.;
            IssmDouble dpdh=rho_ice*g*H/storage_capacity*(drx1-drx2);
            c+= (1.-evr)*dpdh;
        }

        IssmDouble factor=gauss->weight*Jdet;
        for(int i=0;i<numnodes;i++){
            for(int j=0;j<numnodes;j++){
                Ke->values[i*numnodes+j] += factor*(
                    basis[i]*basis[j]/dt
                    + conductivity*c*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j]
                        + dbasis[1*numnodes+i]*dbasis[1*numnodes+j]));
            }
        }
    }

    xDelete<IssmDouble>(xyz_list);
    xDelete<IssmDouble>(basis);
    xDelete<IssmDouble>(dbasis);
    delete gauss;
    return Ke;

}/*}}}*/
ElementVector* HydrologyGlaDS2Analysis::CreatePVector(Element* element){/*{{{*/

    /*Skip if water or ice shelf element*/
    if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;

    /*Intermediaries*/
    int meltflag;
    IssmDouble Jdet,dphi[3],dh[3];
    IssmDouble h,h_old,hg,hw,H,k,melt,G,RO,m;
    IssmDouble vx,vy,ub,frictionheat,alpha2;
    IssmDouble* xyz_list = NULL;
    int numnodes=element->GetNumberOfNodes();
    ElementVector* pe=element->NewElementVector();
    IssmDouble* basis=xNew<IssmDouble>(numnodes);
    IssmDouble* dbasis=xNew<IssmDouble>(2*numnodes);
    element->GetVerticesCoordinates(&xyz_list);

    /*Retrieve inputs and parameters*/
    element->FindParam(&meltflag,HydrologyMeltFlagEnum);
    IssmDouble alpha     = element->FindParam(HydrologySheetAlphaEnum);
    IssmDouble beta      = element->FindParam(HydrologySheetBetaEnum);
    IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
    IssmDouble evr       = element->FindParam(HydrologyEnglacialVoidRatioEnum);
    IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
    IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
    IssmDouble g         = element->FindParam(ConstantsGEnum);
    IssmDouble L         = element->FindParam(MaterialsLatentheatEnum);
    Input* h_input     = element->GetInput(HydrologySheetHeightEnum); _assert_(h_input);
    Input* hold_input  = element->GetInput(HydrologySheetHeightOldEnum); _assert_(hold_input);
    Input* hg_input    = element->GetInput(HydrologyMeanCavityHeightEnum); _assert_(hg_input);
    Input* hw_input    = element->GetInput(HydrologyFlowingSheetHeightEnum); _assert_(hw_input);
    Input* phi_input   = element->GetInput(HydraulicPotentialEnum); _assert_(phi_input);
    Input* H_input     = element->GetInput(ThicknessEnum); _assert_(H_input);
    Input* k_input     = element->GetInput(HydrologySheetConductivityEnum); _assert_(k_input);
    Input* melt_input  = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(melt_input);
    Input* G_input     = element->GetInput(BasalforcingsGeothermalfluxEnum); _assert_(G_input);

    /*Friction, needed for meltflag==0 (geothermal+frictional heat)*/
    Friction* friction=new Friction(element,2);

    Gauss* gauss=element->NewGauss(2);
    while(gauss->next()){

        element->JacobianDeterminant(&Jdet,xyz_list,gauss);
        element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
        element->NodalFunctions(basis,gauss);
        h_input->GetInputValue(&h,gauss);
        hold_input->GetInputValue(&h_old,gauss);
        hg_input->GetInputValue(&hg,gauss);
        hw_input->GetInputValue(&hw,gauss);
        H_input->GetInputValue(&H,gauss);
        k_input->GetInputValue(&k,gauss);
        melt_input->GetInputValue(&melt,gauss);
        h_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
        phi_input->GetInputDerivativeValue(&dphi[0],xyz_list,gauss);

		/*Compute melt m, following md.hydrology.melt_flag*/
		if(meltflag == 0){
			G_input->GetInputValue(&G,gauss);
			friction->GetBasalSlidingSpeeds(&vx, &vy ,gauss);
			ub = sqrt(vx*vx + vy*vy);
			friction->GetAlpha2(&alpha2,gauss);
			frictionheat=alpha2*ub*ub;
			m = (G + frictionheat)/(rho_ice*L);
		}
		else if(meltflag == 1){
			m = melt;
		}
		else{
			Input* RO_input = element->GetInput(SmbRunoffEnum);_assert_(RO_input);
			RO_input->GetInputValue(&RO,gauss);
			m = melt + RO;
		}


        IssmDouble normgradphi=sqrt(dphi[0]*dphi[0]+dphi[1]*dphi[1]);
        if(normgradphi<DBL_EPSILON) normgradphi=DBL_EPSILON;
        IssmDouble conductivity=k*pow(max(hw,DBL_EPSILON),alpha)*pow(normgradphi,beta-2.);

        IssmDouble c=rho_water*g;
        IssmDouble storage_capacity=evr*(rho_ice/rho_water)*H;
        if(storage_capacity>DBL_EPSILON){
            IssmDouble x=(h-hg)/storage_capacity;
            IssmDouble delta=0.01;
            IssmDouble drx1=0.;
            IssmDouble drx2=0.;
            if(x>-delta && x<delta) drx1=(x+delta)/(2.*delta);
            else if(x>=delta) drx1=1.;
            if(x-1.>-delta && x-1.<delta) drx2=(x-1.+delta)/(2.*delta);
            else if(x-1.>=delta) drx2=1.;
            IssmDouble dpdh=rho_ice*g*H/storage_capacity*(drx1-drx2);
            c+=(1.-evr)*dpdh;
        }

        IssmDouble residual_x=dphi[0]-c*dh[0];
        IssmDouble residual_y=dphi[1]-c*dh[1];
        IssmDouble factor=gauss->weight*Jdet;
        for(int i=0;i<numnodes;i++){
            pe->values[i] += factor*(basis[i]*(h_old/dt+m)
                - conductivity*(dbasis[0*numnodes+i]*residual_x
                    + dbasis[1*numnodes+i]*residual_y));
        }
    }

    xDelete<IssmDouble>(xyz_list);
    xDelete<IssmDouble>(basis);
    xDelete<IssmDouble>(dbasis);
    delete gauss;
    delete friction;
    return pe;
}/*}}}*/
void           HydrologyGlaDS2Analysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	element->GetSolutionFromInputsOneDof(solution,HydrologySheetWaterHeightEnum);

    /*Compute hydrology vx and vy for timestepping purposes, store sheet discharge for mean cavity height eq.*/

    /*Intermediaries*/
    IssmDouble  dphi[3],hw,k,phi;
	IssmDouble  h_r;
	IssmDouble  oceanLS,iceLS;
	IssmDouble* xyz_list = NULL;

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize water sheet velocity and discharge*/
	IssmDouble* vx = xNew<IssmDouble>(numvertices);
	IssmDouble* vy = xNew<IssmDouble>(numvertices);
	IssmDouble* d = xNew<IssmDouble>(numvertices);

    /*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numvertices;iv++) vx[iv] = 0.;
		for(int iv=0;iv<numvertices;iv++) vy[iv] = 0.;
		for(int iv=0;iv<numvertices;iv++) d[iv] = 0.;
		element->AddInput(HydrologyWaterVxEnum,vx,P1DGEnum);
		element->AddInput(HydrologyWaterVyEnum,vy,P1DGEnum);
		element->AddInput(HydrologySheetDischargeEnum,d,P1DGEnum);
		xDelete<IssmDouble>(vx);
		xDelete<IssmDouble>(vy);
		xDelete<IssmDouble>(d);
		return;
	}

    /*Retrieve all inputs and parameters*/
	bool istransition;
	element->FindParam(&istransition,HydrologyIsTransitionEnum);
	IssmDouble alpha     = element->FindParam(HydrologySheetAlphaEnum);
	IssmDouble beta      = element->FindParam(HydrologySheetBetaEnum);
	/*IssmDouble omega     = element->FindParam(HydrologyOmegaEnum);*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble mu_water  = element->FindParam(MaterialsMuWaterEnum);
	Input *k_input       = element->GetInput(HydrologySheetConductivityEnum); _assert_(k_input);
	Input *phi_input     = element->GetInput(HydraulicPotentialEnum);         _assert_(phi_input);
	Input *hr_input      = element->GetInput(HydrologyBumpHeightEnum);        _assert_(hr_input);
	Input *hw_input       = element->GetInput(HydrologyFlowingSheetHeightEnum);    _assert_(hw_input);
	Input *oceanLS_input = element->GetInput(MaskOceanLevelsetEnum);          _assert_(oceanLS_input);
	Input *iceLS_input   = element->GetInput(MaskIceLevelsetEnum);            _assert_(iceLS_input);


    /* Start looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		/*Get input values at gauss points*/
      phi_input->GetInputDerivativeValue(&dphi[0],xyz_list,gauss);
      phi_input->GetInputValue(&phi,gauss);
      h_input->GetInputValue(&h,gauss);
      hr_input->GetInputValue(&h_r,gauss); 
      k_input->GetInputValue(&k,gauss);
		oceanLS_input->GetInputValue(&oceanLS,gauss);
		iceLS_input->GetInputValue(&iceLS,gauss);

		/*Set to zero if floating or no ice*/
		if(oceanLS<0. || iceLS>0.){
			vx[iv] = 0.;
         vy[iv] = 0.;
			d[iv] = 0.;
		}
		else{

         /*Get norm of gradient of hydraulic potential and make sure it is >0*/
         IssmDouble normgradphi = sqrt(dphi[0]*dphi[0] + dphi[1]*dphi[1]);
         if(normgradphi < DBL_EPSILON) normgradphi = DBL_EPSILON;

         IssmDouble coeff;
         /*If omega is zero, use standard model, otherwise transition model*/
         /*IssmDouble nu = mu_water/rho_water;
			IssmDouble coeff;
			if(istransition==1 && omega>=DBL_EPSILON){
				IssmDouble hratio = fabs(h/h_r);
				IssmDouble coarg = 1. + 4.*pow(hratio,3-2*alpha)*omega*k*pow(h,3)*normgradphi/nu;
				coeff = nu/2./omega*pow(hratio,2*alpha-3) * (-1 + pow(coarg, 0.5))/normgradphi;  // coeff gives discharge; divide by h to get speed instead of discharge
			}
			else {*/
			
            coeff = k*pow(hw,alpha)*pow(normgradphi,beta-2.);  // coeff gives discharge; divide by h to get speed instead of discharge
			}

			vx[iv] = -coeff/max(DBL_EPSILON,hw)*dphi[0];
			vy[iv] = -coeff/max(DBL_EPSILON,hw)*dphi[1];

			d[iv] = coeff*normgradphi;
	}

	element->AddInput(HydrologyWaterVxEnum,vx,P1DGEnum);
	element->AddInput(HydrologyWaterVyEnum,vy,P1DGEnum);
	element->AddInput(HydrologySheetDischargeEnum,d,P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(d);
	delete gauss;

}/*}}}*/
void           HydrologyGlaDS2Analysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyGlaDS2Analysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	element->InputUpdateFromSolutionOneDof(solution,HydrologySheetWaterHeightEnum);
}/*}}}*/

void HydrologyGlaDS2Analysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Update active elements based on ice levelset and ocean levelset*/
	GetMaskOfIceVerticesLSMx(femmodel,true);
	SetActiveNodesLSMx(femmodel,true);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&ls_active[0],HydrologyMaskNodeActivationEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]>0. && ls_active[in]==1.){
				node->Activate(); //Not sure if we need this!
			}
			else{
				/*This analysis solves for h (sheet water height), so floating/inactive nodes are constrained to h=0*/
				node->Deactivate();// Not sure if we need this
				node->ApplyConstraint(0,0.);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(ls_active);
	}

	return;
}/*}}}*/

/*GlaDS specifics*/
void HydrologyGlaDS2Analysis::UpdateWaterPressure(FemModel* femmodel){/*{{{*/

    for(Object* & object : femmodel->elements->objects){
        Element* element=xDynamicCast<Element*>(object);
        UpdateWaterPressure(element);
    }

}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateWaterPressure(Element* element){/*{{{*/

    /*Intermediary*/
    IssmDouble h, hg, dh;
    IssmDouble H; 
    IssmDouble pi, x, delta, rx1, rx2;
    IssmDouble oceanLS,iceLS;

    /*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

    /*Initialize new water pressure*/
    IssmDouble* pw_new = xNew<IssmDouble>(numvertices);

    /*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numvertices;iv++) {
			pw_new[iv] = 0.;
		}
		element->AddInput(HydrologyWaterPressureEnum,pw_new,P1Enum);
		xDelete<IssmDouble>(pw_new);
		return;
	}
    
    /*Retrieve all inputs and parameters*/
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
    IssmDouble evr       = element->FindParam(HydrologyEnglacialVoidRatioEnum);
    Input* H_input   = element->GetInput(ThicknessEnum); _assert_(H_input);
    Input* h_input   = element->GetInput(HydrologySheetWaterHeightEnum); _assert_(h_input);
    Input* hg_input  = element->GetInput(HydrologyMeanCavityHeightEnum); _assert_(hg_input);
    Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum); _assert_(oceanLS_input);
    Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum); _assert_(iceLS_input);
    
    /* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

        /*Get input values at gauss points*/
        h_input->GetInputValue(&h,gauss);
        hg_input->GetInputValue(&hg,gauss);
        H_input->GetInputValue(&H,gauss);
		oceanLS_input->GetInputValue(&oceanLS,gauss);
		iceLS_input->GetInputValue(&iceLS,gauss);

        /*Set water pressure to 0 if floating or no ice*/
        if(oceanLS<0. || iceLS>0.){
            pw_new[iv] = 0.;
        }
        else{
            /*Compute max water storage*/
            dh = evr*(rho_ice/rho_water)*H;
            if (dh <= 0.){
                pw_new[iv] = 0.;
            }
            else{

                /*Get ice pressure*/
                pi = rho_ice*g*H;

                /*Compute water pressure*/
                /*note this is a C^1 continuous regularised version of piecewise relation (eq.4) in Wells et al., 2026 */
                x = (h-hg)/dh;
                delta = 0.01;
                
                if (x <= -delta) {
                    rx1 = 0.0;
                }
                else if (x < delta) {
                    rx1 = (x + delta) * (x + delta) / (4.0 * delta);
                }
                else {
                    rx1 = x;
                }

                if (x -1.0 < -delta) {
                    rx2 = 0.0;
                }
                else if (x -1.0 < delta) {
                    rx2 = ((x - 1.0) + delta) * ((x - 1.0) + delta) / (4.0 * delta);
                }
                else {
                    rx2 = x - 1.0;
                }
                pw_new[iv] = pi * (rx1 - rx2);

            }
            }
        }
    
        element->AddInput(HydrologyWaterPressureEnum,pw_new,P1Enum);
        /*Clean up and return*/
        delete gauss;
        xDelete<IssmDouble>(pw_new);

}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateFlowingSheetHeight(FemModel* femmodel){/*{{{*/

    for(Object* & object : femmodel->elements->objects){
        Element* element=xDynamicCast<Element*>(object);
        UpdateFlowingSheetHeight(element);
    }

}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateFlowingSheetHeight(Element* element){/*{{{*/

    /*Intermediares*/
    IssmDouble h, pw;
    IssmDouble oceanLS,iceLS;

    int numvertices = element->GetNumberOfVertices();

    /*Initialize new water pressure*/
    IssmDouble* hw_new = xNew<IssmDouble>(numvertices);
    /*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numvertices;iv++) {
			hw_new[iv] = 0.;
		}
		element->AddInput(HydrologyFlowingSheetHeightEnum,hw_new,P1Enum);
		xDelete<IssmDouble>(hw_new);
		return;
	}
    
    /*Retrieve all inputs and parameters*/
    IssmDouble evr       = element->FindParam(HydrologyEnglacialVoidRatioEnum);
    IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
    Input* h_input   = element->GetInput(HydrologySheetHeightEnum); _assert_(h_input);
    Input* pw_input  = element->GetInput(HydrologyWaterPressureEnum); _assert_(pw_input);
    Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum); _assert_(oceanLS_input);
	Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum); _assert_(iceLS_input);

    Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);
        /*Get input values at gauss points*/
        h_input->GetInputValue(&h,gauss);
        pw_input->GetInputValue(&pw,gauss);
        oceanLS_input->GetInputValue(&oceanLS,gauss);
        iceLS_input->GetInputValue(&iceLS,gauss);

        /*Set flowing water height to zero if floating or no ice*/
        if (element->IsAllFloating() || !element->IsIceInElement()) {
            hw_new[iv] = 0.;
        }
        else {
            /*Calculate regularised flowing water height*/
            hw_new[iv] = h - evr * pw / (rho_water * g);
            if (hw_new[iv] < 0.) {
                hw_new[iv] = 0.;
            }
        }
    }
    element->AddInput(HydrologyFlowingSheetHeightEnum,hw_new,P1Enum);
    xDelete<IssmDouble>(hw_new);
    delete gauss;

}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateHydraulicPotential(FemModel* femmodel){/*{{{*/

    for(Object* & object : femmodel->elements->objects){
        Element* element=xDynamicCast<Element*>(object);
        UpdateHydraulicPotential(element);
    }

}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateHydraulicPotential(Element* element){/*{{{*/
    
    /*Intermediares*/
    IssmDouble zb, hw, pw;
    IssmDouble oceanLS,iceLS;

    int numvertices = element->GetNumberOfVertices();

    /*Initialize new hydraulic potential*/
    IssmDouble* phi_new = xNew<IssmDouble>(numvertices);

    /*set to 0 if inactive element*/
    if(element->IsAllFloating() || !element->IsIceInElement()){
        for(int iv=0;iv<numvertices;iv++) {
            phi_new[iv] = 0.;
        }
        element->AddInput(HydraulicPotentialEnum,phi_new,P1Enum);
        xDelete<IssmDouble>(phi_new);
        return;
    }

    /*Retrieve all inputs and parameters*/
    IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
    IssmDouble g         = element->FindParam(ConstantsGEnum);
    Input* zb_input  = element->GetInput(BedEnum); _assert_(zb_input);
    Input* hw_input   = element->GetInput(HydrologyFlowingSheetHeightEnum); _assert_(hw_input);
    Input* pw_input  = element->GetInput(HydrologyWaterPressureEnum); _assert_(pw_input);
    Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum); _assert_(oceanLS_input);
    Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum); _assert_(iceLS_input);

    /* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

        /*Get input values at gauss points*/
        zb_input->GetInputValue(&zb,gauss);
        hw_input->GetInputValue(&hw,gauss);
        pw_input->GetInputValue(&pw,gauss);
        oceanLS_input->GetInputValue(&oceanLS,gauss);
        iceLS_input->GetInputValue(&iceLS,gauss);

        /*Set hydraulic potential to 0 if floating or no ice*/
        if(oceanLS<0. || iceLS>0.){
            phi_new[iv] = 0.;
        }
        else{
            phi_new[iv] = pw + rho_water*g*(zb+hw);
        }
    }

    element->AddInput(HydraulicPotentialEnum,phi_new,P1Enum);
    xDelete<IssmDouble>(phi_new);
    delete gauss;

}/*}}}*/

void HydrologyGlaDS2Analysis::UpdateMeanCavityHeight(FemModel* femmodel){/*{{{*/

    for(Object* & object : femmodel->elements->objects){
        Element* element=xDynamicCast<Element*>(object);
        UpdateMeanCavityHeight(element);
    }

}/*}}}*/
void HydrologyGlaDS2Analysis::UpdateMeanCavityHeight(Element* element){/*{{{*/

    /*Intermediaries */
	IssmDouble  vx,vy,ub,hg_old,N,h_r,H,b;
	IssmDouble  A,B,n,phi,phi_0;
	IssmDouble  alpha,beta;
	IssmDouble  oceanLS,iceLS;
    
    /*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize new sheet thickness*/
	IssmDouble* hg_new = xNew<IssmDouble>(numvertices);

    /*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numvertices;iv++) {
			hg_new[iv] = 0.;
		}
		element->AddInput(HydrologyMeanCavityHeightEnum,hg_new,P1Enum);
		xDelete<IssmDouble>(hg_new);
		return;
	}
    /*Retrieve all inputs and parameters*/
	IssmDouble  dt       = element->FindParam(TimesteppingTimeStepEnum);
	IssmDouble  l_r      = element->FindParam(HydrologyCavitySpacingEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	Input* hr_input = element->GetInput(HydrologyBumpHeightEnum);         _assert_(hr_input);
	Input* vx_input = element->GetInput(VxBaseEnum);                      _assert_(vx_input);
	Input* vy_input = element->GetInput(VyBaseEnum);                      _assert_(vy_input);
	Input* H_input = element->GetInput(ThicknessEnum);                    _assert_(H_input);
	Input* b_input = element->GetInput(BedEnum);                          _assert_(b_input);
	Input* hgold_input = element->GetInput(HydrologyMeanCavityHeightOldEnum);_assert_(hgold_input);
	Input* B_input = element->GetInput(HydrologyRheologyBBaseEnum);       _assert_(B_input);
	Input* n_input = element->GetInput(MaterialsRheologyNEnum);           _assert_(n_input);
	Input* phi_input = element->GetInput(HydraulicPotentialEnum);         _assert_(phi_input);
	Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum);      _assert_(oceanLS_input);
	Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum);          _assert_(iceLS_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		/*Get input values at gauss points*/
		phi_input->GetInputValue(&phi,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		hgold_input->GetInputValue(&hg_old,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		hr_input->GetInputValue(&h_r,gauss);
		b_input->GetInputValue(&b,gauss);
		H_input->GetInputValue(&H,gauss);
		oceanLS_input->GetInputValue(&oceanLS,gauss);
		iceLS_input->GetInputValue(&iceLS,gauss);

		/*Set sheet thickness to zero if floating or no ice*/
		if(oceanLS<0. || iceLS>0.){
			hg_new[iv] = 0.;
		}
		else{

		/*Get values for a few potentials*/
		phi_0   = rho_water*g*b + rho_ice*g*H;
		N = phi_0 - phi;

		/*Get basal velocity*/
		ub = sqrt(vx*vx + vy*vy);

		/*Get A from B and n*/
		A = pow(B,-n);
        alpha = -ub/l_r - 2./pow(n,n)*A*pow(fabs(N),n-1.)*N;
        beta  = ub*h_r/l_r;
		
		/*Get new sheet thickness*/
		hg_new[iv] = ODE1(alpha,beta,hg_old,dt,1);

		/*Make sure it is positive*/
		if(hg_new[iv]<DBL_EPSILON) hg_new[iv] = DBL_EPSILON;
		
		}

	}

	element->AddInput(HydrologyMeanCavityHeightEnum,hg_new,P1Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(hg_new);
	delete gauss;
}/*}}}*/





