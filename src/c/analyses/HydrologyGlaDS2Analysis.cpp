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
        iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basal_melting_rate",BasalforcingsBasalMeltingRateEnum);
    }
    if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sheet_conductivity",HydrologySheetConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.channel_conductivity",HydrologyChannelConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);
    iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",HydrologySheetWaterHeightEnum);
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
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ischannels",HydrologyIschannelsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.melt_flag",HydrologyMeltFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.channel_sheet_width",HydrologyChannelSheetWidthEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.channel_alpha",HydrologyChannelAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.channel_beta",HydrologyChannelBetaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.sheet_alpha",HydrologySheetAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.sheet_beta",HydrologySheetBetaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.omega",HydrologyOmegaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.istransition",HydrologyIsTransitionEnum));
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
    Input* h_input   = element->GetInput(HydrologySheetWaterHeightEnum); _assert_(h_input);
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
    IssmDouble Jdet,dphi[3],dh[3];
    IssmDouble* xyz_list = NULL;
    int numnodes=element->GetNumberOfNodes();
    ElementVector* pe=element->NewElementVector();
    IssmDouble* basis=xNew<IssmDouble>(numnodes);
    IssmDouble* dbasis=xNew<IssmDouble>(2*numnodes);
    element->GetVerticesCoordinates(&xyz_list);

    /*Retrieve inputs and parameters*/
    IssmDouble alpha     = element->FindParam(HydrologySheetAlphaEnum);
    IssmDouble beta      = element->FindParam(HydrologySheetBetaEnum);
    IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
    IssmDouble evr       = element->FindParam(HydrologyEnglacialVoidRatioEnum);
    IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
    IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
    IssmDouble g         = element->FindParam(ConstantsGEnum);
    Input* h_input     = element->GetInput(HydrologySheetWaterHeightEnum); _assert_(h_input);
    Input* hold_input  = element->GetInput(HydrologySheetWaterHeightOldEnum); _assert_(hold_input);
    Input* hg_input    = element->GetInput(HydrologyMeanCavityHeightEnum); _assert_(hg_input);
    Input* hw_input    = element->GetInput(HydrologyFlowingSheetHeightEnum); _assert_(hw_input);
    Input* phi_input   = element->GetInput(HydraulicPotentialEnum); _assert_(phi_input);
    Input* H_input     = element->GetInput(ThicknessEnum); _assert_(H_input);
    Input* k_input     = element->GetInput(HydrologySheetConductivityEnum); _assert_(k_input);
    Input* melt_input  = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(melt_input);

    Gauss* gauss=element->NewGauss(2);
    while(gauss->next()){
        IssmDouble h,h_old,hg,hw,H,k,melt;
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
            pe->values[i] += factor*(basis[i]*(h_old/dt+melt)
                - conductivity*(dbasis[0*numnodes+i]*residual_x
                    + dbasis[1*numnodes+i]*residual_y));
        }
    }

    xDelete<IssmDouble>(xyz_list);
    xDelete<IssmDouble>(basis);
    xDelete<IssmDouble>(dbasis);
    delete gauss;
    return pe;
}/*}}}*/
void           HydrologyGlaDS2Analysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	element->GetSolutionFromInputsOneDof(solution,HydrologySheetWaterHeightEnum);

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

	IssmDouble rho_ice   = femmodel->parameters->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = femmodel->parameters->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = femmodel->parameters->FindParam(ConstantsGEnum);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *bed       = xNew<IssmDouble>(numnodes);
		IssmDouble *thickness = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&bed[0],BaseEnum);
		element->GetInputListOnNodes(&thickness[0],ThicknessEnum);
		element->GetInputListOnNodes(&ls_active[0],HydrologyMaskNodeActivationEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]>0. && ls_active[in]==1.){
				node->Activate(); //Not sure if we need this!
			}
			else{
				IssmDouble phi =  rho_ice*g*thickness[in] + rho_water*g*bed[in]; //FIXME this is correct!
				node->Deactivate();// Not sure if we need this
				node->ApplyConstraint(0,phi);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(bed);
		xDelete<IssmDouble>(thickness);
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
    Input* h_input   = element->GetInput(HydrologySheetWaterHeightEnum); _assert_(h_input);
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


