#include "./HydrologyGlaDSAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

#define AEPS 2.2204460492503131E-015

/*Model processing*/
void HydrologyGlaDSAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologyGlaDSEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spcphi",HydrologyGlaDSAnalysisEnum,P1Enum);

	bool islakes;
	iomodel->FindConstant(&islakes,"md.hydrology.islakes");
	if(islakes){
		/*Fetch lake datalake*/
		IssmDouble* lake_data = NULL;
		IssmDouble* spcphi    = NULL;
		IssmDouble* phi       = NULL;
		iomodel->FetchData(&lake_data,NULL,NULL,"md.hydrology.lake_mask");
		iomodel->FetchData(&spcphi,NULL,NULL,"md.hydrology.spcphi");
		iomodel->FetchData(&phi,NULL,NULL,"md.initialization.hydraulic_potential");

		/*Create spcs from x,y,z, as well as the spc values on those spcs: */
		int count = constraints->Size();
		for(int i=0;i<iomodel->numberofvertices;i++){
			if(iomodel->my_vertices[i]){
				if(xIsNan(spcphi[i]) && lake_data[i]>0.){
					constraints->AddObject(new SpcDynamic(count+1,i+1,0,phi[i], HydrologyGlaDSAnalysisEnum));
					count++;
				}
			}
		}
	
		/*Cleanup*/
		xDelete<IssmDouble>(lake_data);
		xDelete<IssmDouble>(spcphi);
		xDelete<IssmDouble>(phi);
	}




}/*}}}*/
void HydrologyGlaDSAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want GlaDS?*/
	if(hydrology_model!=HydrologyGlaDSEnum) return;

	/*Add channels?*/
	int Ka,La;
	int Kd,Ld;
	bool ischannels;
	IssmDouble* channelarea;
	IssmDouble* channeldischarge;
	iomodel->FindConstant(&ischannels,"md.hydrology.ischannels");
	iomodel->FetchData(&channelarea,&Ka,&La,"md.initialization.channelarea");
	iomodel->FetchData(&channeldischarge,&Kd,&Ld,"md.initialization.channel_discharge");
	if(ischannels){
		/*Get faces (edges in 2d)*/
		CreateFaces(iomodel);
		for(int i=0;i<iomodel->numberoffaces;i++){

			/*Get left and right elements*/
			int element=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]

			/*Now, if this element is not in the partition*/
			if(!iomodel->my_elements[element]) continue;

			/*Add channelarea and discharge from initialization if exists*/
			if(Ka==0 && Kd==0){
				loads->AddObject(new Channel(i+1,0.,0.,i,iomodel));
			}
			if(Ka==0 && Kd==iomodel->numberoffaces){
				loads->AddObject(new Channel(i+1,0.,channeldischarge[i],i,iomodel));
			}
			if(Ka==iomodel->numberoffaces && Kd==0){
				loads->AddObject(new Channel(i+1,channelarea[i],0.,i,iomodel));
			}
			if(Ka==iomodel->numberoffaces && Kd==iomodel->numberoffaces){
				loads->AddObject(new Channel(i+1,channelarea[i],channeldischarge[i],i,iomodel));
			}
			else{
				if (Ka!=0 && Ka!=iomodel->numberoffaces) _error_("Unknown dimension for channel area initialization.");
				if (Kd!=0 && Kd!=iomodel->numberoffaces) _error_("Unknown dimension for channel discharge initialization.");
			}
			iomodel->DeleteData(1,"md.initialization.channelarea");
            iomodel->DeleteData(1,"md.initialization.channel_discharge");
		}

	}
	/*iomodel->DeleteData(1,"md.initialization.channelarea");
	iomodel->DeleteData(1,"md.initialization.channel_discharge");*/
	xDelete<IssmDouble>(channelarea);
    xDelete<IssmDouble>(channeldischarge);


	/*Create discrete loads for Moulins*/
	CreateSingleNodeToElementConnectivity(iomodel);
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(1,"md.mesh.vertexonbase");
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}	
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
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
void HydrologyGlaDSAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want GlaDS?*/
	if(hydrology_model!=HydrologyGlaDSEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologyGlaDSAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  HydrologyGlaDSAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	int    meltflag;	
	iomodel->FindConstant(&meltflag,"md.hydrology.melt_flag");
	bool    islakes;
	iomodel->FindConstant(&islakes,"md.hydrology.islakes");

	/*Now, do we really want GlaDS?*/
	if(hydrology_model!=HydrologyGlaDSEnum) return;

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
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonboundary",MeshVertexonboundaryEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sheet_conductivity",HydrologySheetConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.channel_conductivity",HydrologyChannelConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.englacial_void_ratio",HydrologyEnglacialVoidRatioEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.moulin_input",HydrologyMoulinInputEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",HydrologySheetThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.hydraulic_potential",HydraulicPotentialEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.rheology_B_base",HydrologyRheologyBBaseEnum);
	if(islakes){
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.lake_mask",HydrologyLakeMaskEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.lake_area",HydrologyLakeAreaEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.lake_outletQr",HydrologyLakeOutletQrEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.lake_depth",HydrologyLakeHeightEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.lake_Qin",HydrologyLakeQinEnum);
	}
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum);
	}
	else{
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);

}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");


	/*Now, do we really want GlaDS?*/
	if(hydrology_model!=HydrologyGlaDSEnum) return;

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
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.isincludesheetthickness",HydrologyIsIncludeSheetThicknessEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.creep_open_flag",HydrologyCreepOpenFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.islakes",HydrologyLakeFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.num_lakes",HydrologyNumLakesEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.mean_edge_length",HydrologyMeanEdgeLengthEnum));
	

	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyGlaDSAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyGlaDSAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyGlaDSAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyGlaDSAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyGlaDSAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,dphi[3],h,k,ev;
	IssmDouble  h_r;
	IssmDouble  A,B,n,phi_old,phi,phi_0,H,b,v1;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/*Get all inputs and parameters*/
	bool istransition;
	bool isincludesheetthickness;
	bool creep_open_flag;
	element->FindParam(&istransition,HydrologyIsTransitionEnum);
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	element->FindParam(&creep_open_flag,HydrologyCreepOpenFlagEnum);
	IssmDouble alpha     = element->FindParam(HydrologySheetAlphaEnum);
	IssmDouble beta      = element->FindParam(HydrologySheetBetaEnum);
	IssmDouble omega     = element->FindParam(HydrologyOmegaEnum);
	IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble mu_water  = element->FindParam(MaterialsMuWaterEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	Input* ev_input       = element->GetInput(HydrologyEnglacialVoidRatioEnum); _assert_(ev_input);
	Input* hr_input  = element->GetInput(HydrologyBumpHeightEnum);       _assert_(hr_input);
	Input* k_input   = element->GetInput(HydrologySheetConductivityEnum);_assert_(k_input);
	Input* phi_input = element->GetInput(HydraulicPotentialEnum);        _assert_(phi_input);
	Input* h_input   = element->GetInput(HydrologySheetThicknessEnum);   _assert_(h_input);
	Input* H_input   = element->GetInput(ThicknessEnum);                 _assert_(H_input);
	Input* b_input   = element->GetInput(BedEnum);                       _assert_(b_input);
	Input* B_input   = element->GetInput(HydrologyRheologyBBaseEnum);    _assert_(B_input);
	Input* n_input   = element->GetInput(MaterialsRheologyNEnum);        _assert_(n_input);	

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		phi_input->GetInputDerivativeValue(&dphi[0],xyz_list,gauss);
		phi_input->GetInputValue(&phi,gauss);
		h_input->GetInputValue(&h,gauss);
		k_input->GetInputValue(&k,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		ev_input->GetInputValue(&ev,gauss);
		hr_input->GetInputValue(&h_r,gauss);
		b_input->GetInputValue(&b,gauss);
		H_input->GetInputValue(&H,gauss);

		/*Get norm of gradient of hydraulic potential and make sure it is >0*/
		IssmDouble normgradphi = sqrt(dphi[0]*dphi[0] + dphi[1]*dphi[1]);
		if(normgradphi < AEPS) normgradphi = AEPS;

		/*Use transition model if specified*/
		IssmDouble nu = mu_water/rho_water;
		IssmDouble coeff;
		if(istransition==1 && omega>=AEPS){
			IssmDouble hratio = fabs(h/h_r);
			IssmDouble coarg = 1. + 4.*pow(hratio,3-2*alpha)*omega*k*pow(h,3)*normgradphi/nu;
			coeff = nu/2./omega*pow(hratio,2*alpha-3) * (-1 + pow(coarg, 0.5))/normgradphi;
		}
		else {
			/*If omega is zero, use standard model, otherwise transition model*/
			coeff = k*pow(h,alpha)*pow(normgradphi,beta-2.);
		}

		/*Diffusive term*/
		IssmDouble factor = gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += factor*(
							coeff*dbasis[0*numnodes+i]*dbasis[0*numnodes+j]
							+ coeff*dbasis[1*numnodes+i]*dbasis[1*numnodes+j]);
			}
		}

		/*Closing rate term, see Gagliardini and Werder 2018 eq. A2 (v = v1*phi_i + v2(phi_{i+1}))*/
		phi_0   = rho_water*g*b + rho_ice*g*H;
		if(isincludesheetthickness) phi_0 += rho_water*g*h;
		A=pow(B,-n);
		v1 = 2./pow(n,n)*A*h*(pow(fabs(phi_0 - phi),n-1.)*( - n));
		if (!creep_open_flag) {
			if (phi_0-phi<0) {
				v1 = 0;
			}
		}
		factor = gauss->weight*Jdet*(-v1);
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += factor*basis[i]*basis[j];
			}
		}

		/*Transient term if dt>0*/
		if(dt>0.){
			/*Diffusive term*/ 	
			factor = gauss->weight*Jdet*ev/(rho_water*g*dt);
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += factor*basis[i]*basis[j];
				}
			}
		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* HydrologyGlaDSAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         meltflag;
	IssmDouble  Jdet,w,v2,vx,vy,ub,h,h_r,ev;
	IssmDouble  G,m,melt,RO,frictionheat,alpha2;
	IssmDouble  A,B,n,phi_old,phi,phi_0;
	IssmDouble  H,b;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	bool isincludesheetthickness;
	bool creep_open_flag;
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	element->FindParam(&creep_open_flag,HydrologyCreepOpenFlagEnum);
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&meltflag,HydrologyMeltFlagEnum);
	IssmDouble L         = element->FindParam(MaterialsLatentheatEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble l_r       = element->FindParam(HydrologyCavitySpacingEnum);
	IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	Input* ev_input      = element->GetInput(HydrologyEnglacialVoidRatioEnum);           _assert_(ev_input);
	Input* hr_input     = element->GetInput(HydrologyBumpHeightEnum);                _assert_(hr_input);
	Input* h_input      = element->GetInput(HydrologySheetThicknessEnum);            _assert_(h_input);
	Input* H_input      = element->GetInput(ThicknessEnum);                          _assert_(H_input);
	Input* b_input      = element->GetInput(BedEnum);                                _assert_(b_input);
	Input* G_input      = element->GetInput(BasalforcingsGeothermalfluxEnum);        _assert_(G_input);
	Input* melt_input   = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);_assert_(melt_input);
	Input* RO_input     = NULL;
	Input* B_input      = element->GetInput(HydrologyRheologyBBaseEnum);             _assert_(B_input);
	Input* n_input      = element->GetInput(MaterialsRheologyNEnum);                 _assert_(n_input);
	Input* phiold_input = element->GetInput(HydraulicPotentialOldEnum);              _assert_(phiold_input);
	Input* phi_input    = element->GetInput(HydraulicPotentialEnum);                 _assert_(phi_input);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,2);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*Get input values at gauss points*/
		h_input->GetInputValue(&h,gauss);
		G_input->GetInputValue(&G,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		ev_input->GetInputValue(&ev,gauss);
		hr_input->GetInputValue(&h_r,gauss);
		phi_input->GetInputValue(&phi,gauss);
		b_input->GetInputValue(&b,gauss);
		H_input->GetInputValue(&H,gauss);
		melt_input->GetInputValue(&melt,gauss);

		/*Get basal velocity*/
		friction->GetBasalSlidingSpeeds(&vx, &vy ,gauss);
		ub = sqrt(vx*vx + vy*vy);

		/*Compute cavity opening w*/
		w  = 0.;
		if(h<h_r) w = ub*(h_r-h)/l_r;

		/*Compute frictional heat flux*/
		friction->GetAlpha2(&alpha2,gauss);
		frictionheat=alpha2*ub*ub;

		/*Compute melt (if necessary)*/
		if(meltflag == 0){
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

		/*Compute closing rate*/
		phi_0   = rho_water*g*b + rho_ice*g*H;
		if(isincludesheetthickness) phi_0 += rho_water*g*h;
		A=pow(B,-n);
		v2 = 2./pow(n,n)*A*h*(pow(fabs(phi_0 - phi),n-1.)*(phi_0 +(n-1.)*phi));
		if (!creep_open_flag) {
			if (phi_0-phi<0) {
				v2 = 0.;
			}
		}

		IssmDouble factor = - Jdet*gauss->weight*(w-v2-m);
		for(int i=0;i<numnodes;i++) pe->values[i]+= factor*basis[i];

		/*Transient term if dt>0*/
		if(dt>0.){
			phiold_input->GetInputValue(&phi_old,gauss);
			factor = gauss->weight*Jdet*ev/(rho_water*g*dt)*phi_old;
			for(int i=0;i<numnodes;i++) pe->values[i] += factor*basis[i];
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	delete friction;
	return pe;
}/*}}}*/
void           HydrologyGlaDSAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	/*if lakes are activated, updated phi at the lake outletboundary*/
	bool islakes;
    element->FindParam(&islakes,HydrologyLakeFlagEnum);
	/*Update element general phi*/
	element->GetSolutionFromInputsOneDof(solution,HydraulicPotentialEnum);
	/*Update lake outlet phi*/
	if(!islakes) return;
	else if(islakes){
		UpdateLakeOutletPhi(element);
	}
}/*}}}*/
void           HydrologyGlaDSAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyGlaDSAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	/*if lakes are activated, updated phi at the lake outletboundary*/
	bool islakes;
	element->FindParam(&islakes,HydrologyLakeFlagEnum);
	/*Update element general phi*/
	element->InputUpdateFromSolutionOneDof(solution,HydraulicPotentialEnum);
	/*Update lake outlet phi
	if(islakes){
		UpdateLakeOutletPhi(element);
	}*/
	/*Compute Hydrology Vx and Vy for time stepping purposes and Sheet Discharge as an optional output (These inputs do not affect GlaDS)*/

	/*Intermediaries*/
   	IssmDouble  dphi[3],h,k,phi;
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
	IssmDouble omega     = element->FindParam(HydrologyOmegaEnum);
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble mu_water  = element->FindParam(MaterialsMuWaterEnum);
	Input *k_input       = element->GetInput(HydrologySheetConductivityEnum); _assert_(k_input);
	Input *phi_input     = element->GetInput(HydraulicPotentialEnum);         _assert_(phi_input);
	Input *hr_input      = element->GetInput(HydrologyBumpHeightEnum);        _assert_(hr_input);
	Input *h_input       = element->GetInput(HydrologySheetThicknessEnum);    _assert_(h_input);
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
         	if(normgradphi < AEPS) normgradphi = AEPS;

         	/*If omega is zero, use standard model, otherwise transition model*/
         	IssmDouble nu = mu_water/rho_water;
			IssmDouble coeff;
			if(istransition==1 && omega>=AEPS){
				IssmDouble hratio = fabs(h/h_r);
				IssmDouble coarg = 1. + 4.*pow(hratio,3-2*alpha)*omega*k*pow(h,3)*normgradphi/nu;
				coeff = nu/2./omega*pow(hratio,2*alpha-3) * (-1 + pow(coarg, 0.5))/normgradphi;  // coeff gives discharge; divide by h to get speed instead of discharge
			}
			else {
			coeff = k*pow(h,alpha)*pow(normgradphi,beta-2.);  // coeff gives discharge; divide by h to get speed instead of discharge
			}

			vx[iv] = -coeff/max(AEPS,h)*dphi[0];
			vy[iv] = -coeff/max(AEPS,h)*dphi[1];

			d[iv] = coeff*normgradphi;
			}
		
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
void           HydrologyGlaDSAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

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

void HydrologyGlaDSAnalysis::UpdateLakeOutletPhi(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		UpdateLakeOutletPhi(element);
	}

}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateLakeOutletPhi(Element* element){/*{{{*/

	if(!element->IsAnyLake()){
			return;
	}
	else if(element->IsAnyLake()){
		/*Intermediaries*/
		IssmDouble bed, thickness, lh;
		IssmDouble LakeID;
		IssmDouble phi;

		/*Get number of vertices for this element*/
		int numvertices = element->GetNumberOfVertices();

		/*initialise phiLO*/
		IssmDouble *phiLO = xNew<IssmDouble>(numvertices);

		/*Retrieve all parameters*/
		IssmDouble rho_ice = element->FindParam(MaterialsRhoIceEnum);
		IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
		IssmDouble g = element->FindParam(ConstantsGEnum);
		Input *bed_input = element->GetInput(BaseEnum);_assert_(bed_input);
		Input *thickness_input = element->GetInput(ThicknessEnum);_assert_(thickness_input);
		Input *lh_input = element->GetInput(HydrologyLakeHeightEnum);_assert_(lh_input);
		Input *phi_input = element->GetInput(HydraulicPotentialEnum);_assert_(phi_input);
		Input *LakeID_input = element->GetInput(HydrologyLakeMaskEnum);_assert_(LakeID_input);

		/*Get values at gauss points*/
		Gauss *gauss = element->NewGauss();
		for (int in = 0; in < numvertices; in++){
			gauss->GaussVertex(in);

			/*Get input values at gauss points*/
			bed_input->GetInputValue(&bed, gauss);
			thickness_input->GetInputValue(&thickness, gauss);
			lh_input->GetInputValue(&lh, gauss);
			phi_input->GetInputValue(&phi, gauss);
			LakeID_input->GetInputValue(&LakeID, gauss);

			if (LakeID > 0.)
			{
				phiLO[in] = rho_water * g * bed + rho_water * g * lh;
			}
			else
			{
				phiLO[in] = phi;
			}
		}
		element->AddInput(HydraulicPotentialEnum, phiLO, P1Enum);
		/*clean up*/
		xDelete<IssmDouble>(phiLO);
		delete gauss;
	}		

}/*}}}*/

void HydrologyGlaDSAnalysis::SetChannelCrossSectionOld(FemModel* femmodel){/*{{{*/

	bool ischannels;
	femmodel->parameters->FindParam(&ischannels,HydrologyIschannelsEnum);
	if(!ischannels) return;

	for(int i=0;i<femmodel->loads->Size();i++){
		if(femmodel->loads->GetEnum(i)==ChannelEnum){
			Channel* channel=(Channel*)femmodel->loads->GetObjectByOffset(i);
			channel->SetChannelCrossSectionOld();
		}
	}

}/*}}}*/

void HydrologyGlaDSAnalysis::SetLakeOutletDischargeOld(FemModel* femmodel){/*{{{*/
	bool ischannels;
	femmodel->parameters->FindParam(&ischannels,HydrologyIschannelsEnum);
	bool islakes;
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	if(!ischannels || !islakes) return;
	else if(islakes){
		/*Initialize Qr vector*/
		int nv = femmodel->vertices->NumberOfVertices();
		Vector<IssmDouble>* Qr = new Vector<IssmDouble>(nv);

		/*Loop over all channels, and have them add their discharge*/
		for(int i=0;i<femmodel->loads->Size();i++){
			if(femmodel->loads->GetEnum(i)==ChannelEnum){
	   			Channel* channel=(Channel*)femmodel->loads->GetObjectByOffset(i);
	   		channel->AddDischargeToVector(Qr);
			}
 		}

 		/*Now assemble vector*/
 		Qr->Assemble();

 		/*Update inputs*/
  		InputUpdateFromVectorx(femmodel, Qr, HydrologyLakeOutletQrOldEnum, VertexSIdEnum);
  		delete Qr;
	}
}/*}}}*/

void HydrologyGlaDSAnalysis::UpdateLakeDepth(FemModel* femodel){/*{{{*/

	bool islakes;
	femodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	if(!islakes) return;
	else if(islakes){

		int numlakes,lakeID;
		femodel->parameters->FindParam(&numlakes,HydrologyNumLakesEnum);
		IssmDouble* dt = xNewZeroInit<IssmDouble>(numlakes+1);
		IssmDouble* lake_Qin = xNewZeroInit<IssmDouble>(numlakes+1);
		IssmDouble* lake_area = xNewZeroInit<IssmDouble>(numlakes+1);
		IssmDouble* lake_depth_old = xNewZeroInit<IssmDouble>(numlakes+1);
		IssmDouble* lake_depth = xNewZeroInit<IssmDouble>(numlakes+1);
		
		for (int lake=1;lake<=numlakes;lake++){

			IssmDouble local_qr = 0;
			IssmDouble total_qr;
			IssmDouble qs, lc, le;
			IssmDouble la,lh_old,qr,qin;
			IssmDouble A,B,alpha,beta;
			IssmDouble oceanLS,iceLS,LakeID;

			/*Add contribution of each sum element qr to lake total qr and get relevant values for lake area and Qin*/

			for (Object* & object : femodel->elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				
				/*fetch number of vertices*/
				IssmDouble connectivity;
				int numvertices = element->GetNumberOfVertices();
				/*Skip if inactive element*/
				if(element->IsAllFloating() || !element->IsIceInElement()){
					continue;
				}
				dt[lake]       = element->FindParam(TimesteppingTimeStepEnum);
				lc             = element->FindParam(HydrologyCavitySpacingEnum);
				int le         = element->CharacteristicLength();
				Input* qin_input = element->GetInput(HydrologyLakeQinEnum);           _assert_(qin_input); 
				Input* qr_input = element->GetInput(HydrologyLakeOutletQrEnum);    	  _assert_(qr_input);
				Input* qs_input = element->GetInput(HydrologySheetDischargeEnum);     _assert_(qs_input);
				Input* la_input   = element->GetInput(HydrologyLakeAreaEnum);		  _assert_(la_input);
				Input* lh_old_input = element->GetInput(HydrologyLakeHeightOldEnum);  _assert_(lh_old_input);
				Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum);      _assert_(oceanLS_input);
				Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum);          _assert_(iceLS_input);
				Input* LakeID_input = element->GetInput(HydrologyLakeMaskEnum);	  	  _assert_(LakeID_input);

				/* Start  looping on the number of gaussian points: */
				Gauss* gauss=element->NewGauss();
				for(int iv=0;iv<numvertices;iv++){
					gauss->GaussVertex(iv);

					/*Get input values at gauss points*/
					qs_input->GetInputValue(&qs,gauss);
					qin_input->GetInputValue(&qin,gauss);
					qr_input->GetInputValue(&qr,gauss);
					la_input->GetInputValue(&la,gauss);
					lh_old_input->GetInputValue(&lh_old,gauss);
						
					oceanLS_input->GetInputValue(&oceanLS,gauss);
					iceLS_input->GetInputValue(&iceLS,gauss);
					LakeID_input->GetInputValue(&LakeID,gauss);
					/*Set lake depth to zero if floating or no ice*/
					if(oceanLS<0. || LakeID==0. || iceLS>0.){
						continue;
					}

					if(LakeID==lake){
						/*divide nodal contributions to qr/qs by element connectivity to avoid issues with multiple counts of the same vertices*/
						connectivity=(IssmDouble)element->VertexConnectivity(iv);
						
						local_qr+=qr/connectivity;
						/*add approx sheet contribution*/
						
						if (qr>0.){
							local_qr += qs*(le)/connectivity;
						}
						else if (qr<0.){
							local_qr += -qs*(le)/connectivity;
						}
						lake_depth_old[lake] = lh_old;
						lake_Qin[lake] = qin;
						lake_area[lake] = la;
					}
					
				}
				/*clean up*/
				delete gauss;
			}
		
			/*collapse and recover qr into/out lakeN from everywhere else*/
			ISSM_MPI_Reduce(&local_qr,&total_qr,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
			ISSM_MPI_Bcast(&total_qr,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
			
			/*only solve for lake depth when lake_depth[lakes] (which is numlakes i+1 long) is lakes==1 */
			/*ODE solver for new lake height*/
			A = 1./lake_area[lake];
			B = lake_Qin[lake] - total_qr;

			if(lake_Qin[lake]==0 && lake_area[lake]==0){
				continue;
			}

			/*Define alpha and beta*/
			alpha = 0.;
			beta  = A*B;
			lake_depth[lake] = ODE1(alpha,beta,lake_depth_old[lake],dt[lake],2);
			/*Make sure it is positive*/
			if(lake_depth[lake]<AEPS) lake_depth[lake] = AEPS;
			
			/*store lake depth per vertices*/
			for (Object* & object : femodel->elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				/*fetch number of vertices*/
				int numvertices = element->GetNumberOfVertices();
				
				/*Skip if inactive element*/
				if(element->IsAllFloating() || !element->IsIceInElement()){
					continue;
				}
				/*Initialize new lake depth*/
				IssmDouble* lh_new = xNew<IssmDouble>(numvertices);
				/*Retrieve all inputs and parameters*/
				Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum);      _assert_(oceanLS_input);
				Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum);          _assert_(iceLS_input);
				Input* LakeID_input = element->GetInput(HydrologyLakeMaskEnum);	  _assert_(LakeID_input);
				/* Start  looping on the number of gaussian points: */
				Gauss* gauss=element->NewGauss();
				for(int iv=0;iv<numvertices;iv++){
					gauss->GaussVertex(iv);
					/*Get input values at gauss points*/
					oceanLS_input->GetInputValue(&oceanLS,gauss);
					iceLS_input->GetInputValue(&iceLS,gauss);
					LakeID_input->GetInputValue(&LakeID,gauss);
					/*Set lake depth to zero if floating or no ice*/
					if(oceanLS<0. || LakeID==0. || iceLS>0.){
						lh_new[iv] = 0.;
					}
					else if(LakeID==lake){
						lh_new[iv] = lake_depth[lake];
					}
				}
				/*Add new lake depth to element*/
				element->AddInput(HydrologyLakeHeightEnum,lh_new,P1Enum);
				/*clean up*/
				delete gauss;
				xDelete<IssmDouble>(lh_new);
				
			}

		}

		/*clean up*/
		xDelete<IssmDouble>(dt);
		xDelete<IssmDouble>(lake_Qin);
		xDelete<IssmDouble>(lake_area);
		xDelete<IssmDouble>(lake_depth_old);
		xDelete<IssmDouble>(lake_depth);
	}
	
}/*}}}*/

void HydrologyGlaDSAnalysis::UpdateSheetThickness(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		UpdateSheetThickness(element);
	}

}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateSheetThickness(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  vx,vy,ub,h_old,N,h_r,H,b;
	IssmDouble  A,B,n,phi,phi_0;
	IssmDouble  alpha,beta;
	IssmDouble  oceanLS,iceLS;

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize new sheet thickness*/
	IssmDouble* h_new = xNew<IssmDouble>(numvertices);

	/*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numvertices;iv++) {
			h_new[iv] = 0.;
		}
		element->AddInput(HydrologySheetThicknessEnum,h_new,P1Enum);
		xDelete<IssmDouble>(h_new);
		return;
	}

	/*Retrieve all inputs and parameters*/
	bool isincludesheetthickness;
	bool creep_open_flag;
	bool ishydrologyslope;
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	element->FindParam(&creep_open_flag,HydrologyCreepOpenFlagEnum);
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
	Input* hold_input = element->GetInput(HydrologySheetThicknessOldEnum);_assert_(hold_input);
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
		hold_input->GetInputValue(&h_old,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		hr_input->GetInputValue(&h_r,gauss);
		b_input->GetInputValue(&b,gauss);
		H_input->GetInputValue(&H,gauss);
		oceanLS_input->GetInputValue(&oceanLS,gauss);
		iceLS_input->GetInputValue(&iceLS,gauss);

		/*Set sheet thickness to zero if floating or no ice*/
		if(oceanLS<0. || iceLS>0.){
			h_new[iv] = 0.;
		}
		else{

		/*Get values for a few potentials*/
		phi_0   = rho_water*g*b + rho_ice*g*H;
		if(isincludesheetthickness) phi_0 += rho_water*g*h_old;
		N = phi_0 - phi;

		/*Get basal velocity*/
		ub = sqrt(vx*vx + vy*vy);

		/*Get A from B and n*/
		A = pow(B,-n);

		/*Define alpha and beta*/
		if(h_old<h_r){
			alpha = -ub/l_r - 2./pow(n,n)*A*pow(fabs(N),n-1.)*N;
			if (!creep_open_flag) {
				if (N<0) {
					alpha = -ub/l_r;
				}
			}
			beta  = ub*h_r/l_r;
		}
		else{
			alpha = - 2./pow(n,n)*A*pow(fabs(N),n-1.)*N;
			if (!creep_open_flag) {
				if (N<0) {
					alpha = 0.;
				}
			}
			beta  = 0.;
		}

		/*Get new sheet thickness*/
		h_new[iv] = ODE1(alpha,beta,h_old,dt,1);

		/*Make sure it is positive*/
		if(h_new[iv]<AEPS) h_new[iv] = AEPS;
		
		}

	}

	element->AddInput(HydrologySheetThicknessEnum,h_new,P1Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(h_new);
	delete gauss;
}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateEffectivePressure(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		UpdateEffectivePressure(element);
	}

}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateEffectivePressure(Element* element){/*{{{*/

	/*Intermediary*/
	IssmDouble phi_0, phi_m, p_i;
	IssmDouble H,b,phi,h;
	IssmDouble oceanLS,iceLS;

	int numnodes = element->GetNumberOfNodes();

	/*Get thickness and base on nodes to apply cap on water head*/
	bool isincludesheetthickness;
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	Input *h_input       = element->GetInput(HydrologySheetThicknessEnum);    _assert_(h_input);
   	IssmDouble* N = xNew<IssmDouble>(numnodes);
	IssmDouble  rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble  g         = element->FindParam(ConstantsGEnum);
	Input* H_input   = element->GetInput(ThicknessEnum); _assert_(H_input);
	Input* b_input   = element->GetInput(BedEnum); _assert_(b_input);
	Input* phi_input = element->GetInput(HydraulicPotentialEnum); _assert_(phi_input);
	Input* oceanLS_input = element->GetInput(MaskOceanLevelsetEnum); _assert_(oceanLS_input);
	Input* iceLS_input = element->GetInput(MaskIceLevelsetEnum); _assert_(iceLS_input);

	/*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numnodes;iv++) N[iv] = 0.;
		element->AddInput(EffectivePressureEnum,N,P1Enum);
		xDelete<IssmDouble>(N);
		return;
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(element->FiniteElement(),iv);

		/*Get input values at gauss points*/
		H_input->GetInputValue(&H,gauss);
		b_input->GetInputValue(&b,gauss);
		phi_input->GetInputValue(&phi,gauss);
		oceanLS_input->GetInputValue(&oceanLS,gauss);
		iceLS_input->GetInputValue(&iceLS,gauss);
		h_input->GetInputValue(&h,gauss);

		/*Elevation potential*/
		phi_m = rho_water*g*b;

		/*Compute overburden pressure*/
		p_i = rho_ice*g*H;

		/*Compute overburden potential*/
		phi_0 = phi_m + p_i;
		/*phi_0   = rho_water*g*b + rho_ice*g*H*/;
		if(isincludesheetthickness) phi_0 += rho_water*g*h;

		/*Calculate effective pressure*/
		N[iv] = phi_0 - phi;

		/*Make sure that all floating ice and ice free areas have zero effective pressure*/
		if(oceanLS<0.0) N[iv] = 0.0;
		if(iceLS>0.0) N[iv] = 0.0;

		if(xIsNan<IssmDouble>(N[iv])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(N[iv])) _error_("Inf found in solution vector");
	}

	element->AddInput(EffectivePressureEnum,N,element->FiniteElement());

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(N);
}/*}}}*/
void HydrologyGlaDSAnalysis::UpdateChannelCrossSection(FemModel* femmodel){/*{{{*/

	bool ischannels;
	femmodel->parameters->FindParam(&ischannels,HydrologyIschannelsEnum);
	if(!ischannels) return;

	for(int i=0;i<femmodel->loads->Size();i++){
		if(femmodel->loads->GetEnum(i)==ChannelEnum){
			Channel* channel=(Channel*)femmodel->loads->GetObjectByOffset(i);
			channel->UpdateChannelCrossSection();
		}
	}

}/*}}}*/

void HydrologyGlaDSAnalysis::UpdateLakeOutletDischarge(FemModel* femmodel){/*{{{*/
	bool ischannels;
	femmodel->parameters->FindParam(&ischannels,HydrologyIschannelsEnum);
	bool islakes;
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	if(!ischannels || !islakes){ 
		return;
	}
	else if(islakes){
		/*Initialize Qr vector*/
		int nv = femmodel->vertices->NumberOfVertices();
		Vector<IssmDouble>* Qr = new Vector<IssmDouble>(nv);

		/*Loop over all channels, and have them add their discharge*/
		for(int i=0;i<femmodel->loads->Size();i++){
			if(femmodel->loads->GetEnum(i)==ChannelEnum){
	   			Channel* channel=(Channel*)femmodel->loads->GetObjectByOffset(i);
	   		channel->AddDischargeToVector(Qr);
			}
 		}

 		/*Now assemble vector*/
 		Qr->Assemble();

 		/*Update inputs*/
  		InputUpdateFromVectorx(femmodel, Qr, HydrologyLakeOutletQrEnum, VertexSIdEnum);
  		delete Qr;
	}
}/*}}}*/
