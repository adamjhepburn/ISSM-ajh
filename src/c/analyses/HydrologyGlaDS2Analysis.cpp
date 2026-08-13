#include <float.h> /* defines DBL_EPSILON*/
#include "./HydrologyGlaDSAnalysis.h"
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

}/*}}}*/