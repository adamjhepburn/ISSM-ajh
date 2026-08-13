/*! \file HydrologyGlaDS2Analysis.h 
*  \brief: header file for generic external result object
*/

#ifndef _HydrologyGlaDS2Analysis_
#define _HydrologyGlaDS2Analysis_

/*Headers*/
#include "./Analysis.h"

class HydrologyGlaDS2Analysis: public Analysis{

    public:
        /*Model processing*/
        void CreateConstraints(Constraints* constraints,IoModel* iomodel);
        void CreateLoads(Loads* loads, IoModel* iomodel);
        void CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr=false);
		int  DofsPerNode(int** doflist,int domaintype,int approximation);
		void UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type);
		void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum);
        
        /*Finite element Analysis*/

        /*Specific to GlaDS2*/
};