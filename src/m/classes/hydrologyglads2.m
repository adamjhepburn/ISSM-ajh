%HYDROLOGYGLADS-2 class definition
%
%   Usage:
%      hydrologyglads2=hydrologyglads2();

classdef hydrologyglads2
    properties (SetAccess=public) 
        %Sheet
        pressure_melt_coefficient = 0.;
        sheet_conductivity        = NaN;
        cavity_spacing            = 0.;
        bump_height               = NaN;
        %omega                     = 0; 
        sheet_alpha               = NaN; 
        sheet_beta                = NaN; 
        rheology_B_base           = NaN;
        
        %Channels
        ischannels           = 0;
        channel_conductivity = NaN;
        channel_sheet_width  = 0.;
        channel_alpha        = NaN; 
        channel_beta         = NaN; 

        %Other
        spch               = NaN;
        neumannflux          = NaN;


        %Channels
		%ischannels           = 0;
		%channel_conductivity = NaN;
		%channel_sheet_width  = 0.;
		%channel_alpha        = NaN; 
		%channel_beta         = NaN; 

		%Other
		spcphi               = NaN;
		moulin_input         = NaN;
		neumannflux          = NaN;
		englacial_void_ratio = 0.;
		requested_outputs    = {};
		melt_flag            = 0;
		%istransition         = 0;
	end
    methods
        function self = hydrologyglads2(varargin) % {{{
            switch nargin
                case 0
                    self=setdefaultparameters(self);
                case 1
                    self=structtoobj(self,varargin{1});
                otherwise
                    error('constructor not supported');
            end
        end % }}}
        function list = defaultoutputs(self,md) % {{{
            list = {'HydrologyWaterVx','HydrologyWaterVy','HydrologySheetDischarge','HydrologyWaterPressure','HydraulicPotential','HydrologyFlowingSheetHeight','HydrologyMeanCavityHeight'};
        end % }}}    

        function self = setdefaultparameters(self) % {{{

            %sheet parameters
            self.pressure_melt_coefficient = 7.5e-8; %K/Pa Clapeyron Slope (See table 2 in Wells et al)
            self.cavity_spacing = 2.; %m
            self.sheet_alpha = 5.0/4.0;
            self.sheet_beta = 3.0/2.0;
            self.bump_height = 0.1; %m
            self.sheet_conductivity = 1e-2;

            %other parameters
            self.englacial_void_ratio = 1e-4;
            self.melt_flag = 0;
            self.requested_outputs={'default'};

        end % }}}
        function md = checkconsistency(self,md,solution,analyses) % {{{
            %Early return
			if ~ismember('HydrologyGlads2Analysis',analyses)
				return;
			end

            %sheet
            md = checkfield(md,'fieldname','hydrology.pressure_melt_coefficient','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.sheet_conductivity','size',[md.mesh.numberofvertices 1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.cavity_spacing','numel',[1],'>',0);
			md = checkfield(md,'fieldname','hydrology.bump_height','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
			%md = checkfield(md,'fieldname','hydrology.omega', 'numel', [1], '>=', 0); 
			md = checkfield(md,'fieldname','hydrology.sheet_alpha', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.sheet_beta', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.rheology_B_base','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);

            %other
            md = checkfield(md,'fieldname','hydrology.spch','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.englacial_void_ratio','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.neumannflux','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','hydrology.melt_flag','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','hydrology.istransition','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.creep_open_flag','numel',[1],'values',[0 1]);
			if self.melt_flag==1 || self.melt_flag==2
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			end
        end % }}}

        function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			%Marshall model code first
			WriteData(fid,prefix,'name','md.hydrology.model','data',5,'format','Integer');

						%Sheet
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','pressure_melt_coefficient','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_conductivity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','cavity_spacing','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_height','format','DoubleMat','mattype',1); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_alpha','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_beta','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','rheology_B_base','format','DoubleMat','mattype',1);


			%Others
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','spch','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','neumannflux','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','englacial_void_ratio','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','melt_flag','format','Integer');
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos)
				outputs(pos) = [];  %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}



    end
end