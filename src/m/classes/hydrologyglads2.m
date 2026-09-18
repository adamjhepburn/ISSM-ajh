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
        
        %Other
        spch               	 = NaN;
        neumannflux          = NaN;
		relaxation_omega     = 0;
		englacial_void_ratio = 0.;
		requested_outputs    = {};
		melt_flag            = 0;
		stabilization         = 0; %0: no stabilization, 1: exact upwind

        %Channels
		ischannels           = 0;
		channel_conductivity = NaN;
		channel_sheet_width  = 0.;
		channel_alpha        = NaN; 
		channel_beta         = NaN; 

		%Other

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
            list = {'HydrologyWaterVx','HydrologyWaterVy','HydrologySheetHeight','HydrologySheetDischarge','HydrologyWaterPressure','HydraulicPotential','HydrologyFlowingSheetHeight','HydrologyMeanCavityHeight','ChannelArea','ChannelDischarge','ChannelWaterFilledArea'};
        end % }}}    

        function self = setdefaultparameters(self) % {{{

            %sheet parameters
            self.pressure_melt_coefficient = -7.5e-8; %K/Pa Clapeyron Slope (See table 2 in Wells et al)
            self.cavity_spacing = 2.; %m
            self.sheet_alpha = 5.0/4.0;
            self.sheet_beta = 3.0/2.0;
            self.bump_height = 0.1; %m
            self.sheet_conductivity = 1e-2;

			%channel parameters
			self.ischannels = 1;
			self.channel_conductivity = 1e-2;
			self.channel_sheet_width = 2.; %m
			self.channel_alpha = 5.0/4.0;
			self.channel_beta = 3.0/2.0;

            %other parameters
            self.englacial_void_ratio = 1e-4;
            self.melt_flag = 0;
			self.relaxation_omega = 0.5;
			self.stabilization = 1; %0: no stabilization, 1: exact upwind
            self.requested_outputs={'default'};

        end % }}}
        function md = checkconsistency(self,md,solution,analyses) % {{{
            %Early return
			if ~ismember('HydrologyGlaDS2Analysis',analyses)
				return;
			end

            %sheet
            md = checkfield(md,'fieldname','hydrology.pressure_melt_coefficient','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.sheet_conductivity','size',[md.mesh.numberofvertices 1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.cavity_spacing','numel',[1],'>',0);
			md = checkfield(md,'fieldname','hydrology.bump_height','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
			%md = checkfield(md,'fieldname','hydrology.omega', 'numel', [1], '>=', 0); 
			md = checkfield(md,'fieldname','hydrology.sheet_alpha', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.sheet_beta', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.rheology_B_base','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);

			%channel
			md = checkfield(md,'fieldname','hydrology.ischannels','numel',[1],'values',[0 1]);
			if self.ischannels==1
				md = checkfield(md,'fieldname','hydrology.channel_conductivity','size',[md.mesh.numberofvertices 1],'>',0,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','hydrology.channel_sheet_width','numel',[1],'>',0);
				md = checkfield(md,'fieldname','hydrology.channel_alpha', 'numel', [1], '>', 0); 
				md = checkfield(md,'fieldname','hydrology.channel_beta', 'numel', [1], '>', 0); 
			end

            %other
            md = checkfield(md,'fieldname','hydrology.spch','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.englacial_void_ratio','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.relaxation_omega','numel',[1],'>=',0,'<=',1);
			md = checkfield(md,'fieldname','hydrology.neumannflux','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','hydrology.stabilization','numel',[1],'>=',0,'<=',3);
			md = checkfield(md,'fieldname','hydrology.melt_flag','numel',[1],'values',[0 1 2]);
			if self.melt_flag==1 || self.melt_flag==2
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			end
        end % }}}

		function disp(self) % {{{
			disp(sprintf('   GlaDS2 (hydrologyglads2) solution parameters:'));
			disp(sprintf('      SHEET'));
			fielddisplay(self,'pressure_melt_coefficient','Pressure melt coefficient (c_t) [K Pa^-1]');
			fielddisplay(self,'sheet_conductivity','sheet conductivity (k) [m^(7/4) kg^(-1/2)]');
			fielddisplay(self,'sheet_alpha','First sheet-flow exponent (alpha_s) []'); 
			fielddisplay(self,'sheet_beta','Second sheet-flow exponent (beta_s) []'); 
			fielddisplay(self,'cavity_spacing','cavity spacing (l_r) [m]');
			fielddisplay(self,'bump_height','typical bump height (h_r) [m]');
			fielddisplay(self,'rheology_B_base','Ice rheology factor B at base of ice (B) [Pa s^(-1/3)]');

			disp(sprintf('      CHANNELS'));
			fielddisplay(self,'ischannels','Do we allow for channels? 1: yes, 0: no');
			fielddisplay(self,'channel_conductivity','channel conductivity (k_c) [m^(3/2) kg^(-1/2)]');
			fielddisplay(self,'channel_alpha','First channel-flow exponent (alpha_s) []'); 
			fielddisplay(self,'channel_beta','Second channel-flow exponent (beta_s) []'); 
			fielddisplay(self,'channel_sheet_width','channel sheet width [m]');
			disp(sprintf('      OTHER'));
			fielddisplay(self,'spch','Sheet thickness Dirichlet constraints [m]');
			fielddisplay(self,'neumannflux','water flux applied along the model boundary (m^2/s)');
			fielddisplay(self,'englacial_void_ratio','englacial void ratio (e_v)');
			fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'melt_flag','User specified basal melt? 0: no (default), 1: use md.basalforcings.groundedice_melting_rate');
			fielddisplay(self,'stabilization','Stabilization parameter for advection: 0: no stabilization, 1: exact upwind');
		end % }}}


        function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			%Marshall model code first
			WriteData(fid,prefix,'name','md.hydrology.model','data',8,'format','Integer');

						%Sheet
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','pressure_melt_coefficient','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_conductivity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','cavity_spacing','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_height','format','DoubleMat','mattype',1); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_alpha','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_beta','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','rheology_B_base','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ischannels','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_sheet_width','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_sheet_width','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_alpha','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_beta','format','Double'); 

			%channel 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ischannels','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_conductivity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_sheet_width','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_alpha','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_beta','format','Double');


			%Others
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','spch','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','neumannflux','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','relaxation_omega','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','englacial_void_ratio','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','melt_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','stabilization','format','Double');
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