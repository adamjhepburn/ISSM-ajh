%HYDROLOGYGLADS class definition
%
%   Usage:
%      hydrologyglads=hydrologyglads();

classdef hydrologyglads
	properties (SetAccess=public) 
		%Sheet
		pressure_melt_coefficient = 0.;
		sheet_conductivity        = NaN;
		cavity_spacing            = 0.;
		bump_height               = NaN;
		omega                     = 0; 
		sheet_alpha               = NaN; 
		sheet_beta                = NaN; 
		rheology_B_base           = NaN;
		isincludesheetthickness   = 0;
		creep_open_flag  		  = 1;
		elastic_sheet_flag        = 0;
		elastic_sheet_depth_scale = 0;
		elastic_sheet_exponent    = 0;
		uplift_reg_rate           = 0;
		reg_pressure              = 0;

		%Channels
		ischannels           = 0;
		channel_conductivity = NaN;
		channel_sheet_width  = 0.;
		channel_alpha        = NaN; 
		channel_beta         = NaN; 

		%Other
		spcphi               = NaN;
		moulin_input         = NaN;
		neumannflux          = NaN;
		englacial_void_ratio = 0.;
		requested_outputs    = {};
		melt_flag            = 0;
		islakes              = 0;
		lake_mask		     = 0;
		num_lakes		     = 0;
		lake_area			 = 0;
		lake_Qin			 = 0;
		istransition         = 0;
	end
	methods
		function self = hydrologyglads(varargin) % {{{
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
            if self.islakes
			    list = {'EffectivePressure','HydraulicPotential','HydrologySheetThickness','HydrologySheetDischarge','ChannelArea','ChannelDischarge','HydrologyLakeOutletQr','HydrologyLakeHeight'};
				if self.elastic_sheet_flag
					list = [list {'HydrologyElasticSheetThickness'}];
				end
			else
                list = {'EffectivePressure','HydraulicPotential','HydrologySheetThickness','ChannelArea','ChannelDischarge'};
				if self.elastic_sheet_flag
					list = [list {'HydrologyElasticSheetThickness'}];
				end
            end  
        end % }}}   
		function self = setdefaultparameters(self) % {{{

			%Sheet parameters
			self.pressure_melt_coefficient = 7.5e-8; %K/Pa (See table 1 in Erder et al. 2013)
			self.cavity_spacing = 2.; %m
			self.sheet_alpha = 5.0/4.0;
			self.sheet_beta = 3.0/2.0;
			self.omega = 1./2000.; 
			self.elastic_sheet_depth_scale = 0; %m, see git repo for Stevens et al., 2022
			self.elastic_sheet_exponent = 1; 
			self.uplift_reg_rate = 0.01/1e3/9.81; % m Pa^-1, ~1 m uplift for 100 m excess head
			self.reg_pressure = 1e4; %Pa, see git repo for Stevens et al., 2022

			%Channel parameters
			self.ischannels=false;
			self.channel_conductivity = 5.e-2; %Dow's default, Table uses 0.1
			self.channel_sheet_width = 2.; %m
			self.channel_alpha = 5.0/4.0;
			self.channel_beta = 3.0/2.0;

			%Otherself.omega = 1./2000.; 
			self.englacial_void_ratio = 1.e-5;% Dow's default, Table from Werder et al. uses 1e-3;
			self.requested_outputs={'default'};
			self.melt_flag=0;
			self.islakes=false; % by default no lakes at margin
			self.lake_mask = 0; %By default no lakes
			self.num_lakes = 0; %By default no lakes
			self.lake_area=0; % m2 the area of any ice-marginal lakes; 
			self.lake_Qin=0; %m3s-1, recharge rate of ice marginal lakes Kingslake and Ng 2013 use 1e-2;
			self.istransition = 0; %by default use GlaDS default turbulent code
			self.creep_open_flag = 1;
			self.elastic_sheet_flag = 0; % by default no elastic sheet
			
			
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyGladsAnalysis',analyses)
				return;
			end

			%Sheet
			md = checkfield(md,'fieldname','hydrology.pressure_melt_coefficient','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.sheet_conductivity','size',[md.mesh.numberofvertices 1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.cavity_spacing','numel',[1],'>',0);
			md = checkfield(md,'fieldname','hydrology.bump_height','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.omega', 'numel', [1], '>=', 0); 
			md = checkfield(md,'fieldname','hydrology.sheet_alpha', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.sheet_beta', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.rheology_B_base','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.isincludesheetthickness','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.elastic_sheet_flag','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.elastic_sheet_depth_scale','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.elastic_sheet_exponent','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.uplift_reg_rate','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.reg_pressure','numel',[1],'>=',0);
			%Channels
			md = checkfield(md,'fieldname','hydrology.ischannels','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.channel_conductivity','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.channel_sheet_width','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','hydrology.channel_alpha', 'numel', [1], '>', 0); 
			md = checkfield(md,'fieldname','hydrology.channel_beta', 'numel', [1], '>', 0); 

			%Other
			md = checkfield(md,'fieldname','hydrology.spcphi','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.englacial_void_ratio','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.moulin_input','>=',0,'timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.neumannflux','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','hydrology.melt_flag','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','hydrology.islakes','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.lake_mask','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.num_lakes','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.istransition','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.creep_open_flag','numel',[1],'values',[0 1]);
			if self.melt_flag==1 || self.melt_flag==2
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			end
			if self.islakes==1
				md = checkfield(md,'fieldname','hydrology.lake_mask','Inf',1,'NaN',1,'timeseries',1);
				md = checkfield(md,'fieldname','hydrology.lake_area','size',[md.mesh.numberofvertices 1],'>=',0,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','hydrology.lake_Qin','timeseries',1,'>=',0,'NaN',1,'Inf',1);
            elseif self.islakes==0
                md = checkfield(md,'fieldname','hydrology.lake_mask','size',[md.mesh.numberofvertices 1]);
            end
 
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   GlaDS (hydrologyglads) solution parameters:'));
			disp(sprintf('      SHEET'));
			fielddisplay(self,'pressure_melt_coefficient','Pressure melt coefficient (c_t) [K Pa^-1]');
			fielddisplay(self,'sheet_conductivity','sheet conductivity (k) [m^(7/4) kg^(-1/2)]');
			fielddisplay(self,'sheet_alpha','First sheet-flow exponent (alpha_s) []'); 
			fielddisplay(self,'sheet_beta','Second sheet-flow exponent (beta_s) []'); 
			fielddisplay(self,'cavity_spacing','cavity spacing (l_r) [m]');
			fielddisplay(self,'bump_height','typical bump height (h_r) [m]');
			fielddisplay(self,'omega','transition parameter (omega) []'); 
			fielddisplay(self,'rheology_B_base','Ice rheology factor B at base of ice (B) [Pa s^(-1/3)]');
			fielddisplay(self,'isincludesheetthickness','Do we add rho_w*g*h in effective pressure calculation? 1: yes, 0: no');
			fielddisplay(self,'creep_open_flag','Do we allow cavities to open by creep when N<0? 1: yes, 0: no');
			fielddisplay(self,'elastic_sheet_flag','Does sheet thickness include we an elastic sheet? 1: yes, 0: no');
			fielddisplay(self,'elastic_sheet_depth_scale','Elastic sheet depth scale (c_e) [m]');
			fielddisplay(self,'elastic_sheet_exponent','Elastic sheet exponent (\gamma) []');
			fielddisplay(self,'uplift_reg_rate','Uplift regularization rate (h_{\varepsilon}) [m Pa^-1]');
			fielddisplay(self,'reg_pressure','Regularizing pressure for uplift regularisation (N_{\varepsilon}) [Pa]');

			disp(sprintf('      CHANNELS'));
			fielddisplay(self,'ischannels','Do we allow for channels? 1: yes, 0: no');
			fielddisplay(self,'channel_conductivity','channel conductivity (k_c) [m^(3/2) kg^(-1/2)]');
			fielddisplay(self,'channel_alpha','First channel-flow exponent (alpha_s) []'); 
			fielddisplay(self,'channel_beta','Second channel-flow exponent (beta_s) []'); 
			fielddisplay(self,'channel_sheet_width','channel sheet width [m]');
			disp(sprintf('      OTHER'));
			fielddisplay(self,'spcphi','Hydraulic potential Dirichlet constraints [Pa]');
			fielddisplay(self,'neumannflux','water flux applied along the model boundary (m^2/s)');
			fielddisplay(self,'moulin_input','moulin input (Q_s) [m^3/s]');
			fielddisplay(self,'englacial_void_ratio','englacial void ratio (e_v)');
			fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'melt_flag','User specified basal melt? 0: no (default), 1: use md.basalforcings.groundedice_melting_rate');
			fielddisplay(self,'islakes','User specified lake? 0: no (default), 1: use md.hydrology.lake_mask to identify lake outlets');
			fielddisplay(self,'lake_mask','lake mask (0: for no lake, 1,2,...n for n lakes)');
			fielddisplay(self,'num_lakes','Number of lakes (0: no lakes, 1: one lake, ... n: n lakes)');
			fielddisplay(self,'lake_area','Lake area at vertex (Qr) [m^2]');
			fielddisplay(self,'lake_Qin','Lake refill rate (Qin) [m^3/s]');
			fielddisplay(self,'istransition','do we use standard [0, default] or transition model [1]');
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
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','omega','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_alpha','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_beta','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','rheology_B_base','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','isincludesheetthickness','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','creep_open_flag','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','elastic_sheet_flag','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','elastic_sheet_depth_scale','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','elastic_sheet_exponent','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','uplift_reg_rate','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','reg_pressure','format','Double');

			%Channels
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ischannels','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_conductivity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_sheet_width','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_alpha','format','Double'); 
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_beta','format','Double'); 

			%Others
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','spcphi','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','neumannflux','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','moulin_input','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','englacial_void_ratio','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','melt_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','islakes','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','lake_mask','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','num_lakes','format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','lake_area','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','lake_Qin','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','istransition','format','Boolean');
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];  %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
	end
end
