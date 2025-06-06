%issmSTANAGE class definition
%
%   Usage:
%      cluster=stanage();
%      cluster=stanage('np',3);
%      cluster=stanage('np',3,'login','username');

classdef stanage
    properties (SetAccess=public)  
		 % {{{
		 name           = ''
		 login          = '';
		 numtasks       = 1;
		 cpuspertask    = 8;
		 port           = 0;
		 projectaccount = '';
		 codepath       = '';
		 executionpath  = '';
		 time           = 24*60;
		 memorypernode         = 2;
		 email          = '';
		 mailtype       = '';
	 end
	 %}}}
	 methods
		 function cluster=stanage(varargin) % {{{

			 %initialize cluster using default settings if provided
			 if (exist('stanage_settings')==2), stanage_settings; end

			 %use provided options to change fields
			 cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		 end
		 %}}}
		 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
			 disp(sprintf('    port: %i',cluster.port));
			 disp(sprintf('    numtasks: %i',cluster.numtasks));
			 disp(sprintf('    cpuspertask: %i',cluster.cpuspertask));
			 %disp(sprintf('    projectaccount: %s',cluster.projectaccount));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    memory: %i',cluster.memorypernode));
			 %disp(sprintf('    email: %s', cluster.email));
			 %disp(sprintf('    mailtype: %s', cluster.mailtype));
			 
		 end
		 %}}}
		 function numprocs=np(cluster) % {{{
			 %compute number of processors
			 numprocs=cluster.numtasks*cluster.cpuspertask;
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{
			 if isempty(cluster.name), md = checkmessage(md,'name empty'); end
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			 if ~(cluster.numtasks > 0), md = checkmessage(md,'numtasks must be > 0'); end
			 if ~(cluster.cpuspertask > 0), md = checkmessage(md,'cpuspertask must be > 0'); end
			 if ~(cluster.port >= 0), md = checkmessage(md,'port must be >= 0'); end
			 %if isempty(cluster.email), md = checkmessage(md,'email empty'); end
			 %if isempty(cluster.mailtype), md = checkmessage(md,'mailtype empty'); end
			 %if isempty(cluster.projectaccount), md = checkmessage(md,'projectaccount empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
			 if ~(cluster.time > 0), md = checkmessage(md,'time must be > 0'); end
			 if ~(cluster.memorypernode > 0), md = checkmessage(md,'memory must be > 0'); end
		 end
		 %}}}
		 function BuildKrigingQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{
			 error('not implemented yet');
		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#!/bin/bash\n');
			 fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			 %fprintf(fid,'#SBATCH --account=%s \n',cluster.projectaccount);
             fprintf(fid,'#SBATCH --ntasks=%i  \n',cluster.numtasks);
			 fprintf(fid,'#SBATCH --cpus-per-task=%i\n',cluster.cpuspertask);
			 fprintf(fid,'#SBATCH --time=%i\n',cluster.time); %walltime is in minutes
			 fprintf(fid,'#SBATCH --mem=%igb\n',cluster.memorypernode); %memory in in gigabytes
			 %fprintf(fid,'#SBATCH --mail-user=%s\n',cluster.email); %email
			 %fprintf(fid,'#SBATCH --mail-type=%s\n',cluster.mailtype); 
			 fprintf(fid,'#SBATCH --output=%s.outlog \n',modelname);
			 fprintf(fid,'#SBATCH --error=%s.errlog \n',modelname);
             %fprintf(fid,'#SBATCH --export=ALL \n\n');
             fprintf(fid,'module purge\n');
             fprintf(fid,'module load intel CMake/3.24.3-GCCcore-12.2.0 Autotools\n\n');
			 %fprintf(fid,'echo $PATH\n\n');
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
             fprintf(fid,'source /users/gg4ahe/ISSM-JOKULLHAUP/etc/environment.sh\n\n');
             %fprintf(fid,'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n\n');
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			 fprintf(fid,'srun -n %i --export=ALL %s/issm.exe %s %s %s\n',cluster.np(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname); 
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);
		 end %}}}
		 function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{

			 %compress the files into one zip.
			 compressstring=['tar -zcf ' dirname '.tar.gz '];
			 for i=1:numel(filelist),
				 compressstring = [compressstring ' ' filelist{i}];
			 end
			 system(compressstring);

			 disp('uploading input file and queuing script');
			 issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		 end %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			 disp('launching solution sequence on remote cluster');
			 if ~isempty(restart)
				 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && sbatch ' modelname '.queue '];
			 else
				 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && sbatch ' modelname '.queue '];
			 end
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		 end %}}}
		 function Download(cluster,dirname,filelist)% {{{

			 %copy files from cluster to current directory
			 directory=[cluster.executionpath dirname '/'];
			 %issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
             fileliststr='{';
		     for i=1:numel(filelist)-1,
			    fileliststr=[fileliststr filelist{i} ','];
		    end
		    fileliststr=[fileliststr filelist{end} '}'];

            eval(['!scp ' cluster.login '@' cluster.name ':' directory fileliststr ' ./'])
		 end %}}}
	end
end
