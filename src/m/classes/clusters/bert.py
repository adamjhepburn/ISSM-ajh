import subprocess

import numpy as np

from fielddisplay import fielddisplay
try:
    from bert_settings import bert_settings
except ImportError:
    print('You need bert_settings.py to proceed, check presence and sys.path')
from helpers import *
from pairoptions import pairoptions
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from QueueRequirements import QueueRequirements


class bert(object):
    """BERT cluster class definition

    This is a SLURM queue on the BERT cluster (Aberystwyth University)

    Jobs can be submitted with various resource requirements:
    - numtasks: number of tasks to run
    - cpuspertask: number of CPUs per task
    - memorypernode: memory per node in GB

    Usage:
        cluster = bert()
        cluster = bert('numtasks', 4, 'cpuspertask', 8)
        cluster = bert('numtasks', 4, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = 'bert'
        self.login = ''
        self.numtasks = 1
        self.cpuspertask = 8
        self.port = 0
        self.projectaccount = ''
        self.codepath = ''
        self.executionpath = ''
        self.time = 24 * 60
        self.memorypernode = 2
        self.email = ''
        self.mailtype = ''
        
        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        self = bert_settings(self)

        # OK get other fields
        self = options.AssignObjectFields(self)
        self.np = self.numtasks * self.cpuspertask
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = "class bert object:"
        s = "%s\n%s" % (s, fielddisplay(self, 'name', 'name of the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'login', 'login'))
        s = "%s\n%s" % (s, fielddisplay(self, 'port', 'SSH port'))
        s = "%s\n%s" % (s, fielddisplay(self, 'numtasks', 'number of tasks'))
        s = "%s\n%s" % (s, fielddisplay(self, 'cpuspertask', 'number of CPUs per task'))
        s = "%s\n%s" % (s, fielddisplay(self, 'projectaccount', 'project account'))
        s = "%s\n%s" % (s, fielddisplay(self, 'codepath', 'code path on the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'executionpath', 'execution path on the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'time', 'walltime requested in minutes'))
        s = "%s\n%s" % (s, fielddisplay(self, 'memorypernode', 'memory per node in GB'))
        s = "%s\n%s" % (s, fielddisplay(self, 'email', 'email for job notifications'))
        s = "%s\n%s" % (s, fielddisplay(self, 'mailtype', 'mail type (BEGIN,END,FAIL,ALL)'))
        return s
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Check consistency of cluster configuration
        if not self.name:
            md = md.checkmessage('name empty')
        if not self.login:
            md = md.checkmessage('login empty')
        if not (self.numtasks > 0):
            md = md.checkmessage('numtasks must be > 0')
        if not (self.cpuspertask > 0):
            md = md.checkmessage('cpuspertask must be > 0')
        if not (self.port >= 0):
            md = md.checkmessage('port must be >= 0')
        if not self.codepath:
            md = md.checkmessage('codepath empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath empty')
        if not (self.time > 0):
            md = md.checkmessage('time must be > 0')
        if not (self.memorypernode > 0):
            md = md.checkmessage('memory must be > 0')
        return md
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        if isvalgrind:
            print('valgrind not supported by cluster, ignoring...')
        if isgprof:
            print('gprof not supported by cluster, ignoring...')

        # Write queuing script
        fid = open(modelname + '.queue', 'w')

        fid.write('#!/bin/bash\n')
        fid.write('#SBATCH --job-name=%s\n' % modelname)
        fid.write('#SBATCH --ntasks=%i\n' % self.numtasks)
        fid.write('#SBATCH --cpus-per-task=%i\n' % self.cpuspertask)
        fid.write('#SBATCH --time=%i\n' % self.time)  # walltime is in minutes
        fid.write('#SBATCH --mem=%igb\n' % self.memorypernode)  # memory in gigabytes
        fid.write('#SBATCH --output=%s.outlog\n' % modelname)
        fid.write('#SBATCH --error=%s.errlog\n' % modelname)
        fid.write('module load lmod mpich2/smpd/ge/gcc/64/1.3.2p1 autoconf/2.69 slurm/22.05.2-1\n\n')
        fid.write('export ISSM_DIR="%s/../"\n' % self.codepath)
        fid.write('source $ISSM_DIR/etc/environment.sh\n\n')
        fid.write('cd %s/%s\n\n' % (self.executionpath, dirname))
        fid.write('mpirun -n $SLURM_NTASKS %s/issm.exe %s %s %s\n' % (self.codepath, solution, self.executionpath + '/' + dirname, modelname))
        fid.write('sefffunction "$SLURM_JOB_ID"\n')
        
        if not io_gather:  # concatenate the output files
            fid.write('cat %s.outbin.* > %s.outbin\n' % (modelname, modelname))
        
        fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress the files into one tar.gz
        compressstring = 'tar -zcf %s.tar.gz ' % dirname
        for file in filelist:
            compressstring += ' {}'.format(file)
        subprocess.call(compressstring, shell=True)

        print('uploading input file and queuing script')
        print('Step 1: transferring files to central')
        
        central = 'central.aber.ac.uk'
        path = '/aber/%s' % self.login
        file = dirname + '.tar.gz'
        
        # Transfer to central server
        commandcentral = 'rsync -rav --progress %s %s@%s:%s' % (file, self.login, central, path)
        subprocess.call(commandcentral, shell=True)
        
        print('Step 2: transferring files from central to BERT')
        # Transfer from central to BERT
        commandbert = 'ssh -A %s@%s "rsync -rav --progress --remove-source-files %s/%s %s@%s:%s/%s"' % (
            self.login, central, path, file, self.login, self.name, self.executionpath, file)
        subprocess.call(commandbert, shell=True)
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        print('launching solution sequence on remote cluster')
        
        if not isempty(restart):
            launchcommand = 'cd %s && cd %s && hostname && sbatch %s.queue' % (
                self.executionpath, dirname, modelname)
        else:
            launchcommand = 'cd %s && rm -rf ./%s && mkdir %s && cd %s && mv ../%s.tar.gz ./ && tar -zxf %s.tar.gz && hostname && sbatch %s.queue' % (
                self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)
        
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        directory = '%s/%s/' % (self.executionpath, dirname)
        
        # Convert filelist to scp format
        fileliststr = ' '.join([directory + f for f in filelist])
        
        # Use scp to download files
        command = 'scp %s@%s:%s ./' % (self.login, self.name, fileliststr)
        subprocess.call(command, shell=True)
    # }}}
