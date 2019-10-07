function CSDTI_recon_on_cluster_scratch(scanner,runno,study,series,target_machine,fermi_filter,CS_recon_params,chunk_size)
%% Recent updates: 15 February 2017, BJ Anderson
%  Recompiled all dependent "work setup" executables with corrected
%  write_headfile.  Before this date, gradient table would suffer from
%  notable rounding errors.

%% Usage Notes (Added 23 November 2016, BJ Anderson)
%  scanner          :
%  runno            :
%  study            :
%  series           :
%  target_machine   :
%  fermi_filter     : Default (0). Can be turned on with default values of
%                     w1 = 0.15 and w2 = 0.75 with (1).  w1 and w2 can be
%                     specified with an underscore seperated (and quoted),
%                     string such as '1_0.15_0.75'
%                     ({fermi_filter_flag},w1,w2)
%  CS_recon_params  : This must be a string or otherwise will use default
%                     values.  Unlike fermi_filter, this quoted string must
%                     is comma-separated and case-sensitive, and has the format
%                     'variable_1={value_1},variable_2={value_2}' etc.  Any
%                     or all of the variables can be specified.  These are:
%                     'TVWeight' (default of 0.0012) - total variation
%                     'xfmWeight' (default of 0.006) - transform weight
%                     'Itnlim' (number of iterations, default of 48)
%  chunk_size       : Default(6).  The number of slices per cluster job.
%
%
%  BJ, ADD DESCRIPTION OF BASIC PROCESS OF DEPLOYING THE THREE MATLAB
%  EXECUTABLES



if ~exist('target_machine','var')
    target_machine = 'delos';
end

if ~exist('fermi_filter','var')
    fermi_filter=0;
end


if ~exist('chunk_size','var')
    chunk_size=6; % 25 November 2016, temporarily (?) changed to 6
end

if ~exist('CS_recon_params','var')
    CS_recon_params=''; % Added 23 December 2016 or otherwise default values won't work. --BJA
end


%%  cd /cm/shared/workstation_code_dev/shared/pipeline_utilities
addpath(genpath('/home/rmd22/Documents/MATLAB/'));

%% Get all necessary code for reconstruction
addpath('/cm/shared/workstation_code/shared/mathworks/slurm_shared');
addpath('/cm/shared/workstation_code/shared/civm_matlab_common_utils');
addpath('/cm/shared/workstation_code/shared/mathworks/dirrec');
addpath('/cm/shared/workstation_code/recon/mat_recon_pipe/mat_wrappers');
addpath('/cm/shared/workstation_code/shared/mathworks/dlmcell');
addpath('/home/rmd22/Documents/MATLAB/MATLAB_scripts_rmd/CS/');
addpath('/home/rmd22/Documents/MATLAB/MATLAB_scripts_rmd/CS/sparseMRI_v0.2');
addpath('/home/rmd22/Documents/MATLAB/MATLAB_scripts_rmd/CS/Wavelab850');
WavePath2;

%% Input checks
%Changed recon_path from input variable to hardcoded default for this code...will only need to change if there is amajor change to cluster configuration
%if ~exist('recon_path','var');
recon_path = fullfile('/glusterspace',[runno '.work']);
%end

if ~exist(recon_path,'dir');
    mkdir(recon_path);
end

% recon_path = [recon_path '/']; %this may not be necessary, just in case
if ~ischar(series)
    series = num2str(series);
    if length(series) < 2
        series = ['0' series];
    end
    series = ['ser' series];
end


%% Populate DTI series general header information
databuffer.engine_constants = load_engine_dependency();
databuffer.scanner_constants = load_scanner_dependency(scanner);
databuffer.headfile.U_runno = runno;
databuffer.headfile.U_scanner = scanner;
databuffer.input_headfile = struct;
optstruct.testmode = false;
optstruct.debug_mode = 0;
optstruct.warning_pause = 0;

optstruct.param_file =[runno '.param']; % Commented out on 18 Aug 2017,
%optstruct.param_file = boolean(0);%[runno '.param']; % Commented out on 18 Aug 2017,
%BJA
gui_info_collect(databuffer,optstruct);



%% 28 Aug 2017, BJA: Moved gui_info_collect related work to precede fid pulling, etc.

%% Pull fid and procpar, load reconstruction parameter data
% %CSreconfile = agilent2glusterspaceCS_wn(scanner,runno,study,series,recon_path);

reconfile = agilent2glusterspaceCS(scanner,runno,study,series,recon_path);
%reconfile=[recon_path '/' runno 'recon.mat'];
%fidpath=[recon_path '/' runno '.fid'];
%procpar_path_1=[recon_path '/' runno '.' ];
%procpar_path=[procpar_path_1 'procpar'];
%procpar = readprocparCS(procpar_path);
load(reconfile)
%% Pull CS table (bypass reading it from procpar directly)
full_CS_table_path = procpar.petableCS;
CS_table_parts = strsplit(full_CS_table_path{1},'/');
CS_table_name = CS_table_parts{end};

%target_folder = [fileparts(procparpath) '/'];
table_target = [recon_path '/' CS_table_name];

pull_table_cmd = ['puller_simple  -o -f file ' scanner ' ''../../../../home/vnmr1/vnmrsys/tablib/' CS_table_name ''' ' recon_path]
if ~exist(table_target,'file')
    system(pull_table_cmd)
end




%% Save new general header info to recon.mat file
m = matfile(reconfile,'Writable',true);
m.databuffer = databuffer;
m.optstruct = optstruct;

scale_file = [recon_path '/' runno '_4D_scaling_factor.float'];

%% Perform CS reconstruction and send data to naxosspace
%parpool('local',16) % Unsure if parpool will work...beware!
%par
for k = 1:nvols
    %CS_recon_cluster_bj_multithread(reconfile,k,workpath,scale_file,target_machine,chunk_size);
    CS_recon_cluster_bj_multithread_v2_scratch(reconfile,k,recon_path,scale_file,target_machine,fermi_filter,chunk_size,CS_recon_params); %v2
    %     im_res = CS_recon_cluster(reconfile,k,workpath);
end

%delete(gcp('nocreate'))
end