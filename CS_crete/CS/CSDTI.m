function [CS_recon_job,b0_job] = CSDTI(scanner,runno,study,series,recon_path)

%% Get all necessary code for reconstruction
addpath('/cm/shared/workstation_code/shared/mathworks/slurm_shared');
addpath('/cm/shared/workstation_code/shared/civm_matlab_common_utils');
addpath('/cm/shared/workstation_code/shared/mathworks/dirrec');
addpath('/cm/shared/workstation_code/recon/mat_recon_pipe/mat_wrappers');
addpath('/cm/shared/workstation_code/shared/mathworks/dlmcell');
WavePath2;
% addpath(genpath('/home/rmd22/Documents/MATLAB/'));

%% Input checks
if ~exist('recon_path','var');
    recon_path = fullfile('/glusterspace',runno);
end
mkdir(recon_path);
% recon_path = [recon_path '/']; %this may not be necessary, just in case
if ~ischar(series)
    series = num2str(series);
    if length(series) < 2
        series = ['0' series];
    end
    series = ['ser' series];
end

%% Pull fid and procpar, load reconstruction parameter data
reconfile = agilent2glusterspaceCS(scanner,runno,study,series,recon_path);
load(reconfile)

%% Populate DTI series general header information
databuffer.engine_constants = load_engine_dependency();
databuffer.scanner_constants = load_scanner_dependency(scanner);
databuffer.headfile.U_runno = runno;
databuffer.headfile.U_scanner = scanner;
databuffer.input_headfile = struct;
optstruct.testmode = false;
optstruct.debug_mode = 0;
optstruct.warning_pause = 0;
optstruct.param_file = [runno '.param'];
gui_info_collect(databuffer,optstruct);

%% Save new general header info to recon.mat file
m = matfile(reconfile,'Writable',true);
m.databuffer = databuffer;
m.optstruct = optstruct;

%% Send first image volume to cluster node (should be a b0 image to compute group scaling)
k = 1; 
if ~exist('cc1s','var')
    cc1s=parcluster('civmcluster1_slurm');
end
b0_job=createJob(cc1s);
b0_job.createTask(@CS_recon_cluster,0,{reconfile,k,workpath});
b0_job.submit()
% b0_job.wait()
% b0_job.delete()

% Wait until scaling is in the recon.mat file (I'm sure there's a better
% way to code this...)
mycount = 0;
while mycount < 60*5 % time out after 5 minutes
    try
        fprintf('\nScaling calculated: %f\n',m.scaling);
        break
    catch
        pause(1);
        mycount = mycount+1;
        fprintf('%d ',mycount)
        if mod(mycount,30) == 0
            fprintf('\n');
        end
    end
end

%% Perform CS reconstruction and send data to naxosspace
CS_recon_job=createJob(cc1s);
for k = 2:nvols
    CS_recon_job.createTask(@CS_recon_cluster,0,{reconfile,k,workpath});
%     im_res = CS_recon_cluster(reconfile,k,workpath);
end
CS_recon_job.submit();

% fprintf('Reconstruction of %d image volumes sent to cluster nodes.\n',nvols);
% fprintf('Outputs will be sent to /Volumes/naxosspace/%s_m*/\n',runno);
% fprintf('The fid, procpar, and recon parameters will be sent to /Volumes/naxosspace/%s.work\n',runno);
% pw = '4.signa!';
% workdir = ['/Volumnaxoes/naxosspace/' runno '.work'];
% system(['sshpass -p ' pw ' ssh omega@naxos.duhs.duke.edu ''mkdir -p ' workdir '''']);
% copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(recon_path,[runno '.fid']) ...
%             ' omega@naxos.duhs.duke.edu:' fullfile(workdir,[runno '.fid'])];
% system(copy_cmd);       
% copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(recon_path,[runno '.procpar']) ...
%     ' omega@naxos.duhs.duke.edu:' fullfile(workdir,[runno '.procpar'])];
% system(copy_cmd);
% copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(recon_path,[runno 'recon.mat']) ...
%     ' omega@naxos.duhs.duke.edu:' fullfile(workdir,[runno 'recon.mat'])];
% system(copy_cmd);
    
fprintf('Reconstruction of %d image volumes sent to cluster nodes.\n',nvols);
fprintf('Outputs will be sent to /Volumes/delosspace/%s_m*/\n',runno);
fprintf('The fid, procpar, and recon parameters will be sent to /Volumes/delosspace/%s.work\n',runno);
pw = '4.signa!';
workdir = ['/Volumes/delosspace/' runno '.work'];
system(['sshpass -p ' pw ' ssh omega@delos.duhs.duke.edu ''mkdir -p ' workdir '''']);
copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(recon_path,[runno '.fid']) ...
            ' omega@delos.duhs.duke.edu:' fullfile(workdir,[runno '.fid'])];
system(copy_cmd);       
copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(recon_path,[runno '.procpar']) ...
    ' omega@delos.duhs.duke.edu:' fullfile(workdir,[runno '.procpar'])];
system(copy_cmd);
copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(recon_path,[runno 'recon.mat']) ...
    ' omega@delos.duhs.duke.edu:' fullfile(workdir,[runno 'recon.mat'])];
system(copy_cmd);


end