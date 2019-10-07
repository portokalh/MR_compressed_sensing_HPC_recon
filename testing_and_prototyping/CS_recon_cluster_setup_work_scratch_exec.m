function [ ] = CS_recon_cluster_setup_work_scratch_exec(variables_file,volume_number)
%CS_RECON_CLUSTER_SETUP_WORK_EXEC An executable MATLAB script for setting
%up each volume of CS reconstruction in order to avoid saturating the
%master node (in the context of DTI, with many many volumes to recon.

%%   Import Variables
load(variables_file);
reconfile = variables.reconfile;
procpar_path = variables.procparpath;
outpath = variables.outpath;
scale_file = variables.scale_file;
target_machine = variables.target_machine;

if isfield(variables,'TVWeight') %exist('variables.TVWeight','var') Does not work!!!
    TVWeight = variables.TVWeight;
else 
    TVWeight = 0.0012;
end

if isfield(variables,'xfmWeight') % exist('variables.xfmWeight','var') Does not work!!!
    xfmWeight = variables.xfmWeight;
else
    xfmWeight = 0.006;
end

if isfield(variables,'Itnlim') %exist('variables.Itnlim','var') Does not work!!!
    Itnlim = variables.Itnlim;
else
    Itnlim = 48;
end


if isfield(variables,'sampling_fraction') % exist('variables.sampling_fraction','var') Does not work!!!
    sampling_fraction = variables.sampling_fraction;
end
wavelet_dims = variables.wavelet_dims;

if isfield(variables,'wavelet_type') 
    wavelet_type= variables.wavelet_type;
else
    wavelet_type = 'Daubechies';
end


%% Load data
load(reconfile);
mask = skipint2skiptable(procpar_path); %sampling mask
mask0 = mask; %% BJ

volume_number=str2double(volume_number);


%% Immediately check to see if volume has already been reconstructed
if  ~exist( 'do_russway','var')
    myvolstr =sprintf(['%0' num2str(numel(num2str(nvols-1))) 'i' ],volume_number-1);
else
    myvolstr = num2str(volume_number-1);
    if volume_number-1 < 10
        myvolstr = ['00' myvolstr];
    elseif volume_number-1 < 100
        myvolstr = ['0' myvolstr];
    end
end

dir1 = [runno '_m' myvolstr];
voldir = fullfile(outpath,dir1,[dir1 'images']);
if ~exist(voldir,'dir')
    mkdir(voldir);
end


volume_message = ['Current volume = ' myvolstr '.']

%disp(volume_message)
fidpath=[outpath '/' runno '_m' myvolstr '.fid'] %%NEW
temp_volume_number=1; %%NEW
data = double(load_fidCS(fidpath,max_blocks,ntraces,npoints,bitdepth,temp_volume_number,voldims));
data0=data;

dims = size(data);

% Get sampling PDF (this is not the sampling mask)
% my_sampling_fraction = 0.5;
if ~exist('sampling_fraction','var')
    sampling_fraction = sum(mask(:))/numel(mask); % estimate if no sampling fraction provided
end

petableCS = procpar.petableCS;
pieces = strsplit(petableCS{1},'_pa');
pa_pb = strsplit(pieces{2},'_pb');
pa=str2double(pa_pb{1})/10;
pb=str2double(pa_pb{2})/10;
%pa = 1.8;
%pb = 3.6;

[mypdf,~] = genPDF_wn_v2(dims(2:3),pa,sampling_fraction,pb,false);

mypdf0 = mypdf;

%
% % pad if non-square or non-power of 2

% % padarry to 512*512 for better recon?
dyadic_idx = 2.^(1:12); %dyadic_idx = 2.^[1:12];
pidx = find(max(dims(2:3))<=dyadic_idx,1);
p = 2^pidx;

data = padarray(data,[0 p-dims(2) p-dims(3)]/2,0,'both');
mask = padarray(mask,[p-voldims(2) p-voldims(3)]/2,0,'both');
mypdf = padarray(mypdf,[p-dims(2) p-dims(3)]/2,1,'both'); %pad with 1's since we don't want to divide by zero later

dims1 = size(data);
% Calculate group scaling from first b0 image
if (~exist('scale_file','file') && (volume_number==1))
    %if ~exist('scaling','var')
    
    data0 = fftshift(ifft(fftshift(data0,1),[],1),1); % take ifft in the fully sampled dimension
    data0 = permute(data0,[2,3,1]);
    for n = 1:dims(1)
        data0(:,:,n) = fftshift(ifftn(fftshift(data0(:,:,n)./mypdf0))); % this compensates the intensity for the undersampling
    end
    
    thresh = .999999999; % only clip out very noisy voxels (e.g. zippers, artifacts), and not the eyes
    
    q = quantile(abs(data0(:)),thresh);
    scaling = (2^16-1)/q; % we plan on writing out uint16 not int16, though it won't show up well in the database
    m = matfile(reconfile,'Writable',true);
    m.scaling = scaling;
    
    % Write scaling factor to scale file
    fid = fopen(scale_file,'w');
    fwrite(fid,scaling,'float');
    fclose(fid);
    
    %%% Save the CS sampling mask in the recon.mat as well.
    % This will need to be changed once we start using different masks for
    % different image volume acquisitions. Since we are stuck using just
    % one mask for all acquisitions, it doesn't hurt to store it for later
    % reference. The unconverted mask (aka skiptable), still in Agilent
    % format, can be found in the procpar file
    m.CSmask = mask0;
else
    count=0;
    while (~exist(scale_file,'file'))
        if count > 40
            error_message = 'Waited more than 10 minutes for scaling file to be created; aborting.';
            error(error_message);
        else
            pause(15)
            count=count+1;
        end
    end
    fid=fopen(scale_file,'r');
    scaling = fread(fid,inf,'*float');
end

%% Prep data for reconstruction

%%%%%% One could potentially zero-pad non-square data here in order to use
%%%%%% the Wavelet transform version of CS.

% Convert to a convenient format for CS recon
data = fftshift(ifft(fftshift(data,1),[],1),1); % take ifft in the fully sampled dimension
data = permute(data,[2,3,1]);

% Calculate scaling
im_res = zeros(size(data),'like',data);
im_zfwdc = im_res;
for n=1:dims(1)
    % ifft2c (Lustig) function multiplies by sqrt(length(x(:))). Scaling is
    % required for the CS recon to work. This scaling is later undone using
    % sqrt(numel(mask)) near the end of this script.
    im_zfwdc(:,:,n) = ifft2c(data(:,:,n)./mypdf); % this compensates the intensity for the undersampling
end

% scale data such that the maximum image pixel in zf-w/dc is around 1
% this way, we can use similar lambda for different problems
myscale = max(abs(im_zfwdc(:)));
data = data/myscale;


% Define recon parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = dims1(2:3); 		% image Size
DN = N;		 	% data Size
% TVWeight = 0.0012; 	% Weight for TV penalty - only TV is on, but I encourage you to try wavelets as well.
% xfmWeight = 0.006;	% Weight for Transform L1 penalty
% Itnlim = 48;	%10	% Number of iterations 
OuterIt = length(TVWeight);
% 
% if (exist('CS_recon_params','var') && ischar(CS_recon_params))
%     CS_recon_params_2 = strjoin(strsplit(CS_recon_params,' '),'');
%     CS_recon_opts = strsplit(CS_recon_params_2,',');
%     for ii=1:length(CS_recon_opts)
%         var_and_val = strsplit(CS_recon_opts{ii},'=');
%         if length(var_and_val) == 2
%             variable_i = var_and_val{1};
%             value_i = var_and_val{2};
%             value_ii = str2double(value_i);
%             switch variable_i
%                 case 'TVWeight'
%                     if ((value_ii > 0) && (value_ii <= 1))
%                         TVWeight = value_ii;
%                     end
%                 case 'xfmWeight'
%                     if ((value_ii > 0) && (value_ii <= 1))
%                         xfmWeight = value_ii;
%                     end
%                 case 'Itnlim'
%                     if (value_ii) % Zero iterations not allowed!
%                         Itnlim = ceil(abs(value_ii));
%                     end
%                 otherwise
%             end
%         end
%     end
%     msg = ['Custom CS recon parameter(s) specified: TVWeight=' num2str(TVWeight) ...
%         ', xfmWeight=' num2str(xfmWeight) ', and Itnlim=' num2str(Itnlim) '.'];
%     
%     disp(msg);
%     disp('Please check that these values are correct.');
%     disp(['(Your input was: ''' CS_recon_params '''.']);
% end


phmask = zpad(hamming(64)*hamming(64)',N(1),N(2)); %mask to grab center frequency
phmask = phmask/max(phmask(:));			 %for low-order phase estimation and correction

%% Make and save header data

% DTI volume-specific header data
gui_info = read_headfile(fullfile(databuffer.engine_constants.engine_recongui_paramfile_directory,[runno '.param']));
agi_struct = combine_struct(struct,procpar,'z_Agilent_');
struct1 = combine_struct(agi_struct,gui_info,'U_');
struct1.state = struct1.U_state;
struct1.U_runno = dir1;
struct1 = combine_struct(struct1,databuffer.scanner_constants);
struct1.bw = procpar.sw/2;
struct1.A_channels = 1;
struct1.A_dti_vols = nvols;
struct1.A_echoes = nechoes;

if isfield('procpar','max_bval')
    struct1.A_max_bval = procpar.max_bval; 
end

struct1.B_recon_type = 'CS_recon_cluster';
struct1.U_stored_file_format = 'raw';

struct1.CS_TVWeight = TVWeight;
struct1.CS_xfmWeight = xfmWeight;
struct1.CS_Itnlim = Itnlim;
struct1.CS_OuterIt = OuterIt;

struct1.S_PSDname = char(procpar.seqfil);
struct1.S_scanner_tag = 'A_';
struct1.S_tag = 'A_';
struct1.S_tesla = struct1.scanner_tesla_image_code;
struct1.alpha = procpar.flip1;
struct1.dim_X = dims(1);
struct1.dim_Y = dims(2);
struct1.dim_Z = dims(3);
struct1.fovx = fov(1);
struct1.fovy = fov(2);
struct1.fovz = fov(3);
struct1.hfpmcnt = 1; % not sure where this comes from or if it needs to change
struct1.matlab_functioncall = 'enter_later'; % enter this later
struct1.ne = nechoes;
struct1.scanner_vendor = 'agilent';
struct1.te = procpar.te*1e3;
struct1.tr = procpar.tr*1e6;

struct1.work_dir_path = ['/' target_machine 'space/' runno '.work'];
struct1.CS_sampling_fraction = sampling_fraction;
struct1.CS_acceleration = 1/sampling_fraction;
struct1.group_image_scale = scaling;

write_headfile(fullfile(voldir,[dir1 '.headfile']),struct1,'',0);

tic

work_folder = [outpath '/' dir1 '/work/']; % Will need to delete after work is completed

if ~exist(work_folder,'dir')
    dir_cmd2 = ['mkdir ' work_folder '; chmod 777 ' work_folder];
    system(dir_cmd2);
else
    chmod_cmd = ['chmod 777 ' work_folder];
    system(chmod_cmd);
end
temp_file = [work_folder '/' dir1 '.tmp'];


%% Perform CS reconstruction
% Define auxillary parameters to pass to compiled job
aux_param.mask=mask;
aux_param.maskSize=numel(mask);
aux_param.DN = DN;
aux_param.TVWeight=TVWeight;
aux_param.xfmWeight=xfmWeight;
aux_param.OuterIt=OuterIt;
aux_param.myscale=myscale;
aux_param.scaleFile=scale_file;

aux_param.tempFile=temp_file;
aux_param.totalSlices=dims1(1);

aux_param.dims=dims;
aux_param.dims1=dims1;
aux_param.waveletDims=wavelet_dims;
aux_param.waveletType=wavelet_type;

aux_param.mypdf=mypdf;
aux_param.phmask=phmask;

param = init;
param.Itnlim = Itnlim;




for slice = 1:dims1(1)
    slice_data.(sprintf('s_%i',slice))=data(:,:,slice);
end

time_to_create_master_mat_file=toc

%disp(['Time to create work.mat file = ' time_to_create_master_mat_file  '.'])

% work_folder = [outpath '/' dir1 '/work/']; % Will need to delete after work is completed
%
% if ~exist(work_folder,'dir')
%     dir_cmd2 = ['mkdir ' work_folder '; chmod 777 ' work_folder];
%     system(dir_cmd2);
% else
%     chmod_cmd = ['chmod 777 ' work_folder];
%     system(chmod_cmd);
% end


%% Save common variable file
tic
volume_variable_file = [work_folder dir1 '_workspace.mat'];
if ~exist(volume_variable_file,'file')
    save(volume_variable_file,'-struct','slice_data');
    save(volume_variable_file,'aux_param','-append');
    save(volume_variable_file,'param','-append');
end
time_to_write_master_mat_file=toc
%disp(['Time to write work.mat file = ' time_to_write_master_mat_file  '.'])

%% Create temporary volume for intermediate work
% First [dims(1)] bytes form a header indicating which slices have already
% been reconstructed. Also, we assume that dims1(1) = dims(1) always.

header_size = dims(1);

work_done=zeros([header_size 1]);

if ~exist(temp_file,'file')
    fid=fopen(temp_file,'w');
    new_dims = dims(2)*dims(3)*dims(1);
    vol=zeros([new_dims 1]);
    %vol=[vol(:) vol(:)]; % Don't need intermediate to be complex any more!
    fwrite(fid,work_done,'uint8');
    fseek(fid,header_size,-1);
    fwrite(fid,vol,'double');
    fclose(fid);
    chmod_temp_cmd = ['chmod 777 ' temp_file];
    system(chmod_temp_cmd);
    % else
    %     fid=fopen(temp_file,'r');
    %     work_done=fread(fid,dims(1),'*uint8');
    %     fclose(fid);
end
%
% slices_to_process = find(~work_done)';
%
% num_s2p = length(slices_to_process);
%
%
% msg = ['Reconstructing ' num2str(num_s2p) ' slices.'];
%
% disp(msg)
end

