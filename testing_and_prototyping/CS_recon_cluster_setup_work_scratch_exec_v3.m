function [ ] = CS_recon_cluster_setup_work_scratch_exec_v3(variables_file,volume_number)
%CS_RECON_CLUSTER_SETUP_WORK_EXEC An executable MATLAB script for setting
%up each volume of CS reconstruction in order to avoid saturating the
%master node (in the context of DTI, with many many volumes to recon.

%% Update of original version (implied _v1)

if ~isdeployed
    %variables_file = '/glusterspace/N54806_69/N54806_69_m00/work/N54806_69_m00_setup_variables.mat';
    %variables_file = '/glusterspace/N54799/N54799_m0/work/N54799_m0_setup_variables.mat';
    variables_file = '/glusterspace/S67461_69/S67461_69_m0/work/S67461_69_m0_setup_variables.mat';
    volume_number = '1';
    
end

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

n_sampled_points = sum(mask(:));

volume_number=str2double(volume_number);


%% Immediately check to see if volume has already been reconstructed
if  ~exist( 'do_russway','var')
    max_number = nvols;
    if (nechoes > 1)
        max_number = nechoes;
    end
    
    myvolstr =sprintf(['%0' num2str(numel(num2str(max_number-1))) 'i' ],volume_number-1);
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
only_non_zeros = 1;


fidpath=[outpath '/' runno '_m' myvolstr '.fid'] %%NEW
temp_volume_number=1; %%NEW
tic
%data=single(zeros([floor(npoints/2),n_sampled_points]));
double_down = 0; % Need to test this feature, if we want to keep double processing for all fft operations, and beyond.
data = load_fidCS(fidpath,max_blocks,ntraces,npoints,bitdepth,temp_volume_number,voldims,only_non_zeros,double_down);
%data = double(data); %This should be replaced by setting double_down to 1.
disp('fid data has been successfully loaded');
toc

%data=double(data);


if (nechoes > 1)
    data = reshape(data,[npoints/2 nechoes ntraces/nechoes]);
    data = permute(data,[1 3 2]);
end


tic
data = fftshift(ifft(fftshift(data,1),[],1),1); % take ifft in the fully sampled dimension
toc

d_data=size(data);
dims = [d_data(1) size(mask)];%size(data);

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


[mypdf,~] = genPDF_wn_v2(dims(2:3),pa,sampling_fraction,pb,false);

%mypdf0 = mypdf;

%if (~exist('scale_file','file') && (volume_number==1))
%if ~exist('scaling','var')
mask0 = mask;
mypdf0 = mypdf;
%end
%end

%
% % pad if non-square or non-power of 2

dyadic_idx = 2.^(1:14); %dyadic_idx = 2.^[1:12]; %%%% 12->14
pidx = find(max(dims(2:3))<=dyadic_idx,1);
p = 2^pidx;

if (p>max(dims(2:3)))
    mask = padarray(mask,[p-voldims(2) p-voldims(3)]/2,0,'both');
    mypdf = padarray(mypdf,[p-dims(2) p-dims(3)]/2,1,'both'); %pad with 1's since we don't want to divide by zero later
end
dims1 = [d_data(1) size(mask)];%size(data);

max_number = nvols;
if (nechoes > 1)
    max_number = nechoes;
end

thresh = .999999999; % only clip out very noisy voxels (e.g. zippers, artifacts), and not the eyes
for vol_number = 1:nechoes
    if (nechoes > 1)
        volume_number = vol_number;
        current_data=squeeze(data(:,:,volume_number));
    else
        current_data=data;
        clear data;
    end
    
    
    % Test again for the need to do work, as is done in the calling
    % function
    
    if ((nechoes > 1) && (volume_number > 1))
        if  ~exist( 'do_russway','var')
            max_number = nvols;
            if (nechoes > 1)
                max_number = nechoes;
            end
            myvolstr =sprintf(['%0' num2str(numel(num2str(max_number-1))) 'i' ],volume_number-1);
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
    end
    
    finished_slices = dir( [voldir '/*.raw' ]);
    finished_slices_count = length(finished_slices(not([finished_slices.isdir])));
    should_I_do_work=0;
    if (finished_slices_count ~= voldims(3))
        should_I_do_work=1;
    end
    
    I_really_should_do_work = 0;
    if should_I_do_work
         work_folder = [outpath '/' dir1 '/work/'];
        volume_variable_file = [work_folder dir1 '_workspace.mat'];
        try
            dummy = load(volume_variable_file,'aux_param.maskSize'); % Need to try to load an arbitrary variable from the work file
            clear dummy;
        catch
            I_really_should_do_work = 1;
            if exist(volume_variable_file,'file');
                cmd=['rm ' volume_variable_file]
                system(cmd);
                'should be removing something'
            end
        end
    end
    
    if ((I_really_should_do_work) || (~exist('scale_file','file')&& (volume_number==1)))
        %% Calculate group scaling from first b0 image
        if (~exist('scale_file','file') && (volume_number==1))
            
            tic
            current_slice=zeros([size(mask0)],'like',current_data);
            for n = 1:dims(1)
                current_slice(mask0(:))=current_data(n,:);
                temp_data = abs(ifftn(current_slice./mypdf0)); % 8 May 2017, BJA: Don't need to waste computations on fftshift for scaling calculation
                qq(n)= max(temp_data(:));%quantile(temp_data(:),thresh);
            end
            toc
            
            q = quantile(qq,thresh);
            scaling = (2^16-1)/q; % we plan on writing out uint16 not int16, though it won't show up well in the database
            scaling = double(scaling);
            
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
        end
        
        
        
        %% Prep data for reconstruction
        
        %%%%%% One could potentially zero-pad non-square data here in order to use
        %%%%%% the Wavelet transform version of CS.
        
        % Calculate scaling
        if (exist('scaling','var') && (volume_number == 1) && ~(sum((dims1 - dims))))
             myscale = sqrt(dims1(2)*dims1(3))*(2^16-1)/scaling; % We've already done the heavy lifting for this calculation, if array size doesn't change.
        else
            tic
            current_slice=zeros(size(mask),'like',current_data);
            qq=zeros([1 dims(1)]);
            for n = 1:dims(1)
                current_slice(mask(:))=current_data(n,:);
                temp_data = abs(ifftn(current_slice./mypdf)); % 8 May 2017, BJA: Don't need to waste computations on fftshift for scaling calculation
                qq(n)=max(temp_data(:));%quantile(temp_data(:),thresh);
            end
            
            myscale = sqrt(dims1(2)*dims1(3))*quantile(qq,thresh);
  
            fid = fopen(scale_file,'r');
            scaling = fread(fid,inf,'float');
            fclose(fid);
            optimized_for_memory_time = toc
        end
 
        
        % scale data such that the maximum image pixel in zf-w/dc is around 1
        % this way, we can use similar lambda for different problems
        current_data = current_data/myscale;
        
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
        
        phmask = zpad(hamming(32)*hamming(32)',N(1),N(2)); %mask to grab center frequency
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
        
        
        %%
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
        aux_param.originalMask=mask0;
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
        aux_param.originalMypdf=mypdf0;
        
        aux_param.phmask=phmask;
        
        param = init;
        param.Itnlim = Itnlim;  % Should this be a function of necho?
        
        
        %% Save common variable file
        volume_variable_file = [work_folder dir1 '_workspace.mat'];
        
        if ~exist(volume_variable_file,'file')
            
            real_data = real(current_data);
            imag_data = imag(current_data);
            savefast2(volume_variable_file,'real_data','imag_data');
            save(volume_variable_file,'aux_param','-append');
            save(volume_variable_file,'param','-append');
            
            time_to_write_master_mat_file=toc
            
            
            %% Create temporary volume for intermediate work
            % First [dims(1)] bytes form a header indicating which slices have already
            % been reconstructed. Also, we assume that dims1(1) = dims(1) always.
            
            %header_size = dims(1)*64;
            header_size = dims(1)*16; % 15 May 2017, BJA: moving to unit16 for header info
            
            work_done=zeros([header_size 1]);
            
            if ~exist(temp_file,'file')
                fid=fopen(temp_file,'w');
                new_dims = dims(2)*dims(3)*dims(1);
                vol=zeros([new_dims 1]);
                %vol=[vol(:) vol(:)]; % Don't need intermediate to be complex any more!
                %fwrite(fid,work_done,'uint8');
                fwrite(fid,work_done,'uint16');
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
            
        end
    end
end
