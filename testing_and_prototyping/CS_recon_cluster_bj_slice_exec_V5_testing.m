function CS_recon_cluster_bj_slice_exec_V5(matlab_workspace,slice_indices)
if ~isdeployed
    matlab_workspace ='/glusterspace/S67668/S67668_m13/work/S67668_m13_workspace.mat';
    
    slice_indices ='256_to_260';
end
tic
slice_numbers=[];
slice_number_strings = strsplit(slice_indices,'_');
for ss = 1:length(slice_number_strings)
    temp_string=slice_number_strings{ss};
    if strcmp(temp_string,'to')
        begin_slice = str2double(slice_number_strings{ss-1})+1;
        end_slice = str2double(slice_number_strings{ss+1})-1;
        temp_vec = begin_slice:1:end_slice;
    else
        temp_vec = str2double(temp_string);
    end
    slice_numbers = [slice_numbers temp_vec];
end
slice_numbers=unique(slice_numbers)
toc

if ~isdeployed
    addpath(genpath('/home/rmd22/Documents/MATLAB/'));
end

%% Make sure the goddam workspace file exist

for tt = 1:100 % Will check every 15 seconds, up to 1500 seconds (25 min)
    if  (~exist(matlab_workspace,'file'))
        msg = [matlab_workspace ' not found; waiting 15 seconds...'];
        disp(msg)
        pause(15)
    else
        a = who('-file',matlab_workspace,'aux_param');
        if ~size(a)
            msg = [matlab_workspace ' found, but aux_param has not been written yet; waiting 15 seconds...'];
            disp(msg)
            pause(15)
        end
    end
end


%% Load common workspace params
tic
load(matlab_workspace,'aux_param');
load(matlab_workspace,'param');
time_to_load_common_workspace=toc

%% Wait for scale_file to exist and then read
scale_file=aux_param.scaleFile;

while (~exist(scale_file,'file'))
    pause(15)
end

fid_sc = fopen(scale_file,'r');
scaling = fread(fid_sc,inf,'*float');
fclose(fid_sc);




%% Setup common variables
mask_size=aux_param.maskSize;
mask=aux_param.mask;
DN=aux_param.DN;
TVWeight=aux_param.TVWeight;
xfmWeight=aux_param.xfmWeight;
OuterIt=aux_param.OuterIt;
myscale=aux_param.myscale;
temp_file=aux_param.tempFile;
dims=aux_param.dims;
dims1=aux_param.dims1;
mypdf=aux_param.mypdf;
phmask=aux_param.phmask;


wavelet_dims=aux_param.waveletDims;

if isfield(aux_param,'waveletType')
    wavelet_type=aux_param.waveletType;
else
    wavelet_type = 'Daubechies';
end

XFM = Wavelet(wavelet_type,wavelet_dims(1),wavelet_dims(2));

param.XFM = XFM;
param.TV=TVOP;

im_result=zeros(dims(2),dims(3),length(slice_numbers));


%% Reconstruct slice(s)
for index=1:length(slice_numbers)
    slice_index=slice_numbers(index);
    
    
    fid=fopen(temp_file,'r+');
    work_done = fread(fid,dims(1),'*uint8');
    fclose(fid);
    if 1 %if ~work_done(slice_index)
        
        
        % Load slice specific data
        tic
        
        slice_data_var_name=(sprintf('s_%i',slice_index));
        package=load(matlab_workspace,slice_data_var_name);
        package=package.(slice_data_var_name);
        %package=eval(slice_data_var_name);
        
        time_to_load_workspace = toc
        
        
        tic
        param.data=package(:,:,1);%aux_param.data(:,:,slice_index);
        
        im_zfwdc = ifft2c(param.data./mypdf)/myscale; % this compensates the intensity for the undersampling
        
        %res=package(:,:,2);
        
        res=XFM*im_zfwdc;
        
        
         plotting_today_BJ = 1;
        if plotting_today_BJ
            scale_file=aux_param.scaleFile;
            fid_sc = fopen(scale_file,'r');
            scaling = fread(fid_sc,inf,'*float')
            fclose(fid_sc);    
            i_im_res = XFM'*res;
            i_im_res = i_im_res*myscale/sqrt(mask_size); % --->> check this, sqrt of 2D plane elements required for proper scaling
        
        
        %% Crop out extra k-space if non-square or non-power of 2
        %if sum(original_dims == recon_dims) ~= 3
        %    im_res = fftshift(fftn(fftshift(im_res)));
        %    im_res = im_res((recon_dims(2)-original_dims(2))/2+1:end-(recon_dims(2)-original_dims(2))/2, ...
        %       (recon_dims(3)-original_dims(3))/2+1:end-(recon_dims(3)-original_dims(3))/2);
        %    im_res = fftshift(ifftn(fftshift(im_res)));
        
        figure(2000+slice_index)
        im_to_plot = double(abs(i_im_res')*scaling);
        imagesc(im_to_plot)
        colormap gray
        axis xy
        
        end

        %ph=package(:,:,3);
        ph = exp(1i*angle((ifft2c(param.data.*phmask))));
        
        param.FT = p2DFT(mask, DN, ph, 2);
        time_to_set_up = toc
        
        tic
        for n=1:OuterIt
            param.TVWeight =TVWeight(n);     % TV penalty
            param.xfmWeight = xfmWeight(n);  % L1 wavelet penalty
            res = fnlCg(res, param);
        end
        time_to_recon = toc
        
        tic
        im_res = XFM'*res;
        im_res = im_res*myscale/sqrt(mask_size); % --->> check this, sqrt of 2D plane elements required for proper scaling
        
        %% Crop out extra k-space if non-square or non-power of 2
        if sum(dims == dims1) ~= 3
            im_res = fftshift(fftn(fftshift(im_res)));
            im_res = im_res((dims1(2)-dims(2))/2+1:end-(dims1(2)-dims(2))/2, ...
                (dims1(3)-dims(3))/2+1:end-(dims1(3)-dims(3))/2);
            im_res = fftshift(ifftn(fftshift(im_res)));
        end
        im_result(:,:,index)=double(abs(im_res)*scaling);
        
        figure(slice_index)
        im_to_plot = double(abs(im_res')*scaling);
        imagesc(im_to_plot)
        colormap gray
        axis xy
    else
        msg = ['Slice ' num2str(slice_index) ' appears to already have been processed; skipping.'];
        disp(msg)
        
    end

end
%% Apply group scaling
%im_result = double(abs(im_result)*scaling);
time_to_post_process=toc

%% Prepare to write slices and header info
tic
header_size = dims(1);
if 0
for index=1:length(slice_numbers)
    slice_index=slice_numbers(index);
    im_to_write = im_result(:,:,index);
    image_to_write(:,index) = typecast(im_to_write(:),'uint8');
    data_offset(index)= header_size + (8*dims(2)*dims(3)*(slice_index-1));
end
end
time_to_prepare_write_data = toc

%% Write data
tic

fid=fopen(temp_file,'r+');

work_done=fread(fid,dims(1),'*uint8');

for index=1:length(slice_numbers)
    slice_index=slice_numbers(index);
    if ~work_done(slice_index)
        fseek(fid,data_offset(index),-1);
        fwrite(fid,image_to_write(:,index),'uint8'); %'n'
        msg_1 = ['Slice ' num2str(slice_index) ' successfully reconstructed and written to ' temp_file '.'];
        disp(msg_1);
%    end
%end

%for index=1:length(slice_numbers)
%    slice_index=slice_numbers(index);
%    if ~work_done(slice_index)
        fseek(fid,(slice_index-1),-1);
        fwrite(fid,1,'uint8');
        msg_2 = ['Successful reconstruction flag written to header of ' temp_file ' for slice ' num2str(slice_index) '.'];
        disp(msg_2);
    end
end

fclose(fid);

time_to_write_data=toc
end