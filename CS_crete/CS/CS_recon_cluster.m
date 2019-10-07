function [im_res,scaling] = CS_recon_cluster(reconfile,volume_number,outpath,sampling_fraction)
%% Load data
load(reconfile);
mask = skipint2skiptable(procpar_path); %sampling mask
mask0 = mask;
data = double(load_fidCS(fidpath,max_blocks,ntraces,npoints,bitdepth,volume_number,voldims));
% originalUS = abs(ifftnc(data));
data0 = data;
dims = size(data);

% Get sampling PDF (this is not the sampling mask)
% my_sampling_fraction = 0.5;
if ~exist('sampling_fraction','var')
    sampling_fraction = sum(mask(:))/numel(mask); % estimate if no sampling fraction provided
end
[mypdf,~] = genPDF(dims(2:3),4,sampling_fraction,2,0,false);
mypdf0 = mypdf;

% pad if non-square or non-power of 2
dyadic_idx = 2.^[1:10];
pidx = find(max(dims(2:3))<=dyadic_idx,1);
p = 2^pidx;
data = padarray(data,[0 p-dims(2) p-dims(3)]/2,0,'both');
mask = padarray(mask,[p-dims(2) p-dims(3)]/2,0,'both');
dims1 = size(data);
mypdf = padarray(mypdf,[p-dims(2) p-dims(3)]/2,1,'both'); %pad with 1's since we don't want to divide by zero later

% Calculate group scaling from first b0 image
if ~exist('scaling','var')

    data0 = fftshift(ifft(fftshift(data0,1),[],1),1); % take ifft in the fully sampled dimension
    data0 = permute(data0,[2,3,1]);
    for n = 1:dims(1)
%         data0(:,:,n) = ifftnc(data0(:,:,n)./mypdf); % this compensates the intensity for the undersampling
        data0(:,:,n) = fftshift(ifftn(fftshift(data0(:,:,n)./mypdf0))); % this compensates the intensity for the undersampling
    end
    
    thresh = .999; % only clip out very noisy voxels (e.g. zippers, artifacts), and not the eyes
    q = quantile(abs(data0(:)),thresh);
    scaling = (2^16-1)/q; % we plan on writing out uint16 not int16, though it won't show up well in the database
    m = matfile(reconfile,'Writable',true);
    m.scaling = scaling;
    
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
im_zfwdc = im_zfwdc/myscale;    

% Define recon parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = dims1(2:3); 		% image Size
DN = N;		 	% data Size
TVWeight = 0.0009; 	% Weight for TV penalty - only TV is on, but I encourage you to try wavelets as well.
% xfmWeight = 0.005;	% Weight for Transform L1 penalty
xfmWeight = 0.006;	% Weight for Transform L1 penalty
Itnlim = 24;	%10	% Number of iterations
OuterIt = length(TVWeight);

%generate transform operator

% XFM = Wavelet('Daubechies',4,4);	% Wavelet, this transform doesn't work here, why?
XFM = Wavelet('Daubechies',12,12);	% Wavelet, this transform doesn't work here, why?

% XFM = TIDCT(8,4);			% DCT
% XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
% phmask = zpad(hamming(6)*hamming(6)',N(1),N(2)); %mask to grab center frequency
phmask = zpad(hamming(48)*hamming(48)',N(1),N(2)); %mask to grab center frequency
phmask = phmask/max(phmask(:));			 %for low-order phase estimation and correction

param = init;
param.XFM = XFM;
param.TV = TVOP;
param.Itnlim = Itnlim;

%% Make and save header data

myvolstr = num2str(volume_number-1);
if volume_number-1 < 10
    myvolstr = ['00' myvolstr];
elseif volume_number-1 < 100
    myvolstr = ['0' myvolstr];
end
% save_nii(make_nii(abs(im_res),[1 1 1],[0 0 0],16),fullfile(outpath,['img_' myvolstr '.nii']));
dir1 = [runno '_m' myvolstr];
voldir = fullfile(outpath,dir1,[dir1 'images']);
mkdir(voldir);

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
struct1.A_max_bval = procpar.max_bval;
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
% struct1.work_dir_path = ['/naxosspace/' runno '.work'];
struct1.work_dir_path = ['/delosspace/' runno '.work'];
struct1.CS_sampling_fraction = sampling_fraction;
struct1.CS_acceleration = 1/sampling_fraction;
struct1.group_image_scale = scaling;

write_headfile(fullfile(voldir,[dir1 '.headfile']),struct1,'',0);

%% Perform CS reconstruction

for slice = 1:dims1(1)
	param.data = data(:,:,slice);
	ph = exp(1i*angle((ifft2c(data(:,:,slice).*phmask)))); % estimate phase for phase correction
	param.FT = p2DFT(mask, DN, ph, 2); 
    res = XFM*im_zfwdc(:,:,slice); 

	for n=1:OuterIt
		param.TVWeight =TVWeight(n);     % TV penalty a
		param.xfmWeight = xfmWeight(n);  % L1 wavelet penalty
		res = fnlCg(res, param); 
	end
	im_res(:,:,slice) = XFM'*res;
    
    %%% visualize results
    %     figure(100);
    %     subplot(2,1,1), show(cat(2,abs(im_zfwdc(:,:,slice)),abs(im_res(:,:,slice)))), drawnow;
    % 	if mod(slice,8)==0
    % 		figure(100), subplot(2,1,2), show(cat(2,max(abs(permute(im_zfwdc,[3,1,2])),[],3),max(abs(permute(im_res,[3,1,2])),[],3))), drawnow;
    % 	end
end

% Undo CS-speific scaling and permutation to final form
im_res = permute(im_res,[3 1 2 4]);
im_res = im_res*myscale/sqrt(numel(mask)); % --->> check this, sqrt of 2D plane elements required for proper scaling

% crop out extra k-space if non-square or non-power of 2
if sum(dims == dims1) ~= 3
    im_res = fftshift(fftn(fftshift(im_res)));
    im_res = im_res(:,(dims1(2)-dims(2))/2+1:end-(dims1(2)-dims(2))/2, ...
                    (dims1(3)-dims(3))/2+1:end-(dims1(3)-dims(3))/2);
    im_res = fftshift(ifftn(fftshift(im_res)));
end

% Apply group scaling 
im_res = uint16(abs(im_res)*scaling);

%% Save data

% write_archive_tag(dir1,'/naxosspace',dims(3),struct1.U_code, ...
%                   ['.' struct1.U_stored_file_format],struct1.U_civmid,true,voldir);                              
% write_civm_image(fullfile(voldir,[dir1 struct1.scanner_tesla_image_code 'imx']), ...
%                           im_res,struct1.U_stored_file_format,0,1);

write_archive_tag(dir1,'/delosspace',dims(3),struct1.U_code, ...
                  ['.' struct1.U_stored_file_format],struct1.U_civmid,true,voldir);                              
write_civm_image(fullfile(voldir,[dir1 struct1.scanner_tesla_image_code 'imx']), ...
                          im_res,struct1.U_stored_file_format,0,1);


% % Copy data to naxosspace for now
% pw = '4.signa!';
% copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(outpath,dir1) ...
%             ' omega@naxos.duhs.duke.edu:/Volumes/nexosspace/' dir1];
% system(copy_cmd);
% copy_archivetag_cmd = ['sshpass -p ' pw ' scp -p ' fullfile(voldir,['READY_' dir1]) ...
%             ' omega@naxos.duhs.duke.edu:/Volumes/delosspace/Archive_Tags/READY_' dir1];
% system(copy_archivetag_cmd);

% Copy data to naxosspace for now
pw = '4.signa!';
copy_cmd = ['sshpass -p ' pw ' scp -rp ' fullfile(outpath,dir1) ...
            ' omega@delos.duhs.duke.edu:/Volumes/delosspace/' dir1];
system(copy_cmd);
copy_archivetag_cmd = ['sshpass -p ' pw ' scp -p ' fullfile(voldir,['READY_' dir1]) ...
            ' omega@delos.duhs.duke.edu:/Volumes/delosspace/Archive_Tags/READY_' dir1];
system(copy_archivetag_cmd);

end