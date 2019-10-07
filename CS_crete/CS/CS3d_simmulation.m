addpath(genpath('/home/rmd22/Documents/MATLAB/'));

for sampling_fraction = [.25 .33 .4 .5 .6 .7 .8 .9]
    sampling_fraction

%%% this simulation uses data in my nas4 directory
cd /nas4/rmd22/CS_test/3D/
load G
mymask = fullfile(pwd,'mask2.nii');

%%% This is a specific data set of a mouse heart
% img = double(open_nii('nii4D')); %This should probably be the complex data
[RE,IM,NP,NB,NT,HDR] = load_fid_noprocpar(pwd);
img = iftME(reshape(double(complex(RE,IM)),[64 64 64 13]));
% sampling_fraction = 1;

%%% this is code to correct for non-square acquisitions
dims0 = size(img);
dyadic_idx = 2.^[1:10];
pidx = find(max(dims0(2:3))<=dyadic_idx,1);
p = 2^pidx;
img = padarray(img,[0 p-dims0(2) p-dims0(3)]/2,0,'both');
s = max(abs(img(:)));
imgd = img/s;
dims1 = size(img);

imgk = zeros(dims1,'like',img);
data0 = zeros(dims1,'like',imgk);
data1 = data0;
data = zeros(dims1([2 3 1 4]),'like',imgk);
mask3D = true(false);
for k = 1:dims1(4)
   [pdf{k},val{k}] = genPDF(dims1(2:3),2,sampling_fraction,2,0,false);
   [minIntrVec{k},stat{k},actpctg(k)] = genSampling(pdf{k},5,10);
   mask3D(:,:,:,k) = repmat(permute(minIntrVec{k},[3 1 2]),[dims1(1) 1 1]);
   imgk(:,:,:,k) = fftnc(img(:,:,:,k));
   data0(:,:,:,k) = imgk(:,:,:,k).*double(mask3D(:,:,:,k));
end
data1 = data0;
%% Option to fill in k-space using closest orientation(s)
n_o = 0; % number of orientations to use when filling in
G = Gflip1;
DWidx = 2:dims1(4);
for k = DWidx
   ang_rad = acos((G(k,:)*G')/(G(k,:)*G(k,:)'));
   theta=radtodeg(ang_rad);
   theta(theta > 90) = 180-theta(theta > 90);
   [~,idx] = sort(theta,'descend');
   tmp = zeros(dims1(1:3),'like',imgk);
   for i = idx(DWidx(end-n_o:end))
       datai = data0(:,:,:,i);
       tmp(mask3D(:,:,:,i)) = datai(mask3D(:,:,:,i));
   end
   data1(:,:,:,k) = tmp;
end
%%
for k = 1:dims1(4)
   data1(:,:,:,k) = fftshift(ifft(fftshift(data1(:,:,:,k),1),[],1),1); % take ifft in the fully sampled dimension
   data(:,:,:,k) = permute(data1(:,:,:,k),[2,3,1]);
end

%%

im_res = zeros(size(data),'like',imgk);
for k = 1:dims1(4)
    k
dims = size(data(:,:,:,k));
im_zfwdc = zeros(dims(1:3));
for n=1:dims(3)
	im_zfwdc(:,:,n) = ifft2c(data(:,:,n,k)./pdf{k}); % this compensates the intensity for the undersampling
end

% scale data such that the maximum image pixel in zf-w/dc is around 1
% this way, we can use similar lambda for different problems
myscale(k) = max(abs(im_zfwdc(:)));
data(:,:,:,k) = data(:,:,:,k)/myscale(k);
im_zfwdc = im_zfwdc/myscale(k);    

mask = minIntrVec{k};
% data(:,:,:,k) = data(:,:,:,k).*double(repmat(mask,[1 1 dims(3)]));
% pdf = repmat(pdf,[1 1 dims(3)]);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [size(data(:,:,:,k),1),size(data(:,:,:,k),2)]; 		% image Size
DN = N;		 	% data Size
TVWeight = [0.001] ; 	% Weight for TV penalty - only TV is on, but I encourage you to try wavelets as well.
xfmWeight = [0.005];	% Weight for Transform L1 penalty
Itnlim = 8;	%10	% Number of iterations
OuterIt = length(TVWeight);

%generate transform operator

XFM = Wavelet('Daubechies',4,4);	% Wavelet, this transform doesn't work here, why?
% XFM = TIDCT(8,4);			% DCT
% XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
phmask = zpad(hamming(6)*hamming(6)',N(1),N(2)); %mask to grab center frequency
phmask = phmask/max(phmask(:));			 %for low-order phase estimation and correction
res = zeros(N);

param = init;
param.XFM = XFM;
param.TV = TVOP;
param.Itnlim = Itnlim;

tic
for slice = 1:size(data(:,:,:,k),3)

	param.data = data(:,:,slice,k);
	ph = exp(1i*angle((ifft2c(data(:,:,slice,k).*phmask)))); % estimate phase for phase correction
	param.FT = p2DFT(mask, DN, ph, 2); 
    res = XFM*im_zfwdc(:,:,slice); %%% not in original code, set this as x0?

	for n=1:OuterIt
		param.TVWeight =TVWeight(n);     % TV penalty 
		param.xfmWeight = xfmWeight(n);  % L1 wavelet penalty
		res = fnlCg(res, param); %previous image is used as x0 in each step since they are close together
	end
	im_res(:,:,slice,k) = XFM'*res;
%     figure(100);
%     subplot(2,1,1), show(cat(2,abs(im_zfwdc(:,:,slice)),abs(im_res(:,:,slice,k)))), drawnow;
% 	if mod(slice,8)==0
% 		figure(100), subplot(2,1,2), show(cat(2,max(abs(permute(im_zfwdc,[3,1,2])),[],3),max(abs(permute(im_res(:,:,:,k),[3,1,2])),[],3))), drawnow;
% 	end
end

im_res(:,:,:,k) = myscale(k)*im_res(:,:,:,k);

end

% undo ifft2c scaling
im_res = permute(abs(im_res),[3 1 2 4])/sqrt(numel(mask)); %%% add sqrt(numel(mask)) scale factor

out_dir = fullfile(pwd,['LSQR' num2str(n_o) '_it' num2str(Itnlim) '_CS' num2str(100*sampling_fraction) 'pct/']);
mkdir(out_dir)
cd(out_dir);
mat2nii(im_res);
DTI_LSQRc(out_dir,'im_res.nii',Gflip1,bval,1,mymask);


end