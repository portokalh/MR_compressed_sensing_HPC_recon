clear
place=2;%1 Mac desktop omega,2 Mac stick 3 PC work 4 PC home
iterations=10
resxy=360
acceleration =8
%acceleration0 =8
%acceleration =acceleration0*resxy/256

dims_xy=[resxy resxy]
pa=2.3
paf=pa*10
pbtest=resxy*0.05
pb=5.6
%pb0= 5.6
%pb=pb0*resxy/256
%pb=~fraction * 0.05
pbf=pb*10
if acceleration >= 1
sampling_fraction = 1/acceleration;
else
sampling_fraction = acceleration;
acceleration = 1/sampling_fraction;
end
if numel(dims_xy) == 1
dims_xy = [dims_xy dims_xy];
end


filenameout = ['CS' num2str(dims_xy(1)) '_' num2str(acceleration) 'x_pa' num2str(paf) '_' 'pb' num2str(pbf) '.txt'];

sampling_name_png=['CStabletest','.png'];


[mypdf,val] = genPDF_wn_v2(dims_xy,pa,sampling_fraction ,pb,1);

[mask, stat, actpctg]= genSampling(mypdf,iterations,1);
mask=reshape(mask,resxy*resxy,1)
mask=transpose(mask)
mask(end)=1 
if (place==1)
    mytablefname=['/Volumes/CivmUsers/omega/Desktop/CS' num2str(dims_xy(1)) '_' num2str(acceleration) 'x_pa' num2str(paf) '_' 'pa' num2str(pbf) ];%MAC
end
if (place==2)
    mytablefname=['/Volumes/256GBYTE/Lustig_Wang/CStables/CS' num2str(dims_xy(1)) '_' num2str(acceleration) 'x_pa' num2str(paf) '_' 'pa' num2str(pbf) ];%MAC
   
end

if (place==3)
    mytablefname=['I:\Lustig_Wang\CStables\CS' num2str(dims_xy(1)) '_' num2str(acceleration) 'x_pa' num2str(paf) '_' 'pa' num2str(pbf) ];%PC Wo
end
if (place==4)
     mytablefname=['J:\Lustig_Wang\CStables\CS' num2str(dims_xy(1)) '_' num2str(acceleration) 'x_pa' num2str(paf) '_' 'pb' num2str(pbf)  ];%PC Home
end

dlmwrite([mytablefname '.txt'],mask,'')

im8=reshape(mask,dims_xy);


%figure(8)
figure1=figure('Name','GPCmask8')
imagesc(im8)
axis square
filename=[mytablefname, '_ab8.png'];
print(filename,'-dpng', '-r300');
