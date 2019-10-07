xyz = [
-0.76600940	-0.57960717	0.27800203
-0.76600940	0.57960717	0.27800203
-0.57960717	-0.27800203	0.76600940
-0.57960717	0.27800203	0.76600940
-0.27800203	-0.76600940	0.57960717
-0.27800203	0.76600940	0.57960717
0.27800203	-0.76600940	0.57960717
0.27800203	0.76600940	0.57960717
0.57960717	-0.27800203	0.76600940
0.57960717	0.27800203	0.76600940
0.76600940	-0.57960717	0.27800203
0.76600940	0.57960717	0.27800203
-0.85065081	-0.52573111	0.00000000
-0.85065081	0.52573111	0.00000000
-0.52573111	0.00000000	0.85065081
0.00000000	-0.85065081	0.52573111
0.00000000	0.85065081	0.52573111
0.52573111	0.00000000	0.85065081
-0.93417236	0.00000000	0.35682209
0.93417236	0.00000000	0.35682209
0.00000000	-0.35682209	0.93417236
0.00000000	0.35682209	0.93417236
-0.35682209	-0.93417236	0.00000000
-0.35682209	0.93417236	0.00000000
0.57735027	-0.57735027	0.57735027
0.57735027	0.57735027	0.57735027
-0.57735027	-0.57735027	0.57735027
-0.57735027	0.57735027	0.57735027
-1.00000000	0.00000000	0.00000000
0.00000000	0.00000000	1.00000000
0.00000000	1.00000000	0.00000000
-0.75142188	0.00000000	0.65982207
-0.65982207	-0.75142188	0.00000000
-0.65982207	0.75142188	0.00000000
0.00000000	-0.65982207	0.75142188
0.00000000	0.65982207	0.75142188
0.75142188	0.00000000	0.65982207



]';

npts = length(xyz);

% rotate

[xSphere,ySphere,zSphere] = sphere(100);

pt1 = [0 0 1];

alphaVal = .5;


nframes = 50000;
primeplus = 360*(3-sqrt(5))/2;

% nframes=floor(nframes/2)*2+1; %Must have odd number of frames
cview=floor(nframes/2)+1;     %Center frame

is = 0:nframes-1; %In Gary's code, i=acview_start

z_coords = abs(1-(is/cview)); %In Gary's code f=fThing
angs = primeplus.*is.*(pi/180); %azimuthal Angle in radians
ds = sqrt(1-(z_coords.^2));
x_coords = ds.*cos(angs);
y_coords = ds.*sin(angs);

%Handle negatives
z_coords = z_coords.*((2*(is<=cview))-1);

%normalize
ivec_lengths = 1./sqrt((x_coords.^2) + (y_coords.^2) + (z_coords.^2));
xs = x_coords.*ivec_lengths;
ys = y_coords.*ivec_lengths;
zs = z_coords.*ivec_lengths;

figure(2);
hold off;

thetas = acos(zs);
phis = atan2(y_coords,x_coords);
cond_number = zeros(size(thetas));
for iFrame=1:nframes
    [outXYZ] = rotateVectors(thetas(iFrame), phis(iFrame), xyz);
%     figure(1);
%     surfl(xSphere, ySphere, zSphere);
%     alpha(alphaVal);
%     shading interp;
%     hold on;
%     plot3(outXYZ(1,:),outXYZ(2,:),outXYZ(3,:),'.','MarkerSize',15);
    
%     hold off;
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     drawnow;
    %     axis([-1 1 -1 1 -1 1]);
    
    cond_number(iFrame) = calc_cond_number(outXYZ(1,:),outXYZ(2,:),outXYZ(3,:));
%     figure(2);
%     plot(iFrame,cond_number(iFrame),'.');
%     axis([0 nframes+1 0 10]);
%     xlabel('Rotation Iteration');
%     ylabel('Condition number');
%     hold on;
    test =1;
end
    figure(2);
    plot(1:nframes,cond_number,'.');
    %axis([0 nframes+1 1.5810 1.5813]);
    axis([0 nframes+1 1.3 2.6]);
    xlabel('Rotation Iteration');
    ylabel('Condition number');
    sampling_name = 'TestName '
    title([sampling_name, ' Npts=' num2str(npts) ', Mean= ' num2str(mean(cond_number)) ', std= ' num2str(std(cond_number))]);
    hold on;
    
    disp(['Mean= ' num2str(mean(cond_number)) ', std= ' num2str(std(cond_number))]);
    
    
