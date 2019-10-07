xyz = [
-0.22625179716522253	-0.4383142333597305	0.8698797371555327	
0.09100990098617207	0.9093013608570719	0.40606432134079995	
-0.8616261890680835	0.4178807577373823	0.28805552003883267	
0.6932118651890893	0.3788684744678597	0.613119881438745	
0.7125571131797109	-0.4009377414426551	0.5757701693765137	
-0.8392129337719915	-0.4855957070756299	0.24478247702710665	
-0.2455970870344242	0.34149197961414235	0.9072294630904685	
0.13684499205081127	-0.9383081496715114	0.31757088092367336		
]';


npts = length(xyz);
% rotate
basename = 'APS';
sampling_name=['APS',num2str(npts)];
sampling_name_jpg=['APS',num2str(npts),'.jpg'];

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

%caculating  zero rotation condition number
for iNpts=6:npts
    % Calculate condition numbers
    cond_number0(1) = calc_cond_number(x(1:1), y(1:1), z(1:1));
end

    figure(2);
    plot(1:nframes,cond_number,'.');
    axis([0 nframes+1 1.3 2.6]);
    xlabel('Rotation Iteration');
    ylabel('Condition number');
    title([sampling_name,' ConditionNumber0=' num2str(cond_number(1)), ' Npts=' num2str(npts) ', Mean=' num2str(mean(cond_number)) ', std=' num2str(std(cond_number))]);
    hold on;
    saveas(figure(2),sampling_name_jpg);

    disp(['Sampling Name = ' num2str(sampling_name) , npts]);
    disp(['Mean = ' num2str(mean(cond_number)) ', std = ' num2str(std(cond_number))]);
    
%caculating  zero rotation condition number
for iNpts=6:npts
    % Calculate condition numbers
    cond_number0(1) = calc_cond_number(x(1:1), y(1:1), z(1:1));
end
disp(['Condition Number0= ' num2str(cond_number(1))]);

