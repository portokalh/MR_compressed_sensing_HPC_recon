function [outXYZ] = rotateVectors(theta, phi, xyz)
% phi_rad = phi*pi/180;
% theta_rad = theta*pi/180;
Mphi = [1 0 0
    0 cos(phi) -sin(phi)
    0 sin(phi) cos(phi)];
Mthet = [cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];
outXYZ = Mthet*Mphi*xyz;
% [theta_in phi_in, r_in] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
% [xout yout zout] = sph2cart(theta_in+theta_rad, phi_in+phi_rad, r_in);
end