function [Q, Q2] = zandrhointegralsmod2(sourcelowv, sourcehighv,rho, phi, z)
% This method takes in the field coordinates, rho, phi and z and outputs
% the z-integral as Q, the first output and the rho integral as Q2, the
% second output. Note that these are not the field components, but the
% integrals involved in computing the components. The primed coordinates
% are the source coodrinates and the unprimed represent the field
% positions.

format long
lowv = sourcelowv;%10*pi/numchan; % phi' lower limit
highv = sourcehighv;%11*pi/numchan; % phi' upper limit
lowu = 1e-6; % rho' lower limit. this should be the radius of the node 
highu = 2e-6; % rho' upper limit. this should be the value of p where the current decays to 0.
accuracy = 1e-8;
integrand = @(u,v)(((2.*u - highu))./((u.^2 + rho.^2 - 2.*rho.*u.*cos(phi-v) + z.^2).^(3/2)));
Q = dblquad(integrand, lowu, highu, lowv, highv, accuracy, @quadl);

% Correction
integrandc = @(v)((1)./((lowu.^2 + rho.^2 - 2.*rho.*lowu.*cos(phi-v) + z.^2).^(3/2)));
Qc = integral(integrandc, lowv, highv);
Qc = -((highu - lowu)*lowu)*Qc;


Q = Q+Qc;
% Second integral
% 
lowv2 = sourcelowv;%10*pi/numchan; % phi' lower limit
highv2 = sourcehighv;%11*pi/numchan; % phi' upper limit
lowu2 = 1e-6; % rho' lower limit
highu2 = 2e-6; % rho' upper limit
accuracy = 1e-8;
integrand2 = @(u,v)(((2.*u - highu2).*(rho - u.*cos(phi-v)))./((u.^2 + rho.^2 - 2.*rho.*u.*cos(phi-v) + z.^2).^(3/2)));
Q2 = dblquad(integrand2, lowu2, highu2, lowv2, highv2, accuracy, @quadl);

% Correction
integrand2c = @(v)((1).*(rho - lowu.*cos(phi-v))./((lowu.^2 + rho.^2 - 2.*rho.*lowu.*cos(phi-v) + z.^2).^(3/2)));
Q2c = integral(integrand2c, lowv2, highv2);
Q2c = -((highu2 - lowu2)*lowu2)*Q2c;

Q2 = Q2 + Q2c;

