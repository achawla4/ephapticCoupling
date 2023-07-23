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
lowu = 1e-6; % rho' lower limit
highu = 2e-6; % rho' upper limit
accuracy = 1e-8;
integrand = @(u,v)(((highu - u).*u)./((u.^2 + rho.^2 - 2.*rho.*u.*cos(phi-v) + z.^2).^(3/2)));
Q = dblquad(integrand, lowu, highu, lowv, highv, accuracy, @quadl);

% Second integral
% 
lowv2 = sourcelowv;%10*pi/numchan; % phi' lower limit
highv2 = sourcehighv;%11*pi/numchan; % phi' upper limit
lowu2 = 1e-6; % rho' lower limit
highu2 = 2e-6; % rho' upper limit
accuracy = 1e-8;
integrand2 = @(u,v)(((highu2 - u).*u.*(rho - u.*cos(phi-v)))./((u.^2 + rho.^2 - 2.*rho.*u.*cos(phi-v) + z.^2).^(3/2)));
Q2 = dblquad(integrand2, lowu2, highu2, lowv2, highv2, accuracy, @quadl);

