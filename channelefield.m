function [Efieldx, Efieldy, Efieldz] = channelefield(rprimeplus, thetaprimeplus, phiprimeplus, rprimeminus, thetaprimeminus, phiprimeminus, r, theta, phi, choice)

deltaqdoubleprime = 12e-6; % known charge transferred through channel in cgs units (statCoulombs).
rnor = 0.5e-6*100; % radius of node of ranvier in cgs units
d = 1e-9*100; %thickness of cell membrane at node in cgs units
g = rnor/1000000; % arbitrary dumbell radius in cgs units
s = d+(3*g/4); % separation of dipolar charges in cgs units

rminusprimemag = rprimeminus; 
rplusprimemag = rprimeplus;
fieldposmag = r; 

phiprimeplus = phiprimeplus;
phiprimeminus = phiprimeminus; 
thetaminusprime = thetaprimeminus; 

% Find the required angle
rminusprimedotr = r*rminusprimemag*(sin(phi)*sin(phiprimeminus)*cos(theta - thetaminusprime) + cos(phi)*cos(phiprimeminus));

reqangle = acos(rminusprimedotr/(rminusprimemag*r)); % assuming r is nonnegative.


[rminusprimex, rminusprimey, rminusprimez] = spher2cart(rprimeminus, thetaprimeminus, phiprimeminus)
[rplusprimex, rplusprimey, rplusprimez] = spher2cart(rprimeplus, thetaprimeplus, phiprimeplus);
[x,y,z] = spher2cart(r, theta, phi)

G = [x,y,z]-[rminusprimex, rminusprimey, rminusprimez];
[rminusr, rminustheta, rminusphi] = cart2spher(G(1), G(2), G(3));
rminusmag = sqrt(rminusprimemag^2 + fieldposmag^2 - 2*fieldposmag*rminusprimemag^2*cos(reqangle));

% Thus we can find the unit vector n

N = (1/rminusmag).*G;
if choice == 1
% Alternate dipole moment method
rprimemag = 0.5*(rplusprimemag + rminusprimemag);
thetaprime = 0.5*(thetaprimeplus + thetaprimeminus);
phiprime = 0.5*(phiprimeplus + phiprimeminus);
f = rprimemag*cos(thetaprime);
[AA, BB, CC] = spher2cart(rprimemag, thetaprime, phiprime)
[DD, EE, FF] =  spher2cart(f, 0, 0)
II = [AA, BB, CC]-[DD, EE, FF];
requnit = sqrt(II(1)^2 + II(2)^2 + II(3)^2);
% Thus
P = (deltaqdoubleprime/requnit).*[II(1), II(2), II(3)];


else 
% Find the dipole moment
P = deltaqdoubleprime.*([rplusprimex, rplusprimey, rplusprimez] - [rminusprimex, rminusprimey, rminusprimez]);
end
% Find the dot product of unit vector n with the dipole moment
ndotp = P(1)*N(1) + P(2)*N(2) + P(3)*N(3);

% Find the denominator

DEN = [x, y, z] - 0.5.*([rminusprimex, rminusprimey, rminusprimez] + [rplusprimex, rplusprimey, rplusprimez])
%magnitude of denominator
magden = sqrt(DEN(1)^2 + DEN(2)^2 + DEN(3)^2);
denfinal = magden^3; 

% find the field
Efield = (3*ndotp.*N - P)./denfinal;

Efieldx = Efield(1);
Efieldy = Efield(2);
Efieldz = Efield(3);
