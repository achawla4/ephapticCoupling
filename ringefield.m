function [Efieldx, Efieldy, Efieldz] = ringefield(rprime, thetaprime, r, theta, phi, choice)
deltaqdoubleprime = 12e-6; % known charge transferred through channel in cgs units
rnor = 0.5e-6*1e2; % radius of node of ranvier in cgs units
d = 1e-9*1e2; %thickness of cell membrane at node in cgs units
g = rnor/1000000; % arbitrary dumbell radius in cgs units
s = d+(3*g/4); % separation of dipolar charges in cgs units



rprimemag = abs(rprime);
f = rprimemag*cos(thetaprime);


rprimeplusmag = sqrt(f^2 + (rnor + s/2)^2);
thetaprimeplus = atan((rnor + s/2)/f);
thetaprimeminus = atan((rnor-s/2)/f);

rprimeminusmag = sqrt(f^2 + (rnor -s/2)^2);



Summedx = 0;
Summedy = 0;
Summedz = 0;


for phiprimeplus = 2*pi/16:(2*pi/16):32*pi/16
    

[myEfieldx, myEfieldy, myEfieldz] = channelefield(rprimeplusmag, thetaprimeplus, phiprimeplus, rprimeminusmag, thetaprimeminus, phiprimeplus-pi/8, r, theta, phi, choice);

Summedx = Summedx + myEfieldx; 
Summedy = Summedy + myEfieldy;
Summedz = Summedz + myEfieldz;
end

Efieldx = Summedx;
Efieldy = Summedy;
Efieldz = Summedz;

