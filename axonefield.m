function [Efieldx, Efieldy, Efieldz] = axonefield(numnodes, r, theta, phi, choice, offsetrmag, offsettheta)

internodallength = 50e-6*100;%50e-6*100;

rnor = 0.5e-6*1e2; % radius of node of ranvier in cgs units
d = 1e-9*1e2; %thickness of cell membrane at node in cgs units
g = rnor/1000000; % arbitrary dumbell radius in cgs units
s = d+(3*g/4); % separation of dipolar charges in cgs units

Summedx = 0;
Summedy = 0;
Summedz = 0;

for f = internodallength:internodallength:numnodes*internodallength
    thetaprime = atan(rnor/f)+offsettheta;
    rprime = sqrt(f^2 + rnor^2)+offsetrmag;
    
    
    
    [Efieldx, Efieldy, Efieldz] = ringefield(rprime, thetaprime, r, theta, phi, choice);
    
    %if f == 5*internodallength % there is only one node in an axon of ''length'' 10 nodes
    Summedx = Summedx + Efieldx;
    Summedy = Summedy  + Efieldy;
    Summedz = Summedz + Efieldz;
    %end
end

Efieldx = Summedx;
Efieldy = Summedy;
Efieldz = Summedz;

