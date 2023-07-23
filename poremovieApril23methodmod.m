%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Aman Chawla%
%% Rev2: March 31, 2016%
%% Rev: April 16, 2016%
%% Rev4: April 23, 2016: Extension to PWM-PPM %%
%%%%%%%%%%%%%%%%%%%%%%%

%%%% NOTE
%%% This file is based on Figure 1 of the following reference: 
%%%% Kallen, Roland G., Sidney A. Cohen, and Robert L. Barchi. 
%%% "Structure, function and expression of voltage-dependent sodium channels." 
%%% Molecular neurobiology 7.3-4 (1993): 383-428.
%%% This figure shows that on application of a command voltage, the patch
%%% which contains a single sodium channel, opens at a random
%%% post-initiation time (Pulse Position Modulation) and for a random duration
%%% (Pulse Width Modulation).
%% We incorporate this by positioning the trapezoid to have a random start time
%% and a random width. 

%% This file creates a movie that shows the electric field emanating from a node of ranvier
%% of negligible extension (in the axial z direction), inner diameter 1 micron,
%% and outer diameter 2 microns. The field is shown for upto 3 microns along the axon (in the 
%% z direction) and between 2.5 and 5.5 microns perpendicular distance (rho)
%% to the surface of the axonal membrane. With regards to time, the movie depicts 
%% the situation in which the ion channel time profile follows a trapezoidal opening
%% and closing pattern.

function returncellarray = poremovieApril23methodmod(start, width, trialnum,b, loops,fieldzero, conductivity)

% This program gives the field due to a ring at z' = 0. To find the field
% due to the entire node which goes from z' = 0 to z' = L, we simply add the
% fields at z, z-L/10, z-2*L/10, z-3*L/10, ..., z-L, for 10 rings. 


disp(start)
disp(width)
% We will plot the \hat{z} and \hat{rho} components of the E-field in two
% dimensions. Since they are orthogonal components, we can just use regular
% two dimensional vector field plotting methods. Since the electric field
% has no \hat{phi} component, we can ignore that component completely.
% Further, since the field is (rotationally) symmetric about the z-axis, we will just plot
% it in the y-z plane. Assuming the infinitesimal ring of current to be
% placed on the x-y plane. We will not vary phi. 

% Looking in a two dimensional picture, z will act as our y-coordinate and
% rho will act as our x coordinate. 
xstepsize = 0.1e-6;
ystepsize = 0.1e-6;
 x=2.5e-6:xstepsize:5.5e-6; % this is the cylindrical rho coordinate range, from the surface 
 % of the node of ranvier at 2.5 microns to about 5.5 microns away from the
 % surface. [[Field points]]
 y=0.1e-6:ystepsize:3.1e-6; % this is the cylindrical z coordinate range, which goes 
 % from 0.1 micro meters to about 3 microns along the length of the axon.
 % [[Field points]]
        
 mylength = length(x);
for i = 1:mylength
    for j = 1:mylength
        grid{i,j}  = {x(i),y(j)};
    end
end
% Extract the x-coordinates
for i = 1:mylength
    for j = 1:mylength
        px(i,j) = grid{i,j}(1);
    end
end
% Extract the y-coordinates
for i = 1:mylength
    for j = 1:mylength
        py(i,j) = grid{i,j}(2);
    end
end
pxnum = cell2mat(px);
pynum = cell2mat(py);
% conductivity = 1; %conductivity of extracellular medium Siemens per meter.

rho1prime = 1e-6; % lowu, lowu2, in microns - where the current source is at its peak strength.

Inot = 2e-12/(2*pi*rho1prime*0.1e-6); % peak current strength at surface of node is 2 pico Amperes. Thickness of 
% nodal cylindrical patch that acts as the source is taken to be 0.1
% microns roughly. radius is rho1prime. Thus Inot is the peak current per
% unit area emanating from the source cylindrical annulus at the node.
rho2prime = 2e-6; % highu, highu2, in microns - where the current source zeros out.
%For each (x,y) pair in the above matrix, compute the two integrals and
%store them as the px and py, respectively rho and z integrals. 
numchan = 20; %number of ion channels in a ring = 20; assume half of them are in the half ring. 
numchan = numchan/2;% just for the half-ring (line of sight influence assumption)
sourcelowv = [0:9]*pi/numchan; % single ion channel position in ring
sourcehighv = sourcelowv+(pi/numchan);% phi-width of the channel is added to get the upper limit.
returncellarray = cell(1,3);
field = cell(1,loops);
for j = 1:loops
    field{1,j} = 0;
end

for channelnum = 0:9

for i = 1:mylength
    for j = 1:mylength
        [Q, Q2] = zandrhointegralsmod(sourcelowv(channelnum+1), sourcehighv(channelnum+1),pxnum(i,j),pi/2,pynum(i,j));
        p{i,j} = {Q2*((Inot/(rho2prime-rho1prime))/(4*pi*conductivity)),Q*(pynum(i,j)*(Inot/(rho2prime-rho1prime))/(4*pi*conductivity))};
    end
end

% Extract the x-coordinates of the E-field vector, in units of Volt per
% meter (since all units are in S.I.)
for i = 1:mylength
    for j = 1:mylength
        vecx(i,j) = p{i,j}(1);
    end
end
% Extract the y-coordinates of the vector
for i = 1:mylength
    for j = 1:mylength
        vecy(i,j) = p{i,j}(2);
    end
end
vecxnum = cell2mat(vecx);
vecynum = cell2mat(vecy);

%axis tight manual;
ax.NextPlot = 'replaceChildren';
deltat = b/loops; % size of time steps


F(loops) = struct('cdata', [], 'colormap', []);
x = cell(1,loops);
y = cell(1,loops);

for j = 1:loops
        x{1,j} = pxnum;
        y{1,j} = pynum;
        field{1,j} = field{1,j} + randomtrapezoid(fieldzero,start, width,1,b,b*j/loops).*sqrt(vecxnum.^2 + vecynum.^2);
end

returncellarray{1,1} = x;
returncellarray{1,2} = y;



end
returncellarray{1,3} = field;