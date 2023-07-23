%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1 Axon Optimized Program %%
%%% Author: Aman Chawla      %%
%%% Date: October 31, 2015   %%
%%% Rev: October 31, 2015    %%
%%% Rev2: Nov 12, 2015 
%%% Rev3: Dec 3, 2015       %%
%%% Rev4: Dec 10, 2015  
%%% Rev5: Dec 27, 2015    
%%% Tweak1: April 9, 2016  
%%% 3D: Dec 24, 2016       %%
%%% Rev8: May 29, 2020: Added Conduction Velocity computation correctly %%
%%% Rev9: May 31, 2020: Study of effect of removing ephaptic coupling, on
%%% conduction velocity reduction/non-increasing nature. %%
%%% May 1, 2022: Goal is to see g-wave + ion channel noise together %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; 
clear all;
clf;
clc;
%% Note 2: This file adds conduction  velocity computation to mainaxonApril29.m
%% which has the below description.
%%% Note: This file calls the corrected activation update (which has the
%%% correct betah formula from Frankenhaeuser 1960) and also has a larger
%%% step size in time. With these two corrections synchronization is
%%% observed as expected. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  LTminus1 P_Na P_p P_K g_l

M = 100;
N =30; % Number of axons in bundle; multiple of 3 for ease of entry of W matrix as a repmat of N=3 case; N < M
d = 1e-3; %1e-3; %1e-3; % axon diameter (in centimeters), which is therefore 1e-5 meters.
ell = 2.5e-4; %2.5e-4; % length of Ranvier node (in centimeters.See Reutskiy Table 2.)
L = .2; % Original value: 0.2. length of a myelinated internodal segment in centi meters
delX = L/10; % Original value: L/10. length of a mesh point
delT = 2e-6; %2e-6; % length of a time unit in seconds
timesteps = 500; % number of time steps. duration = timesteps*delT
endlength = 250;%840; % total number of positions
temperature = 20;% temperature in degrees Celsius See Frankenhaeuser Huxley 1960
nodalspacing = 20;
Ranvier_array = 20:nodalspacing:endlength-20; % positions of Ranvier nodes
% For A only 
injectnumberarraya = [1,1,1,1,2,3,3,ones(1,N-15+8)];

% For V only

injectnumberarrayv = [20,21,22,23,24,25,27,28, 20*ones(1,N-15+7)];

% For AV -- joint stimulation of A and V
injectnumberarrayav = [1,1,1,1,2,20,21, 22, 23, 24, 25, ones(1,N-15),27,28,3,3];

% Equal angles case
theta = 0; %degrees
% oneightyminustheta = 180-theta;
oldcostheta = cosd(theta);
h = 0.1; % eventually time-varying
gcostheta = (1-h)*cosd(180-theta);
% olddelT = delT/cosd(theta); % why? figure out whether this is done so as to maintain stability?
%delT = delT/gcostheta;
%delT = olddelT;


ratio = .52;%0;%.5e1 % weakening the coupling to see the effect. %for very small ratios (1e-3) there is no interaction depolarization. for 1e-1 we begin to see interaction depolarization and for 1e0 there is full interaction depolarization
%% so we should stay less than 1. for 2e-1, we begin to see that the second axon, which is stimulated later, shows a faster AP condvel.
%% To see this more clearly, we reduce the nodal spacing and indeed see the effect more clearly. Next, we go to a ratio of 3e-1, in hopes of 
%% strengthening the velocity change and we do see the strengthening. So we next increase the simulation duration and axon lengths to see a better
%% effect. It seems like an edge effect (right edge). But there is a clear hyperpolarization on the second axon. So let us strengthen this hyperpolarization.
%% This may have an effect on the speed on the second one. Next we will try changing the fiber resistance itself. It is presently 1e5. We lower it to 1e4.
%% This has the effect of no propagation on either axon. So we raise it 1e6. We see an AP on the unstimulated axon before it is stimulated. So we come back to 1e5.
%% Next we will also try changing N to 3. We come back to N = 2. We check the delaydelta now. It is 50. We will reduce it first to 30. We raise it to 90. We raise
%% the ratio to 2.5e-1 from 2e-1. We come back to 2e-1.
rf = 1.27e8/cosd(theta) %1e7; %fiber resistance per unit length; Ohm/cm (Waxman); see Snider's page. rf must be less than 1e8.
% rf=1e4, r0 = 5e6 works with N=1. This is actually a two degree of freedom
% holonomic
% control problem. The control space is different from
% the physical or 3D space. Probably can be solved analytically. Change
% them alternately (ratio and rf).
r0 = ratio*rf;%0.5e3;%e8; %0.5e8; %r0 = 5*1e6 for N=2 with rf = 1e8 works; %  medium resistance per unit length; working value for coupling to be observed: 1e8
r0 = r0; % To account for the geometry.
cm =  1.87e-11/cosd(theta); % F/cm (both)
gm = 5.6e-9/cosd(theta); % per Ohm-cm (both)
cnd = 3.14e-9; % F/cm (from Waxman - Reutskiy has wrong dimensions)
gnd = 0;  % Not used in Reutskiy, but included in code here
InjectionDuration = 1e-5; % duration of current injected in seconds.
injectnumber1 = 1; % axon number to be injected
stimpos1 = 20; % nodal position to be injected.
delaydelta = 0; % delay between injections. 
injectnumber2 =1;
stimpos2 = 20;
injectnumber3 = 1;
stimpos3 = 20;
delaydelta2 = 0;
stimpos4 = 20;
delaydelta3 = 0;
injectnumber4 = 1;
% injectnumberarray = [1,1,1,1,2,20,21, 22, 23, 24, 25, ones(1,N-15),27,28,3,3];
stimposarray = [20,20,20,20,20.*ones(1,N-4)];
delaydeltarray = zeros(1,N);
Eonoff = 0; % State of applied field. 0: off, 1: on
Estrength = 2.5e-11; % Strength of applied field
tic;
diagterm = 1;
offdiagterm = 0.89625;
%K = [diagterm offdiagterm offdiagterm; offdiagterm diagterm offdiagterm; offdiagterm offdiagterm diagterm];
%K = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%K = [0.95 1 1; 1 0.95 1; 1 1 0.95];
% nontrivial geometry with Kdiag=1
%K = [1, .2, 1; .2, 1, .9; 1, .9, 1];
% The above still has 1 to 3 entry as 1, but it should be strictly less
% than 1. 
%K = [1, .5, 0.99; .5, 1, .5; 0.99, .5, 1];
%K1 = [1, 0.4, 0.35, 0.3; 0.4, 1, 0.4, 0.35; 0.35, 0.4, 1, 0.35; 0.3, 0.35, 0.35, 1];
%K = K1;
K2 = [1, 0.8, 0.5, 0.6; 0.8, 1, 0.8, 0.8; 0.5, 0.8, 1, 0.6; 0.6, 0.8, 0.6, 1];
K2prime = [1, 0.99, 0.5, 0.6; 0.99, 1, 0.99, 0.99; 0.5, 0.99, 1, 0.6; 0.6, 0.99, 0.6, 1];
K = K2prime;
% The below are the dimensionless distances from axon i to axon j The distances must obey the
% triangle inequality.
d11_1 = 0;
d12_1 = 0.12108;
d13_1 = 0.0101;
d23_1 = 0.111;
d22_1 = 0;
d33_1 = 0;
d21_1 = d12_1;
d31_1 = d13_1;
d32_1 = d23_1;
p1 = -1; % This is the power or drop-off factor.
% The below is the resultant K-matrix.
% K3_1 = [(1+d11_1)^(p1), (1+d12_1)^(p1), (1+d13_1)^(p1); (1+d12_1)^(p1), (1+d22_1)^(p1), (1+d23_1)^(p1); (1+d13_1)^(p1), (1+d23_1)^(p1), (1+d33_1)^(p1)];
% K3_1 = repmat(K3_1,N/3);

%
%K3_1 = [1, 1, 1; 1, 1, 1; 1, 1, 1];
%Knew = [1, .892, .99; .892, 1, .9; .99, .9, 1];
d12_2 = .1;
d13_2 = .2;
d23_2 = .25;
p2 = 4;
K3_2 = [1, d12_2^(p2), d13_2^(p2); d12_2^(p2), 1, d23_2^(p2); d13_2^(p2), d23_2^(p2), 1];
%%No coupling K4 matrix
K4 = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% modified nontrivial geometry (reflection of below K)
%K = [0, 1.5, .5; 1.5, 0, 1; .5, 1,  0];
%K = [0, .5, 1.5; .5, 0, 1; 1.5, 1, 0]; %nontrivial triangular geometry
%K = [0, 1, 1; 1, 0, 1; 1, 1, 0]; % trivial triangular geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXECUTE CODE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K3_1 = Wmatrixcreator(N,M);

K3_1a = [
    1.0000    0.0161    0.0118    0.0707    0.0092    0.0252    0.0228    0.0135    0.0108    0.0180    0.0165    0.0325    0.0078    0.0181    0.0115    0.0105    0.0120
    0.0161    1.0000    0.0383    0.0200    0.0178    0.0409    0.0261    0.0324    0.0285    0.1239    0.0182    0.0301    0.0137    0.0152    0.0216    0.0179    0.0269
    0.0118    0.0383    1.0000    0.0140    0.0176    0.0218    0.0158    0.0438    0.0994    0.0321    0.0128    0.0183    0.0146    0.0111    0.0314    0.0149    0.0207
    0.0707    0.0200    0.0140    1.0000    0.0101    0.0371    0.0261    0.0164    0.0125    0.0232    0.0172    0.0556    0.0085    0.0177    0.0135    0.0114    0.0134
    0.0092    0.0178    0.0176    0.0101    1.0000    0.0131    0.0148    0.0130    0.0166    0.0160    0.0153    0.0119    0.0476    0.0123    0.0114    0.0384    0.0388
    0.0252    0.0409    0.0218    0.0371    0.0131    1.0000    0.0314    0.0257    0.0184    0.0575    0.0187    0.1006    0.0105    0.0169    0.0186    0.0142    0.0185
    0.0228    0.0261    0.0158    0.0261    0.0148    0.0314    1.0000    0.0154    0.0139    0.0282    0.0436    0.0316    0.0114    0.0341    0.0124    0.0192    0.0232
    0.0135    0.0324    0.0438    0.0164    0.0130    0.0257    0.0154    1.0000    0.0385    0.0323    0.0119    0.0217    0.0111    0.0107    0.0595    0.0120    0.0156
    0.0108    0.0285    0.0994    0.0125    0.0166    0.0184    0.0139    0.0385    1.0000    0.0250    0.0115    0.0158    0.0145    0.0101    0.0333    0.0137    0.0183
    0.0180    0.1239    0.0321    0.0232    0.0160    0.0575    0.0282    0.0323    0.0250    1.0000    0.0185    0.0381    0.0125    0.0157    0.0214    0.0166    0.0238
    0.0165    0.0182    0.0128    0.0172    0.0153    0.0187    0.0436    0.0119    0.0115    0.0185    1.0000    0.0187    0.0118    0.0595    0.0100    0.0230    0.0230
    0.0325    0.0301    0.0183    0.0556    0.0119    0.1006    0.0316    0.0217    0.0158    0.0381    0.0187    1.0000    0.0098    0.0177    0.0165    0.0132    0.0165
    0.0078    0.0137    0.0146    0.0085    0.0476    0.0105    0.0114    0.0111    0.0145    0.0125    0.0118    0.0098    1.0000    0.0100    0.0102    0.0230    0.0219
    0.0181    0.0152    0.0111    0.0177    0.0123    0.0169    0.0341    0.0107    0.0101    0.0157    0.0595    0.0177    0.0100    1.0000    0.0092    0.0170    0.0169
    0.0115    0.0216    0.0314    0.0135    0.0114    0.0186    0.0124    0.0595    0.0333    0.0214    0.0100    0.0165    0.0102    0.0092    1.0000    0.0103    0.0129
    0.0105    0.0179    0.0149    0.0114    0.0384    0.0142    0.0192    0.0120    0.0137    0.0166    0.0230    0.0132    0.0230    0.0170    0.0103    1.0000    0.0492
    0.0120    0.0269    0.0207    0.0134    0.0388    0.0185    0.0232    0.0156    0.0183    0.0238    0.0230    0.0165    0.0219    0.0169    0.0129    0.0492    1.0000
    0.0124    0.0420    0.1667    0.0148    0.0166    0.0238    0.0163    0.0548    0.0748    0.0357    0.0130    0.0197    0.0138    0.0113    0.0347    0.0145    0.0200
    0.0195    0.0419    0.0207    0.0237    0.0173    0.0397    0.0624    0.0191    0.0175    0.0449    0.0304    0.0340    0.0129    0.0232    0.0147    0.0211    0.0299
    0.0126    0.0438    0.0324    0.0147    0.0292    0.0226    0.0225    0.0213    0.0264    0.0342    0.0189    0.0191    0.0194    0.0149    0.0167    0.0265    0.0541
    0.0221    0.0162    0.0115    0.0214    0.0119    0.0193    0.0408    0.0114    0.0104    0.0171    0.0465    0.0207    0.0096    0.0905    0.0097    0.0157    0.0164
    0.1208    0.0175    0.0127    0.1218    0.0095    0.0293    0.0230    0.0148    0.0115    0.0198    0.0161    0.0397    0.0080    0.0171    0.0125    0.0107    0.0123
    0.0144    0.1208    0.0514    0.0175    0.0190    0.0320    0.0224    0.0348    0.0351    0.0671    0.0167    0.0249    0.0146    0.0139    0.0231    0.0180    0.0274
    0.0089    0.0156    0.0150    0.0096    0.0922    0.0120    0.0141    0.0116    0.0143    0.0143    0.0153    0.0111    0.0471    0.0123    0.0103    0.0426    0.0328
    0.0096    0.0167    0.0149    0.0104    0.0553    0.0130    0.0163    0.0117    0.0139    0.0154    0.0186    0.0121    0.0298    0.0145    0.0102    0.0885    0.0415
    0.0143    0.0152    0.0114    0.0145    0.0146    0.0154    0.0290    0.0105    0.0104    0.0153    0.0790    0.0154    0.0115    0.0541    0.0090    0.0226    0.0203
    0.0095    0.0156    0.0139    0.0102    0.0445    0.0125    0.0159    0.0111    0.0130    0.0145    0.0186    0.0116    0.0277    0.0146    0.0097    0.0846    0.0357
    0.0128    0.0148    0.0114    0.0131    0.0161    0.0143    0.0248    0.0103    0.0104    0.0147    0.0526    0.0141    0.0125    0.0358    0.0089    0.0267    0.0219
    0.0107    0.0231    0.0203    0.0119    0.0630    0.0160    0.0189    0.0148    0.0183    0.0204    0.0191    0.0144    0.0280    0.0147    0.0125    0.0495    0.0909
    0.0180    0.0345    0.0191    0.0210    0.0185    0.0309    0.0653    0.0171    0.0164    0.0349    0.0369    0.0277    0.0135    0.0256    0.0135    0.0243    0.0341]

%   Columns 18 through 30
K3_1b = [
    0.0124    0.0195    0.0126    0.0221    0.1208    0.0144    0.0089    0.0096    0.0143    0.0095    0.0128    0.0107    0.0180
    0.0420    0.0419    0.0438    0.0162    0.0175    0.1208    0.0156    0.0167    0.0152    0.0156    0.0148    0.0231    0.0345
    0.1667    0.0207    0.0324    0.0115    0.0127    0.0514    0.0150    0.0149    0.0114    0.0139    0.0114    0.0203    0.0191
    0.0148    0.0237    0.0147    0.0214    0.1218    0.0175    0.0096    0.0104    0.0145    0.0102    0.0131    0.0119    0.0210
    0.0166    0.0173    0.0292    0.0119    0.0095    0.0190    0.0922    0.0553    0.0146    0.0445    0.0161    0.0630    0.0185
    0.0238    0.0397    0.0226    0.0193    0.0293    0.0320    0.0120    0.0130    0.0154    0.0125    0.0143    0.0160    0.0309
    0.0163    0.0624    0.0225    0.0408    0.0230    0.0224    0.0141    0.0163    0.0290    0.0159    0.0248    0.0189    0.0653
    0.0548    0.0191    0.0213    0.0114    0.0148    0.0348    0.0116    0.0117    0.0105    0.0111    0.0103    0.0148    0.0171
    0.0748    0.0175    0.0264    0.0104    0.0115    0.0351    0.0143    0.0139    0.0104    0.0130    0.0104    0.0183    0.0164
    0.0357    0.0449    0.0342    0.0171    0.0198    0.0671    0.0143    0.0154    0.0153    0.0145    0.0147    0.0204    0.0349
    0.0130    0.0304    0.0189    0.0465    0.0161    0.0167    0.0153    0.0186    0.0790    0.0186    0.0526    0.0191    0.0369
    0.0197    0.0340    0.0191    0.0207    0.0397    0.0249    0.0111    0.0121    0.0154    0.0116    0.0141    0.0144    0.0277
    0.0138    0.0129    0.0194    0.0096    0.0080    0.0146    0.0471    0.0298    0.0115    0.0277    0.0125    0.0280    0.0135
    0.0113    0.0232    0.0149    0.0905    0.0171    0.0139    0.0123    0.0145    0.0541    0.0146    0.0358    0.0147    0.0256
    0.0347    0.0147    0.0167    0.0097    0.0125    0.0231    0.0103    0.0102    0.0090    0.0097    0.0089    0.0125    0.0135
    0.0145    0.0211    0.0265    0.0157    0.0107    0.0180    0.0426    0.0885    0.0226    0.0846    0.0267    0.0495    0.0243
    0.0200    0.0299    0.0541    0.0164    0.0123    0.0274    0.0328    0.0415    0.0203    0.0357    0.0219    0.0909    0.0341
    1.0000    0.0215    0.0309    0.0118    0.0134    0.0556    0.0143    0.0143    0.0115    0.0134    0.0114    0.0193    0.0196
    0.0215    1.0000    0.0327    0.0253    0.0205    0.0335    0.0159    0.0182    0.0227    0.0174    0.0211    0.0232    0.1208
    0.0309    0.0327    1.0000    0.0151    0.0133    0.0499    0.0236    0.0252    0.0164    0.0226    0.0168    0.0454    0.0329
    0.0118    0.0253    0.0151    1.0000    0.0206    0.0146    0.0117    0.0135    0.0367    0.0135    0.0274    0.0142    0.0270
    0.0134    0.0205    0.0133    0.0206    1.0000    0.0156    0.0090    0.0098    0.0139    0.0096    0.0125    0.0110    0.0186
    0.0556    0.0335    0.0499    0.0146    0.0156    1.0000    0.0164    0.0171    0.0142    0.0159    0.0140    0.0244    0.0293
    0.0143    0.0159    0.0236    0.0117    0.0090    0.0164    1.0000    0.0748    0.0150    0.0624    0.0168    0.0457    0.0172
    0.0143    0.0182    0.0252    0.0135    0.0098    0.0171    0.0748    1.0000    0.0184    0.1667    0.0211    0.0526    0.0203
    0.0115    0.0227    0.0164    0.0367    0.0139    0.0142    0.0150    0.0184    1.0000    0.0188    0.0958    0.0175    0.0264
    0.0134    0.0174    0.0226    0.0135    0.0096    0.0159    0.0624    0.1667    0.0188    1.0000    0.0219    0.0419    0.0193
    0.0114    0.0211    0.0168    0.0274    0.0125    0.0140    0.0168    0.0211    0.0958    0.0219    1.0000    0.0191    0.0247
    0.0193    0.0232    0.0454    0.0142    0.0110    0.0244    0.0457    0.0526    0.0175    0.0419    0.0191    1.0000    0.0255
    0.0196    0.1208    0.0329    0.0270    0.0186    0.0293    0.0172    0.0203    0.0264    0.0193    0.0247    0.0255    1.0000]
myK3_1 = [K3_1a K3_1b];
% call the crank nicolson N axon method to compute the axonal voltages.
axon = packagedCNNDec24_2016v6(N, d, ell, L, delT, timesteps, endlength, temperature, Ranvier_array, r0, cm, gm, cnd, gnd, rf, InjectionDuration, injectnumberarraya, stimposarray, delaydeltarray, Eonoff, Estrength, myK3_1, theta);
Timespent = toc
LTminus1 = timesteps;
% 
% % for count = 1:N
% %     [a(count), b(count)] = velocitycomputer(count,axon);
% % end
% % if (a(1) ~= a(2))|| (b(1) ~= b(2))
% %    disp('Success!');
% % elseif (a(1) == a(2))&&(b(1)==b(2))
% %     disp('Failure');
% % end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SAVE OUTPUT %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename = strcat('numax', num2str(N), 'timedur', num2str(timesteps), 'axonlen', num2str(endlength), 'interres', num2str(r0), 'stimpos', num2str(stimpos1), 'date', date);
% 
% save(filename);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT OUTPUT %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
figure(1)

    for plotnumber = 1:N
        for i = 1:endlength
         subplot(N,1,plotnumber)
                outputs = screation(1);

    if ismember(i,Ranvier_array)
        datamy = axon(plotnumber).V(i,1:LTminus1);
        spikeposonly  = find(datamy>0.1);
        spikeval = datamy(spikeposonly);
        yvals = embed(spikeposonly,zeros(1,LTminus1),spikeval);
        
plot([1:LTminus1].*delT, datamy, outputs);
axis([[1, LTminus1]*delT,-.1,.15]);
hold on;
grid on;
    end
    end
    end

    % plot the CAP

%     figure(2)
%     datamy = zeros(1,LTminus1)
%     for plotnumber = 1:N
%         for i = 1:endlength
% %          subplot(N,1,plotnumber)
%                 outputs = screation(1);
% 
%     if ismember(i,Ranvier_array)
%         datamy = datamy+axon(plotnumber).V(i,1:LTminus1);
% %         spikeposonly  = find(datamy>0.1);
% %         spikeval = datamy(spikeposonly);
% %         yvals = embed(spikeposonly,zeros(1,LTminus1),spikeval);
%         
% 
%     end
%     end
%     end
% 
%     plot([1:LTminus1].*delT, datamy, outputs);
% axis([[1, LTminus1]*delT,-.1,.15]);
% 
% grid on;
% 
% title('Voltage vs Time at 1st six Ranvier nodes of axon 1 and 2');
% 
% grid on;
% 
% 
% %%%%COMPUTE CV %%%%%
% %%% CV is computed at two places: initial CV near the injection point (say nodes 40 and 60) and
% %%% final CV near the end of the axon (say nodes 240 and 260). This would allow observation of any
% %%% change in CV due to ephaptic interaction. The CV is computed by finding
% %%% the time at which the voltage crosses the 100mV mark at node 40 and
% %%% subtracting that from the time at which the voltage crosses the 100mV
% %%% mark at node 60 and finally inverting this time difference and
% %%% multiplying by the distance (20*delX). This is repeated for nodes 240
% %%% and 260. This is done for both the fibers. 
% initialnode1 = 40;
% initialnode2 = 60;
% initialnodalgap = initialnode2-initialnode1;
% finalnode1 = 120;
% finalnode2 = 140;
% finalnodalgap = finalnode2 - finalnode1;
% 
% %% CV initial on fiber 1
% temp1 = axon(1).V(initialnode1, 1:LTminus1) > .1;
% difftemp1 =  diff(temp1);
% for i =  1:length(difftemp1)
%     if difftemp1(i) == 1
%         storedtime1 = i+1;
%     end
% end
% 
% temp2 = axon(1).V(initialnode2, 1:LTminus1) > .1;
% difftemp2 =  diff(temp2);
% for i =  1:length(difftemp2)
%     if difftemp2(i) == 1
%         storedtime2 = i+1;
%     end
% end
% 
% timegap = delT.*(storedtime2-storedtime1);
% inversetimegap = 1/timegap;
% CVinitfiber1 = initialnodalgap*delX*inversetimegap
% CVinitfiber1mps = CVinitfiber1/100
% 
% %% CV final on fiber 1
% temp3 = axon(1).V(finalnode1, 1:LTminus1) > .1;
% difftemp3 =  diff(temp3);
% for i =  1:length(difftemp3)
%     if difftemp3(i) == 1
%         storedtime3 = i;
%     end
% end
% 
% temp4 = axon(1).V(finalnode2, 1:LTminus1) > .1;
% difftemp4 =  diff(temp4);
% for i =  1:length(difftemp4)
%     if difftemp4(i) == 1
%         storedtime4 = i;
%     end
% end
% 
% timegap2 = delT.*(storedtime4-storedtime3);
% inversetimegap2 = 1/timegap2;
% CVfinalfiber1 = finalnodalgap*delX*inversetimegap2 % This would be in centimeters per second
% CVfinalfiber1mps = CVfinalfiber1/100
% %% CV initial on fiber 2
% temp5 = axon(2).V(initialnode1, 1:LTminus1) > .1;
% difftemp5 =  diff(temp5);
% for i =  1:length(difftemp5)
%     if difftemp5(i) == 1
%         storedtime5 = i;
%     end
% end
% 
% temp6 = axon(2).V(initialnode2, 1:LTminus1) > .1;
% difftemp6 =  diff(temp6);
% for i =  1:length(difftemp6)
%     if difftemp6(i) == 1
%         storedtime6 = i;
%     end
% end
% 
% timegap3 = delT.*(storedtime6-storedtime5);
% inversetimegap3 = 1/timegap3;
% CVinitfiber2 = initialnodalgap*delX*inversetimegap3
% CVinitfiber2mps = CVinitfiber2/100
% %% CV final on fiber 2
% temp7 = axon(2).V(finalnode1, 1:LTminus1) > .1;
% difftemp7 =  diff(temp7);
% for i =  1:length(difftemp7)
%     if difftemp7(i) == 1
%         storedtime7 = i;
%     end
% end
% 
% temp8 = axon(2).V(finalnode2, 1:LTminus1) > .1;
% difftemp8 =  diff(temp8);
% for i =  1:length(difftemp8)
%     if difftemp8(i) == 1
%         storedtime8 = i;
%     end
% end
% 
% timegap4 = delT.*(storedtime8-storedtime7);
% inversetimegap4 = 1/timegap4;
% CVfinalfiber2 = finalnodalgap*delX*inversetimegap4
% CVfinalfiber2mps = CVfinalfiber2/100
% 
% %% CV initial on fiber 3
% temp9 = axon(3).V(initialnode1, 1:LTminus1) > .1;
% difftemp9 =  diff(temp9);
% for i =  1:length(difftemp9)
%     if difftemp9(i) == 1
%         storedtime9 = i;
%     end
% end
% 
% temp10 = axon(3).V(initialnode2, 1:LTminus1) > .1;
% difftemp10 =  diff(temp10);
% for i =  1:length(difftemp10)
%     if difftemp10(i) == 1
%         storedtime10 = i;
%     end
% end
% 
% timegap5 = delT.*(storedtime10-storedtime9);
% inversetimegap5 = 1/timegap5;
% CVinitfiber3 = initialnodalgap*delX*inversetimegap5
% CVinitfiber3mps = CVinitfiber3/100
% %% CV final on fiber 2
% temp11 = axon(3).V(finalnode1, 1:LTminus1) > .1;
% difftemp11 =  diff(temp11);
% for i =  1:length(difftemp11)
%     if difftemp11(i) == 1
%         storedtime11 = i;
%     end
% end
% 
% temp12 = axon(3).V(finalnode2, 1:LTminus1) > .1;
% difftemp12 =  diff(temp12);
% for i =  1:length(difftemp12)
%     if difftemp12(i) == 1
%         storedtime12 = i;
%     end
% end
% 
% timegap6 = delT.*(storedtime12-storedtime11);
% inversetimegap6 = 1/timegap6;
% CVfinalfiber3 = finalnodalgap*delX*inversetimegap6
% CVfinalfiber3mps = CVfinalfiber3/100
% 
% 
% %
% % 
% % figure(2)
% % for time = 1:(timesteps-1)
% %     for plotnumber = 1:N
% %         subplot(N,1,plotnumber)
% %         plot(1:endlength, axon(plotnumber).V(1:endlength,time));
% %    
% %         plot(1:endlength, axon(plotnumber).Vcoeff(1:endlength,time), 'r');
% %         string = strcat('Position along axon number ', num2str(plotnumber));
% %         xlabel(string);
% %         ylabel(time);
% %         title('Absolute voltage (blue) profile ');
% %         axis([1, endlength, -1, 1]); 
% %         drawnow
% %      
% %     end
% %     
% % end
% % 
% figure(3)
% mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(1).V(:,:));
% xlabel('Time in seconds');
% ylabel('Position in centimeters');
% zlabel('Voltage');
% 
% figure(4)
% mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(2).V(:,:));
% xlabel('Time in seconds');
% ylabel('Position in centimeters');
% zlabel('Voltage');
% 
% figure(5)
% mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(3).V(:,:));
% xlabel('Time in seconds');
% ylabel('Position in centimeters');
% zlabel('Voltage');
% 
% % figure(5)
% % mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(3).V(:,:));
% % xlabel('Time in seconds');
% % ylabel('Position in centimeters');
% % zlabel('Voltage');
% % 
% % 
% % % Ionic currents in milliAmps per cm^2 vs Time (msec)
% % farnode = 180;
% % LTminus1 = timesteps-1;
% % figure(4)
% % plot(1:LTminus1, axon(1).J_Na(farnode, 1:LTminus1),'b')
% % hold on;
% % plot(1:LTminus1, axon(1).J_K(farnode, 1:LTminus1),'r')
% % hold on;
% % plot(1:LTminus1, axon(1).J_p(farnode, 1:LTminus1),'k')
% % hold on;
% % plot(1:LTminus1, axon(1).J_l(farnode, 1:LTminus1),'c')
% % hold on;
% % xlabel('Time');
% % ylabel('Current value');
% % title('Various currents blue sodium, red potassium at position far from stimpos1');
% % 
% % figure(5)
% % plot(1:LTminus1, axon(1).J_Na(stimpos1, 1:LTminus1),'b')
% % hold on;
% % plot(1:LTminus1, axon(1).J_K(stimpos1, 1:LTminus1),'r')
% % hold on;
% % plot(1:LTminus1, axon(1).J_p(stimpos1, 1:LTminus1),'k')
% % hold on;
% % plot(1:LTminus1, axon(1).J_l(stimpos1, 1:LTminus1),'c')
% % hold on;
% % xlabel('Time');
% % ylabel('Current value');
% % title('Various currents blue sodium, red potassium at stimpos1');
% % 
% % %
% % 
% % figure(6)
% % subplot(2,3,1)
% % plot(1:LTminus1, axon(1).J_Na(farnode, 1:LTminus1),'b')
% % hold on;
% % plot(1:LTminus1, axon(1).J_K(farnode, 1:LTminus1),'r')
% % hold on;
% % plot(1:LTminus1, axon(1).J_p(farnode, 1:LTminus1),'k')
% % hold on;
% % plot(1:LTminus1, axon(1).J_l(farnode, 1:LTminus1),'c')
% % hold on;
% % xlabel('Time');
% % ylabel('Current value');
% % title('Various currents blue sodium, red potassium, cyan leakage, black p, at position far from stimpos1');
% % 
% 
% % Ionic current
% figure(9)
% plot((1:LTminus1-1).*delT, axon(1).J_Na(stimpos1, 1:(LTminus1-1)),'b')
% hold on;
% plot((1:LTminus1-1).*delT, axon(1).J_K(stimpos1, 1:(LTminus1-1)),'r')
% xlabel('Time');
% ylabel('Current magnitude');
% title('Ionic currents: blue sodium, red potassium at position of stimulation with time');
% figure(13)
% plot((1:LTminus1-1).*delT, axon(1).J_Na(stimpos1+40, 1:LTminus1-1),'b')
% hold on;
% plot((1:LTminus1-1).*delT, axon(1).J_K(stimpos1+40, 1:LTminus1-1),'r')
% xlabel('Time');
% ylabel('Current magnitude');
% title('Ionic currents: blue sodium, red potassium at position downstream from stimulation position with time');
% % 
% % subplot(2,3,3)
% % mesh(1:(timesteps), 1:endlength, axon(1).V(:,:));
% % xlabel('Time');
% % ylabel('Position');
% % zlabel('Voltage');
% % 
% figure(10)
% farnode = 160;
% plot(1:LTminus1, axon(1).m(farnode, 1:LTminus1),'b')
% hold on;
% plot(1:LTminus1, axon(1).n(farnode,1:LTminus1), 'r')
% hold on;
% plot(1:LTminus1, axon(1).h(farnode,1:LTminus1), 'k')
% hold on;
% plot(1:LTminus1, axon(1).p(farnode,1:LTminus1), 'c')
% hold on;
% xlabel('Time');
% ylabel('Gating variable magnitude');
% title('Gating variables with time at position 160 which is downaxon from the position of injection');
% % 
% figure(11)
% plot(1:LTminus1, axon(1).m(stimpos1, 1:LTminus1),'b')
% hold on;
% plot(1:LTminus1, axon(1).n(stimpos1,1:LTminus1), 'r')
% hold on;
% plot(1:LTminus1, axon(1).h(stimpos1,1:LTminus1), 'k')
% hold on;
% plot(1:LTminus1, axon(1).p(stimpos1,1:LTminus1), 'c')
% hold on;
% xlabel('Time');
% ylabel('Gating variable magnitude');
% title('Gating variables with time at stimulated node');
% % 
% % subplot(2,3,6)
% % plot(1:LTminus1, axon(1).J_ion(farnode, 1:LTminus1),'b');
% % xlabel('Time');
% % ylabel('Ionic current');
% % title('Ionic current at position far from stimpos1');
% % 
% disp(ratio);
% disp(rf);
