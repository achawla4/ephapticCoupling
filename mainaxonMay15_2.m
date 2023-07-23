%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1 Axon Optimized Program %%
%%% Author: Aman Chawla      %%
%%% Date: October 31, 2015   %%
%%% Rev: October 31, 2015    %%
%%% Rev2: Nov 12, 2015 
%%% Rev3: Dec 3, 2015       %%
%%% Rev4: Dec 10, 2015  
%%% Rev5: Dec 27, 2015    
%%% Tweak1: April 9, 2016   %%
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

N = 1; % Number of axons in bundle
d = 1e-3; %1e-3; %1e-3; % axon diameter (in centimeters), which is therefore 1e-5 meters.
ell = 2.5e-4; %2.5e-4; % length of Ranvier node (in centimeters.See Reutskiy Table 2.)
L = .2; % Original value: 0.2. length of a myelinated internodal segment in centi meters
delX = L/10; % Original value: L/10. length of a mesh point
delT = 2e-6; %2e-6; % length of a time unit in seconds
timesteps = 2000; % number of time steps. duration = timesteps*delT
endlength = 400;%840; % total number of positions
temperature = 29.3; %26.6;%20;% temperature in degrees Celsius See Frankenhaeuser Huxley 1960
nodalspacing = 10;
Ranvier_array = [20,30, 60:nodalspacing:endlength-20]; % positions of Ranvier nodes

ratio = 0.5;%0;%.5e1 % weakening the coupling to see the effect. %for very small ratios (1e-3) there is no interaction depolarization. for 1e-1 we begin to see interaction depolarization and for 1e0 there is full interaction depolarization
%% so we should stay less than 1. for 2e-1, we begin to see that the second axon, which is stimulated later, shows a faster AP condvel.
%% To see this more clearly, we reduce the nodal spacing and indeed see the effect more clearly. Next, we go to a ratio of 3e-1, in hopes of 
%% strengthening the velocity change and we do see the strengthening. So we next increase the simulation duration and axon lengths to see a better
%% effect. It seems like an edge effect (right edge). But there is a clear hyperpolarization on the second axon. So let us strengthen this hyperpolarization.
%% This may have an effect on the speed on the second one. Next we will try changing the fiber resistance itself. It is presently 1e5. We lower it to 1e4.
%% This has the effect of no propagation on either axon. So we raise it 1e6. We see an AP on the unstimulated axon before it is stimulated. So we come back to 1e5.
%% Next we will also try changing N to 3. We come back to N = 2. We check the delaydelta now. It is 50. We will reduce it first to 30. We raise it to 90. We raise
%% the ratio to 2.5e-1 from 2e-1. We come back to 2e-1.
rf = 1.27e8 %1e7; %fiber resistance per unit length; Ohm/cm (Waxman); see Snider's page. rf must be less than 1e8.
% rf=1e4, r0 = 5e6 works with N=1. This is actually a two degree of freedom
% holonomic
% control problem. The control space is different from
% the physical or 3D space. Probably can be solved analytically. Change
% them alternately (ratio and rf).
r0 = ratio*rf;%0.5e3;%e8; %0.5e8; %r0 = 5*1e6 for N=2 with rf = 1e8 works; % interfiber medium resistance per unit length; working value for coupling to be observed: 1e8
cm =  1.87e-11; % F/cm (both)
gm = 5.6e-9; % per Ohm-cm (both)
cnd = 3.14e-9; % F/cm (from Waxman - Reutskiy has wrong dimensions)
gnd = 0;  % Not used in Reutskiy, but included in code here
InjectionDuration = 1e-5; % duration of current injected in seconds.
injectnumber1 = 1; % axon number to be injected
stimpos1 = 20; % nodal position to be injected.
delaydelta = 50; % delay between injections. 
injectnumber2 = 1;
stimpos2 = 20;
injectnumber3 = 1;
stimpos3 = 20;
delaydelta2 = 0;
stimpos4 = 20;
delaydelta3 = 0;
injectnumber4 = 1;
Eonoff = 0; % State of applied field. 0: off, 1: on
Estrength = 2.5e-11; % Strength of applied field
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXECUTE CODE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call the crank nicolson N axon method to compute the axonal voltages.
axon = packagedCNNDec28(N, d, ell, L, delT, timesteps, endlength, temperature, Ranvier_array, r0, cm, gm, cnd, gnd, rf, InjectionDuration, injectnumber1, stimpos1, delaydelta, injectnumber2, stimpos2, delaydelta2, injectnumber3, stimpos3, delaydelta3, injectnumber4, stimpos4, Eonoff, Estrength);
Timespent = toc
LTminus1 = timesteps;

% for count = 1:N
%     [a(count), b(count)] = velocitycomputer(count,axon);
% end
% if (a(1) ~= a(2))|| (b(1) ~= b(2))
%    disp('Success!');
% elseif (a(1) == a(2))&&(b(1)==b(2))
%     disp('Failure');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE OUTPUT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat('numax', num2str(N), 'timedur', num2str(timesteps), 'axonlen', num2str(endlength), 'interres', num2str(r0), 'stimpos', num2str(stimpos1), 'date', date);

save(filename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OUTPUT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
figure(1)
for i = 1:endlength
    for plotnumber = 1:N
        subplot(N,1,plotnumber)
    if ismember(i,Ranvier_array)
       
plot([1:LTminus1].*delT, axon(plotnumber).V(i, 1:LTminus1))
axis([[1, LTminus1]*delT,-.1,.15]);
hold on;
grid on;
    end
    end
end

title('Voltage vs Time at 1st six Ranvier nodes of axon 1 and 2');

grid on;


% %%%%COMPUTE CV %%%%%
% %%% CV is computed at two places: initial CV near the injection point (say nodes 40 and 60) and
% %%% final CV near the end of the axon (say nodes 240 and 260). This would allow observation of any
% %%% change in CV due to ephaptic interaction. The CV is computed by finding
% %%% the time at which the voltage crosses the 100mV mark at node 40 and
% %%% subtracting that from the time at which the voltage crosses the 100mV
% %%% mark at node 60 and finally inverting this time difference and
% %%% multiplying by the distance (20*delX). This is repeated for nodes 240
% %%% and 260. This is done for both the fibers. 
% 
% %% CV initial on fiber 1
% temp = axon(1).V(40, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime1 = i;
%     end
% end
% 
% temp = axon(1).V(60, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime2 = i;
%     end
% end
% 
% timegap = delT.*(storedtime2-storedtime1);
% inversetimegap = 1/timegap;
% CVinitfiber1 = 20*delX*inversetimegap
% 
% %% CV final on fiber 1
% temp = axon(1).V(200, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime3 = i;
%     end
% end
% 
% temp = axon(1).V(240, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime4 = i;
%     end
% end
% 
% timegap2 = delT.*(storedtime2-storedtime1);
% inversetimegap2 = 1/timegap2;
% CVfinalfiber1 = 20*delX*inversetimegap2
% 
% %% CV initial on fiber 2
% temp = axon(2).V(40, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime5 = i;
%     end
% end
% 
% temp = axon(2).V(60, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime6 = i;
%     end
% end
% 
% timegap3 = delT.*(storedtime6-storedtime5);
% inversetimegap3 = 1/timegap3;
% CVinitfiber2 = 20*delX*inversetimegap3
% 
% %% CV final on fiber 2
% temp = axon(2).V(200, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime7 = i;
%     end
% end
% 
% temp = axon(2).V(240, 1:LTminus1) > .1;
% difftemp =  diff(temp);
% for i =  1:length(difftemp)
%     if difftemp(i) == 1
%         storedtime8 = i;
%     end
% end
% 
% timegap4 = delT.*(storedtime8-storedtime7);
% inversetimegap4 = 1/timegap4;
% CVfinalfiber3 = 20*delX*inversetimegap4


%
% 
% figure(2)
% for time = 1:(timesteps-1)
%     for plotnumber = 1:N
%         subplot(N,1,plotnumber)
%         plot(1:endlength, axon(plotnumber).V(1:endlength,time));
%    
%         plot(1:endlength, axon(plotnumber).Vcoeff(1:endlength,time), 'r');
%         string = strcat('Position along axon number ', num2str(plotnumber));
%         xlabel(string);
%         ylabel(time);
%         title('Absolute voltage (blue) profile ');
%         axis([1, endlength, -1, 1]); 
%         drawnow
%      
%     end
%     
% end
% 
figure(3)
mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(1).V(:,:));
xlabel('Time in seconds');
ylabel('Position in centimeters');
zlabel('Voltage');

figure(4)
mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(2).V(:,:));
xlabel('Time in seconds');
ylabel('Position in centimeters');
zlabel('Voltage');

% figure(5)
% mesh([1:(LTminus1)].*delT, [1:endlength].*delX, axon(3).V(:,:));
% xlabel('Time in seconds');
% ylabel('Position in centimeters');
% zlabel('Voltage');
% 
% 
% % Ionic currents in milliAmps per cm^2 vs Time (msec)
% farnode = 180;
% LTminus1 = timesteps-1;
% figure(4)
% plot(1:LTminus1, axon(1).J_Na(farnode, 1:LTminus1),'b')
% hold on;
% plot(1:LTminus1, axon(1).J_K(farnode, 1:LTminus1),'r')
% hold on;
% plot(1:LTminus1, axon(1).J_p(farnode, 1:LTminus1),'k')
% hold on;
% plot(1:LTminus1, axon(1).J_l(farnode, 1:LTminus1),'c')
% hold on;
% xlabel('Time');
% ylabel('Current value');
% title('Various currents blue sodium, red potassium at position far from stimpos1');
% 
% figure(5)
% plot(1:LTminus1, axon(1).J_Na(stimpos1, 1:LTminus1),'b')
% hold on;
% plot(1:LTminus1, axon(1).J_K(stimpos1, 1:LTminus1),'r')
% hold on;
% plot(1:LTminus1, axon(1).J_p(stimpos1, 1:LTminus1),'k')
% hold on;
% plot(1:LTminus1, axon(1).J_l(stimpos1, 1:LTminus1),'c')
% hold on;
% xlabel('Time');
% ylabel('Current value');
% title('Various currents blue sodium, red potassium at stimpos1');
% 
% %
% 
% figure(6)
% subplot(2,3,1)
% plot(1:LTminus1, axon(1).J_Na(farnode, 1:LTminus1),'b')
% hold on;
% plot(1:LTminus1, axon(1).J_K(farnode, 1:LTminus1),'r')
% hold on;
% plot(1:LTminus1, axon(1).J_p(farnode, 1:LTminus1),'k')
% hold on;
% plot(1:LTminus1, axon(1).J_l(farnode, 1:LTminus1),'c')
% hold on;
% xlabel('Time');
% ylabel('Current value');
% title('Various currents blue sodium, red potassium, cyan leakage, black p, at position far from stimpos1');
% 

% Ionic current
figure(9)
plot((1:LTminus1).*delT, axon(1).J_Na(stimpos1, 1:LTminus1),'b')
hold on;
plot((1:LTminus1).*delT, axon(1).J_K(stimpos1, 1:LTminus1),'r')
xlabel('Time');
ylabel('Current magnitude');
title('Ionic currents: blue sodium, red potassium at position of stimulation with time');
% 
% subplot(2,3,3)
% mesh(1:(timesteps), 1:endlength, axon(1).V(:,:));
% xlabel('Time');
% ylabel('Position');
% zlabel('Voltage');
% 
figure(10)
farnode = 160;
plot(1:LTminus1, axon(1).m(farnode, 1:LTminus1),'b')
hold on;
plot(1:LTminus1, axon(1).n(farnode,1:LTminus1), 'r')
hold on;
plot(1:LTminus1, axon(1).h(farnode,1:LTminus1), 'k')
hold on;
plot(1:LTminus1, axon(1).p(farnode,1:LTminus1), 'c')
hold on;
xlabel('Time');
ylabel('Gating variable magnitude');
title('Gating variables with time at position 160 which is downaxon from the position of injection');
% 
figure(11)
plot(1:LTminus1, axon(1).m(stimpos1, 1:LTminus1),'b')
hold on;
plot(1:LTminus1, axon(1).n(stimpos1,1:LTminus1), 'r')
hold on;
plot(1:LTminus1, axon(1).h(stimpos1,1:LTminus1), 'k')
hold on;
plot(1:LTminus1, axon(1).p(stimpos1,1:LTminus1), 'c')
hold on;
xlabel('Time');
ylabel('Gating variable magnitude');
title('Gating variables with time at stimulated node');
% 
% subplot(2,3,6)
% plot(1:LTminus1, axon(1).J_ion(farnode, 1:LTminus1),'b');
% xlabel('Time');
% ylabel('Ionic current');
% title('Ionic current at position far from stimpos1');
% 
disp(ratio);
disp(rf);
