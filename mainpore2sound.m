clear all;
clc;
close all;
%function [y1, y2, y3, start, width] = mainpore2sound(option)
% if option is 1, we are interested in sound production for a single patch
% clamp experiment. if option is 2, we are interested in generating the
% movie.
% 
% if option==2
% 
% %% Set a generic width. This corresponds to a normalized duration of command voltage.
% b = 10e-3; % 10 milliseconds
% conductivity = .01;
% fieldzero = 1e-10;
% loops = 100; % resolution in time
% N = 4; % number of patch clamp experiments to record single channel profile.
% %% Generate the random start point
% start = abs((b/2).*randn([N,1]));
% %% Generate the random width
% width = abs(0.5*(b/2).*randn([N,1]));
% %F(loops) = struct('cdata', [], 'colormap', []);
% SummedF(N) = struct('x', [], 'y', [], 'field', []);
% 
% initialreturnedvalue = poremovieApril23method(start(1),width(1),1,b,loops,fieldzero,conductivity);
% SummedF(1).x = initialreturnedvalue{1};
% SummedF(1).y = initialreturnedvalue{2};
% SummedF(1).field = initialreturnedvalue{3};
% if N >= 2
% for i = 2:N
%    
%     returnedvalue = poremovieApril23method(start(i), width(i),i,b,loops,fieldzero,conductivity);
%     for j = 1:loops
%     SummedF(i).x{1,j} = returnedvalue{1}{1,j};
%     SummedF(i).y{1,j} = returnedvalue{2}{1,j};
%     SummedF(i).field{1,j} = returnedvalue{3}{1,j} + SummedF(i-1).field{1,j};
%     end
% end
% end
% 
% 
% % Find the maximum value of the Electric field for the case N=1.
% 
% for j = 1:loops
%     mymax(j) = max(max(SummedF(1).field{1,j}));
% end
% realmax = max(mymax)
% 
% %Display the combined movie.
% 
% for j = 1:loops
%     for i = 1:N
%     [c, h] = contourf(SummedF(i).x{1,j}, SummedF(i).y{1,j},SummedF(N).field{1,j}, 20); 
%     
%       
%        colorbar
%        %clabel(c,h)
% %       caxis('auto');
% %         caxis([0 0.24*N/conductivity])
% caxis([0 realmax])
%        axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
%        xlabel('\rho');
% 
%        ylabel('z');
%        title(num2str(b*j/loops));
%         drawnow
% 
%         F(j) = getframe(gcf);
%     end
% 
% end
% colorbar
% movie(figure, F, 1, 15);
% 
% % 
% filename = strcat('newPoresMay0416_', num2str(N), '.avi');
% movie2avi(F, filename, 'compression', 'None', 'FPS', 5);
% x = 0;
% y= 0;
% field = 0;
%if option==1
    %% Set a generic width. This corresponds to a normalized duration of command voltage.
factor = 2*10000000;
b = 2*10e-2; % 10 milliseconds
conductivity = .01;
fieldzero = 0;
loops = 88200; % resolution in time
N = 2; % number of patch clamp experiments to record single channel profile.
%% Generate the random start point
start = abs((b/2).*randn([N,1]));
%% Generate the random width
width = abs(0.5*(b/2).*randn([N,1]));
%F(loops) = struct('cdata', [], 'colormap', []);
SummedF(N) = struct('x', [], 'y', [], 'field', []);
% i=1;
% initialreturnedvalue = poremovieApril23method(start(1),width(1),1,b,loops,fieldzero,conductivity);
% SummedF(i).x = initialreturnedvalue{1};
% SummedF(i).y = initialreturnedvalue{2};
% SummedF(i).field = initialreturnedvalue{3};
% 
% 
% if N >= 2
% for i = 2:N
%    
%     returnedvalue = poremovieApril23method(start(i), width(i),i,b,loops,fieldzero,conductivity);
%     for j = 1:loops
%     SummedF(i).x{1,j} = returnedvalue{1}{1,j};
%     SummedF(i).y{1,j} = returnedvalue{2}{1,j};
%     SummedF(i).field{1,j} = returnedvalue{3}{1,j} + SummedF(i-1).field{1,j};
%     end
% end
% end
% 
% x = SummedF(i).x;
% y = SummedF(i).y;
% field = SummedF(i).field;




%% Sound of a polynomial

%This Program generates a Frequency Modulated(Chirp) signal using a
%polynomial as message signal

%% Screen Display Specifications

%Measure screen size of the device
%Calculate position values of figure Windows

scrsz = get(0,'ScreenSize');
P1=[40 500 scrsz(3)/3 scrsz(4)/3];
P2=[40 80 scrsz(3)/3 scrsz(4)/3];
P3=[600 500 scrsz(3)/3 scrsz(4)/3];
P4=[600 80 scrsz(3)/3 scrsz(4)/3];
P5=[1000 500 scrsz(3)/3 scrsz(4)/3];
P6=[1000 80 scrsz(3)/3 scrsz(4)/3];

%% Parameter Sepecifications

Fs  = 44100;              % Set the Sampling frequency
T   = 1/Fs;               % Calculate Sample time
f   = 200;                % Set Carrier frequency of 200Hz
s   = 2;                  % Set the duration of signal to be generated
L   = s*Fs;               % Calculate the length of signal to be generated
t   = (0:L-1)*T;          % Time vector
a   = 1;                  % Set a, b, and c coefficients here
b   = 3;
c   = 2;

%% Signal Generation Block

% xstepsize = 0.1e-6;
% ystepsize = 0.1e-6;
%  x=2.5e-6:xstepsize:5.5e-6; % this is the cylindrical rho coordinate range, from the surface 
%  % of the node of ranvier at 2.5 microns to about 5.5 microns away from the
%  % surface.
%  y=0.1e-6:ystepsize:3.1e-6; % this is the cylindrical z coordinate range, which goes 
%  % from 0.1 micro meters to about 3 microns along the length of the axon.
%         
%  mylength = length(x);
% for i = 1:mylength
%     for j = 1:mylength
%         mygrid{i,j}  = {x(i),y(j)};
%     end
% end
% % Extract the x-coordinates
% for i = 1:mylength
%     for j = 1:mylength
%         px(i,j) = mygrid{i,j}(1);
%     end
% end
% % Extract the y-coordinates
% for i = 1:mylength
%     for j = 1:mylength
%         py(i,j) = mygrid{i,j}(2);
%     end
% end
% pxnum = cell2mat(px);
% pynum = cell2mat(py);
% % conductivity = 1; %conductivity of extracellular medium Siemens per meter.
% 
% rho1prime = 1e-6; % lowu, lowu2, in microns - where the current source is at its peak strength.
% 
% Inot = 2e-12/(2*pi*rho1prime*0.1e-6); % peak current strength at surface of node is 2 pico Amperes. Thickness of 
% % nodal cylindrical patch that acts as the source is taken to be 0.1
% % microns roughly. radius is rho1prime. Thus Inot is the peak current per
% % unit area emanating from the source cylindrical annulus at the node.
% rho2prime = 2e-6; % highu, highu2, in microns - where the current source zeros out.
% %For each (x,y) pair in the above matrix, compute the two integrals and
% %store them as the px and py, respectively rho and z integrals. 
% for i = 1:mylength
%     for j = 1:mylength
%         [Q, Q2] = zandrhointegrals(pxnum(i,j),pi/2,pynum(i,j));
%         p{i,j} = {Q2*((Inot/(rho2prime-rho1prime))/(4*pi*conductivity)),Q*(pynum(i,j)*(Inot/(rho2prime-rho1prime))/(4*pi*conductivity))};
%     end
% end
% 
% % Extract the x-coordinates of the E-field vector, in units of Volt per
% % meter (since all units are in S.I.)
% for i = 1:mylength
%     for j = 1:mylength
%         vecx(i,j) = p{i,j}(1);
%     end
% end
% % Extract the y-coordinates of the vector
% for i = 1:mylength
%     for j = 1:mylength
%         vecy(i,j) = p{i,j}(2);
%     end
% end
% vecxnum = cell2mat(vecx);
% vecynum = cell2mat(vecy);

% we will choose the field along the first row and first column of vecxnum and vecynum.
for j = 1:loops
       y1(j) = randomtrapezoid(fieldzero,start(1), width(1),1,b,b*j/loops);
end

for i = 1:N
    y1old = y1;
for j = 1:loops
       y1(j) = randomtrapezoid(fieldzero,start(i), width(i),1,b,b*j/loops);
end
y1 = y1 + y1old;
end

z1 = y1.*sin(2*pi*f*t);
z1 = z1';% Generate sinusoid 
y2 = 1.0*sin(2*pi*f*t.^2)';                 % Generate modulated sinusoid 
y3 = 1.0*sin(2*pi*f*a*(t.^3)+b*(t.^2)+c)';  % Generate modulated sinusoid 

%% Visualisation of the Generated Signals

% % Plotting Carrier Signal
% 
% figure('position', P1);
% figure(1);
% plot(t,y1);
% grid on
% title('Carrier Signal');
% xlabel('Time in seconds');
% ylabel(['Carrier' num2str(f) 'Hz']);

%% Plotting Linear FM Signal

figure('position', P2);
figure(2);
plot(t,y2);
grid on
title('Synthetic FM Signal1');
xlabel('Time in seconds');
ylabel(['Linear FM Signal']);

%% Plotting Non Linear FM Signal

figure('position', P3);
figure(3);
plot(t,y3);
grid on
title('Synthetic FM Signal1');
xlabel('Time in seconds');
ylabel('Non Linear FM Signal');

% %% Spectrogram Display 
% 
% figure('position',P4);
% figure(4);
% spectrogram(y1,256);
% text(0.4,28,'Frame Size 256','FontSize',12);
% 
% figure('position',P5);
% figure(5);
% spectrogram(y2,256);
% text(0.4,28,'Frame Size 256','FontSize',12);
% 
% figure('position',P6);
% figure(6);
% spectrogram(y3,256);
% text(0.4,28,'Frame Size 256','FontSize',12);

%% Playback of the Synthesised Signals
 
%soundsc(y1,Fs)
%soundsc(y2,Fs)
%soundsc(y3,Fs)

%% Save the Resulting Files in Wav Format
audiowrite('Carrier.wav', z1, Fs);
audiowrite('LinearChirp.wav', y2, Fs);
audiowrite('NonLinearChirp.wav', y3, Fs);

%% End of Program
%end




