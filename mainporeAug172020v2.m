function result = mainporeAug172020v2(mnow,nnow,hnow,pnow,t,oldm,oldn,oldh,oldp)

%% Set a generic width. This corresponds to a normalized duration of command voltage.
%b = 10e-3; % 10 milliseconds
% epsilon_not = 8.854e-12% Farads/meter
% er_saltwater = 75; % relative permittivity of salt water: Hasted, Ritson and Collie (1948)
% er_membrane = 9 % at 8 nm thickness, for human red blood cells: Gimsa, Muler, Schelle and Fuhr (1996)
% er_cytoplasm = 50 % for human rbc: Pauly and Schwan (1966)
% The above two values are from Simeonova, Wachner and Gimsa (2002)
permittivity = .01; %epsilon_not*er_saltwater%.01; % actually this is the permittivity - needs to be relabelled throughout.
fieldzero = 1e-40; % any small value
%loops = 100; % resolution in time
N = 1; % number of patch clamp experiments to record single channel profile.
%% Generate the random start point
%start = abs((b/2).*randn([N,1]));
%% Generate the random width
%width = abs(0.5*(b/2).*randn([N,1]));
%F(loops) = struct('cdata', [], 'colormap', []);
SummedF(N) = struct('x', [], 'y', [], 'field', []);

myinitialreturnedvalue = poremovieAug17methodnewv2(1,fieldzero,permittivity, mnow,nnow,hnow,pnow, oldm, oldn, oldh, oldp);
mylength=31;
mylength2 = 31;
SummedF(1).x = myinitialreturnedvalue{1,1};
SummedF(1).y = myinitialreturnedvalue{1,2};
SummedF(1).field = myinitialreturnedvalue{1,3};
SummedF(1).oldsourcelowv = myinitialreturnedvalue{1,4};
if N >= 2
for i = 2:N
   
    myreturnedvalue = poremovieAug17methodnewv2(i,fieldzero,permittivity,mnow,nnow,hnow,pnow,oldm, oldn, oldh, oldp);
    j=1;
%     SummedF(i).x = returnedvalue(1,1);
%     SummedF(i).y = returnedvalue(1,2);
    SummedF(i).field{:,:} = myreturnedvalue{1,3} + SummedF(i-1).field{:,:};
    
end
end
% Find the maximum value of the Electric field for the case N=1.

for j = 1
    mymax(j) = max(max(SummedF(1).field));
end
[realmax, maxindex] = max(mymax);
for j = 1
    mymin(j) = min(min(SummedF(1).field));
end
[realmin, minindex] = min(mymin);
% Compute the potential drop across an axon of radius 5 microns.


% %Display the combined movie.
% vidObj = VideoWriter('Myporefield.avi');
% vidObj.FrameRate = 4;
% vidObj.Quality = 100;
% open(vidObj);


result{1,1} = SummedF(1).x;
result{1,2} = SummedF(1).y;
result{1,3} = SummedF(1).field;
result{1,4} = realmin;


% for j = maxindex
%     [c, h] = contourf(real(SummedF(1).x), real(SummedF(1).y),real(SummedF(1).field), 20); 
%     
%       
%        colorbar
%        %clabel(c,h)
% %       caxis('auto');
%          caxis([200 10000])
% %caxis([0 realmax])
%        axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
%        xlabel('\rho (meters)');
% 
%        ylabel('z (meters)');
%        title(num2str(t));
%         drawnow
        
%         currFrame = getframe(gcf);
%         writeVideo(vidObj,currFrame);

        %F(j) = getframe(gcf);
%     end

end

% close(vidObj);

%colorbar
%movie(figure, F, 1, 15);

% 
% filename = strcat('newPoresMay0416_', num2str(N), '.avi');
% movie2avi(F, filename, 'compression', 'None', 'FPS', 5);
