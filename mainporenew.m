function mainporenew(mnow,nnow,hnow,pnow)
%% Set a generic width. This corresponds to a normalized duration of command voltage.
b = 10e-3; % 10 milliseconds
% epsilon_not = 8.854e-12% Farads/meter
% er_saltwater = 75; % relative permittivity of salt water: Hasted, Ritson and Collie (1948)
% er_membrane = 9 % at 8 nm thickness, for human red blood cells: Gimsa, Muler, Schelle and Fuhr (1996)
% er_cytoplasm = 50 % for human rbc: Pauly and Schwan (1966)
% The above two values are from Simeonova, Wachner and Gimsa (2002)
permittivity = .01; %epsilon_not*er_saltwater%.01; % actually this is the permittivity - needs to be relabelled throughout.
fieldzero = 1e-40; % any small value
loops = 100; % resolution in time
N = 2; % number of patch clamp experiments to record single channel profile.
%% Generate the random start point
%start = abs((b/2).*randn([N,1]));
%% Generate the random width
%width = abs(0.5*(b/2).*randn([N,1]));
%F(loops) = struct('cdata', [], 'colormap', []);
SummedF(N) = struct('x', [], 'y', [], 'field', []);

initialreturnedvalue = poremovieApril23methodnew(1,b,loops,fieldzero,permittivity);
SummedF(1).x = initialreturnedvalue{1};
SummedF(1).y = initialreturnedvalue{2};
SummedF(1).field = initialreturnedvalue{3};
if N >= 2
for i = 2:N
   
    returnedvalue = poremovieApril23methodnew(start(i), width(i),i,b,loops,fieldzero,permittivity);
    for j = 1:loops
    SummedF(i).x{1,j} = returnedvalue{1}{1,j};
    SummedF(i).y{1,j} = returnedvalue{2}{1,j};
    SummedF(i).field{1,j} = returnedvalue{3}{1,j} + SummedF(i-1).field{1,j};
    end
end
end
% Find the maximum value of the Electric field for the case N=1.

for j = 1:loops
    mymax(j) = max(max(SummedF(1).field{1,j}));
end
[realmax, maxindex] = max(mymax)
for j = 1:loops
    mymin(j) = min(min(SummedF(1).field{1,j}));
end
[realmin, minindex] = min(mymin)
% Compute the potential drop across an axon of radius 5 microns.


% %Display the combined movie.
% vidObj = VideoWriter('Myporefield.avi');
% vidObj.FrameRate = 4;
% vidObj.Quality = 100;
% open(vidObj);

for j = maxindex
    for i = 1%:N
    [c, h] = contourf(SummedF(i).x{1,j}, SummedF(i).y{1,j},SummedF(N).field{1,j}, 20); 
    
      
       colorbar
       %clabel(c,h)
%       caxis('auto');
%         caxis([0 0.24*N/permittivity])
caxis([0 realmax])
       axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
       xlabel('\rho (meters)');

       ylabel('z (meters)');
       title(strcat(num2str(b*j/loops), '(seconds)'));
        drawnow
        
%         currFrame = getframe(gcf);
%         writeVideo(vidObj,currFrame);

        %F(j) = getframe(gcf);
    end

end

% close(vidObj);

%colorbar
%movie(figure, F, 1, 15);

% 
% filename = strcat('newPoresMay0416_', num2str(N), '.avi');
% movie2avi(F, filename, 'compression', 'None', 'FPS', 5);
