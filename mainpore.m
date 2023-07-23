close all;
clear all;
clc
clf
%% Set a generic width. This corresponds to a normalized duration of command voltage.
b = 10e-3; % 10 milliseconds
conductivity = .01;
fieldzero = 1e-10;
loops = 100; % resolution in time
N = 1; % number of patch clamp experiments to record single channel profile.
%% Generate the random start point
start = abs((b/2).*randn([N,1]));
%% Generate the random width
width = abs(0.5*(b/2).*randn([N,1]));
%F(loops) = struct('cdata', [], 'colormap', []);
SummedF(N) = struct('x', [], 'y', [], 'field', []);

initialreturnedvalue = poremovieApril23method(start(1),width(1),1,b,loops,fieldzero,conductivity);
SummedF(1).x = initialreturnedvalue{1};
SummedF(1).y = initialreturnedvalue{2};
SummedF(1).field = initialreturnedvalue{3};
if N >= 2
for i = 2:N
   
    returnedvalue = poremovieApril23method(start(i), width(i),i,b,loops,fieldzero,conductivity);
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
realmax = max(mymax)

%Display the combined movie.

for j = 1:loops
    for i = 1:N
    [c, h] = contourf(SummedF(i).x{1,j}, SummedF(i).y{1,j},SummedF(N).field{1,j}, 20); 
    
      
       colorbar
       %clabel(c,h)
%       caxis('auto');
%         caxis([0 0.24*N/conductivity])
caxis([0 realmax])
       axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
       xlabel('\rho');

       ylabel('z');
       title(num2str(b*j/loops));
        drawnow

        F(j) = getframe(gcf);
    end

end
colorbar
movie(figure, F, 1, 15);

% 
filename = strcat('newPoresMay0416_', num2str(N), '.avi');
movie2avi(F, filename, 'compression', 'None', 'FPS', 5);
