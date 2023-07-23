% Find the peak electric field due to an axon at any field position
close all; 
clear all; 
clf; 
clc;
numnodes = 5;
choice = 1; % choice should be set at 0 for the old method (discussed with Dr. Snider) and 
% 1 for the new method. The new method gives the right directions (normal to the axon) and also
% the right ball park magnitude (1e8). 





% % Plot
% 
% r = linspace(1.5e-6*1e2,1000e-6*1e2, 10);
% phi = pi/2; % Since it should be cylindrically symmetric
% theta = linspace(pi/(2*90),pi/2, 10);
% coordstore = cell(length(r), length(theta));
% fieldstore = cell(length(r), length(theta));
% 
% for i = 1:length(r)
%     for j = 1:length(theta)
%         disp(i)
%         disp(j)
%         phi = pi/2;
%         [Efieldx, Efieldy, Efieldz] = axonefield(numnodes, r(i), theta(j), phi);
%         [mycartcoordsx, mycartcoordsy, mycartcoordsz] = spher2cart(r(i), theta(j), phi); % field positions in Cartesian coordinates
%         coordstore{i,j} = [mycartcoordsx, mycartcoordsy, mycartcoordsz];
%         fieldstore{i,j} = [Efieldx, Efieldy, Efieldz];
%     end
% end
% 
% % Develop X, Y, Z, U, V, W
% 
% for i = 1:length(r)
%     for j = 1:length(theta)
%         
%         X(i,j) = 1e-2.*coordstore{i,j}(1);
%         Y(i,j) = 1e-2.*coordstore{i,j}(2);
%         Z(i,j) = 1e-2.*coordstore{i,j}(3);
%         U(i,j) = (1/(1e6*(1/3e10))).*fieldstore{i,j}(1);
%         V(i,j) = (1/(1e6*(1/3e10))).*fieldstore{i,j}(2);
%         W(i,j) = (1/(1e6*(1/3e10))).*fieldstore{i,j}(3);
%     end
% end
% figure(1)
% quiver3(X, Y, Z, U, V, W);

        

 % Plot
numpoints = 10;
r = linspace(1e-6*1e2,500e-6*1e2, numpoints);
phi = linspace(0, 2*pi,numpoints); % Since it should be cylindrically symmetric
theta = linspace(7*2*pi/360,pi-7*2*pi/360, numpoints);
coordstore = cell(length(r), length(theta), length(phi));
fieldstore = cell(length(r), length(theta), length(phi));

for i = 1:length(r)
    for j = 1:length(theta)
       for k= 1:length(phi)
        
        [Efieldx, Efieldy, Efieldz] = axonefield(numnodes, r(i), theta(j), phi(k), choice);
        [mycartcoordsx, mycartcoordsy, mycartcoordsz] = spher2cart(r(i), theta(j), phi(k)); % field positions in Cartesian coordinates
        coordstore{i,j,k} = [mycartcoordsx, mycartcoordsy, mycartcoordsz];
        fieldstore{i,j,k} = [Efieldx, Efieldy, Efieldz];
       end
    end
end

% Develop X, Y, Z, U, V, W

for i = 1:length(r)
    for j = 1:length(theta)
        for k = 1:length(phi)
        X(i,j,k) = 1e-2.*coordstore{i,j,k}(1);
        Y(i,j,k) = 1e-2.*coordstore{i,j,k}(2);
        Z(i,j,k) = 1e-2.*coordstore{i,j,k}(3);
        U(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore{i,j,k}(1); % convert to SI volts per meter
        V(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore{i,j,k}(2);
        W(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore{i,j,k}(3);
        end
        
    end
    
end

    
figure(2)
Emag = sqrt(U.^2 + V.^2 + W.^2);
L = Emag;
 cn=4; % <- #conditions
     cmap=jet(cn); % <- select your own colors...
     tftmpl=false(size(L)); % <- create an all-false template...
for i=1:cn
     tf=tftmpl; % <- start with a clean tf...
if i==1
         tf(L<1e6) = true;

elseif i==2
    tf(L>=1e6 & L < 1e7) = true;
elseif i == 3
    tf(L>=1e7 & L < 5e7) = true;
else
     tf(L>=5e7)=true;
end
% the data
% - modified from HELP QUIVER3
%      [x,y]=meshgrid(-2:.5:2,-2:.5:2);
%      z=x .*exp(-x.^2-y.^2);
%      [u,v,w]=surfnorm(x,y,z);
% the engine
    % tf=U<0; % <- color is different for U<0 and U>=0
% the result
     quiver3(X(tf),Y(tf),Z(tf),U(tf),V(tf),W(tf), 'color',cmap(i,:));
     hold on;
end
%colorbar
%      hold on;
%      tf=~tf;
%      quiver3(X(tf),Y(tf),Z(tf),U(tf),V(tf),W(tf),'color',[0,1,0],'linewidth',1);
%      %surf(X,Y,Z);
%      colormap(bone);


% quiver3(X, Y, Z, U, V, W);
Az = 30;
EL = 30;
view(Az, EL)

% Find the peak field magnitude 

Emag = sqrt(U.^2 + V.^2 + W.^2);
max(max(max(Emag)))
% C = repmat([1 2 3],numel(X),1);
% figure(4)
% scatter3(X(:,1), Y(:,1), Z(:,1), Emag, 'r')


% figure(3)
% for j = 1:length(phi)
% surf(X(:, :, j), Y(:, :, j), Z(:, :, j), Emag(:, :, j))
% hold on;
% end

 % Plot
%  internodallength = 50e-6*100;
% 
% numpoints = 10;
% x = [linspace(-5e-6*1e2,-1.5e-6*1e2, numpoints/2), linspace(1.5e-6*1e2, 5e-6*1e2, numpoints/2)];
% y = x;
% z = linspace(-numnodes*internodallength, numnodes*internodallength, numpoints);
% for i = 1:length(x)
%     [r(i), phi(i), theta(i)] = cart2spher(x(i), y(i), z(i));
% end
% 
% coordstore = cell(length(r), length(theta), length(phi));
% fieldstore = cell(length(r), length(theta), length(phi));
% 
% for i = 1:length(r)
%     for j = 1:length(theta)
%        for k= 1:length(phi)
%         
%         [Efieldx, Efieldy, Efieldz] = axonefield(numnodes, r(i), theta(j), phi(k));
%         [mycartcoordsx, mycartcoordsy, mycartcoordsz] = spher2cart(r(i), theta(j), phi(k)); % field positions in Cartesian coordinates
%         coordstore{i,j,k} = [mycartcoordsx, mycartcoordsy, mycartcoordsz];
%         fieldstore{i,j,k} = [Efieldx, Efieldy, Efieldz];
%        end
%     end
% end
% 
% % Develop X, Y, Z, U, V, W
% 
% for i = 1:length(r)
%     for j = 1:length(theta)
%         for k = 1:length(phi)
%         X(i,j,k) = 1e-2.*coordstore{i,j,k}(1);
%         Y(i,j,k) = 1e-2.*coordstore{i,j,k}(2);
%         Z(i,j,k) = 1e-2.*coordstore{i,j,k}(3);
%         U(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore{i,j,k}(1);
%         V(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore{i,j,k}(2);
%         W(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore{i,j,k}(3);
%         end
%         
%     end
%     
% end
%     
% figure(3)
% 
% quiver3(X, Y, Z, U, V, W);
% 
% Emag = sqrt(U.^2 + V.^2 + W.^2);
% 
% 
% 
% 
