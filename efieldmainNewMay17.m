% Find the peak electric field due to an axon at any field position
close all; 
clear all; 
clf; 
clc;
numnodes = 10;
choice = 1; % choice should be set at 0 for the initial method and 1 for the new method.
% The new method gives the right directions (normal to the axon) and also
% the right ball park magnitude (1e8).

numaxons = 1; 

if numaxons == 1

% Plot
numpoints = 10;
% x = linspace(-50e-6*1e2, 50e-6*1e2, numpoints);
% y = linspace(-50e-6*1e2, 50e-6*1e2, numpoints);
% z = linspace(-50e-6*1e2, 50e-6*1e2, numpoints);
% [x,y,z] = meshgrid(x,y,z);

r = linspace(1e-6*1e2,5e-6*1e2, 10);
phi = linspace(0, 2*pi, 10); % Since it should be cylindrically symmetric
%theta = linspace(0,pi, 9);help 
theta = linspace(30*2*pi/360,pi-30*2*pi/360, 10);
coordstore = cell(length(r), length(theta), length(phi));
fieldstore = cell(length(r), length(theta), length(phi));

for i = 1:length(r)
    for j = 1:length(theta)
       for k= 1:length(phi)
        %[r, theta, phi] = cart2spher(x(i), y(j), z(k));
        [Efieldx, Efieldy, Efieldz] = axonefield(numnodes, r(i), theta(j), phi(k), choice, 0, 0);
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
        s = sqrt(U(i,j,k)^2 + V(i,j,k)^2 + W(i,j,k)^2);
        sU(i,j,k) = U(i,j,k)/s;
        sV(i,j,k) = V(i,j,k)/s;
        sW(i,j,k) = W(i,j,k)/s;
        end
        
    end
    
end
figure(1)
quiver3(X,Y,Z,sU,sV,sW)
    
figure(2)
Emag = sqrt(U.^2 + V.^2 + W.^2);
L = Emag;
 cn=8; % <-  number of conditions
     cmap=jet(cn); % <- select colors.
     tftmpl=false(size(L)); % <- create an all-false template.
for i=1:cn
     tf=tftmpl; % <- start with a clean tf.
if i==1
    tf(L<1e2) = true;
elseif i == 2
    tf(L>=1e2&L<1e3) = true;
elseif i == 3
    tf(L>=1e3 & L<1e4) = true;
elseif i == 4
    tf(L>=1e4 & L<1e5) = true;
elseif i == 5
    tf(L>=1e5 & L<1e6) = true;
elseif i==6
    tf(L>=1e6 & L < 1e7) = true;
elseif i == 7
    tf(L>=1e7 & L < 5e7) = true;
else
     tf(L>=5e7)=true;
end

     coneplot(X(tf),Y(tf),Z(tf),sU(tf),sV(tf),sW(tf));%,'color',cmap(i,:));
     hold on;
end
Az = 30;
EL = 30;
view(Az, EL)

% Find the peak field magnitude 

Emag = sqrt(U.^2 + V.^2 + W.^2);
max(max(max(Emag)))

elseif numaxons == 2
    
     % First work on the first axon, located at the origin.
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
        
        [Efieldx, Efieldy, Efieldz] = axonefield(numnodes, r(i), theta(j), phi(k), choice, 0, 0);
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

% Now work on the second axon located at
axon2rmag = 200e-6;
axon2theta = pi/2;
axon2phi = pi/32;

    

     % Plot
axon2numpoints = 10;
r = linspace(1e-6*1e2,500e-6*1e2, axon2numpoints);
phi = linspace(0, 2*pi,axon2numpoints); % Since it should be cylindrically symmetric
theta = linspace(7*2*pi/360,pi-7*2*pi/360, axon2numpoints);
coordstore2 = cell(length(r), length(theta), length(phi));
fieldstore2 = cell(length(r), length(theta), length(phi));

for i = 1:length(r)
    for j = 1:length(theta)
       for k= 1:length(phi)
        
        [Efieldx2, Efieldy2, Efieldz2] = axonefield(numnodes, r(i), theta(j), phi(k), choice, axon2rmag, axon2theta);
        [mycartcoordsx2, mycartcoordsy2, mycartcoordsz2] = spher2cart(r(i), theta(j), phi(k)); % field positions in Cartesian coordinates
        coordstore2{i,j,k} = [mycartcoordsx2, mycartcoordsy2, mycartcoordsz2];
        fieldstore2{i,j,k} = [Efieldx2, Efieldy2, Efieldz2];
       end
    end
end

% Develop X, Y, Z, U, V, W

for i = 1:length(r)
    for j = 1:length(theta)
        for k = 1:length(phi)
        X2(i,j,k) = 1e-2.*coordstore2{i,j,k}(1);
        Y2(i,j,k) = 1e-2.*coordstore2{i,j,k}(2);
        Z2(i,j,k) = 1e-2.*coordstore2{i,j,k}(3);
        U2(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore2{i,j,k}(1); % convert to SI volts per meter
        V2(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore2{i,j,k}(2);
        W2(i,j,k) = (1/(1e6*(1/3e10))).*fieldstore2{i,j,k}(3);
        end
        
    end
    
end
figure(3)
Emag2 = sqrt(U2.^2 + V2.^2 + W2.^2);
L2 = Emag2;
 cn=4; % <-  number of conditions
     cmap=jet(cn); % <- select colors.
     tftmpl=false(size(L2)); % <- create an all-false template.
for i=1:cn
     tf=tftmpl; % <- start with a clean tf.
if i==1
         tf(L2<1e6) = true;

elseif i==2
    tf(L2>=1e6 & L2 < 1e7) = true;
elseif i == 3
    tf(L2>=1e7 & L2 < 5e7) = true;
else
     tf(L2>=5e7)=true;
end

     quiver3(X2(tf),Y2(tf),Z2(tf),U(tf) + U2(tf),V(tf) + V2(tf),W(tf) + W2(tf), 'color',cmap(i,:));
     
     
     hold on;
end

figure(4)
 quiver3(X2,Y2,Z2,U + U2,V + V2,W + W2);
% hold on
% 
% Emag = sqrt(U.^2 + V.^2 + W.^2);
% L = Emag;
%  cn=4; % <-  number of conditions
%      cmap=jet(cn); % <- select colors.
%      tftmpl=false(size(L)); % <- create an all-false template.
% for i=1:cn
%      tf=tftmpl; % <- start with a clean tf.
% if i==1
%          tf(L<1e6) = true;
% 
% elseif i==2
%     tf(L>=1e6 & L < 1e7) = true;
% elseif i == 3
%     tf(L>=1e7 & L < 5e7) = true;
% else
%      tf(L>=5e7)=true;
% end
% 
%      quiver3(X(tf),Y(tf),Z(tf),U(tf),V(tf),W(tf), 'color',cmap(i,:));
%      hold on;


Az = 30;
EL = 30;
view(Az, EL)

% % Find the peak field magnitude 
% 
% Emag = sqrt(U.^2 + V.^2 + W.^2);
% max(max(max(Emag)))
    else
        % do nothing
    end
    
