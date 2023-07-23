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

function returncellarray = poremovieAug17methodnew(trialnum,fieldzero, permittivity, mnow, nnow,hnow,pnow, oldans, oldm, oldn, oldh, oldp)



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

 x=5.5e-6;%2.5e-6:xstepsize:5.5e-6; % this is the cylindrical rho coordinate range, from the surface 
 % of the node of ranvier at 2.5 microns to about 5.5 microns away from the
 % surface.
 y=3.1e-6;%0.1e-6:ystepsize:3.1e-6; % this is the cylindrical z coordinate range, which goes 
 % from 0.1 micro meters to about 3 microns along the length of the axon.
        
 mylength = length(x);
 mylength2 = length(y);
for i = 1:mylength
    for j = 1:mylength2
        grid{i,j}  = {x(i),y(j)};
    end
end
% Extract the x-coordinates
for i = 1:mylength
    for j = 1:mylength2
        px(i,j) = grid{i,j}(1);
    end
end
% Extract the y-coordinates
for i = 1:mylength
    for j = 1:mylength2
        py(i,j) = grid{i,j}(2);
    end
end
pxnum = cell2mat(px);
pynum = cell2mat(py);
% permittivity = 1; %permittivity of extracellular medium Siemens per meter.
%lnor = 2e-6; % length of node of Ranvier
rho1prime = 1e-6; % lowu, lowu2, in microns - where the current source is at its peak strength.
totnumchan = 10; %number of ion channels in a node, total, sodium
windowsize = ceil(.28*totnumchan);
% if (mnow > 1)||(mnow == NaN)
%     mnow = 1;
% elseif (mnow < 0)||(mnow == NaN)
%     mnow = 0;
% end
% if (nnow > 1)||(nnow == NaN)
%     nnow = 1;
% elseif (nnow < 0)||(nnow == NaN)
%     nnow = 0;
% end
% if (hnow > 1)||(hnow == NaN)
%     hnow = 1;
% elseif (hnow < 0)||(hnow == NaN)
%     hnow = 0;
% end
% hnow
% mnow
% totnumchan
numchannelsactive = ceil(hnow*mnow^3*totnumchan);
if numchannelsactive <1
    numchannelsactive = 1;
elseif numchannelsactive > totnumchan
    numchannelsactive = totnumchan;
else
    %do nothing
%     numchannelsactive = 5;
%     disp('error smudged');
end

% if oldh<0
%     oldh = 0;
% end
% if oldh>1
%     oldh=1;
% end
% if oldm>1
%     oldm = 1;
% end
% if oldm < 0
%     oldm = 0;
% end


numchannelprevactive = ceil(oldh*oldm^3*totnumchan);
if numchannelprevactive <1
    numchannelprevactive = 1;
elseif numchannelprevactive > totnumchan
    numchannelprevactive = totnumchan;
else
     % do nothing
%     numchannelprevactive = 3;
%     disp('error 2 smudged');
end

% numchannelprevactive = ceil(oldh*oldm^3*totnumchan);
if numchannelsactive <1
    numchannelsactive=1;
end
if numchannelsactive == numchannelprevactive
    temp = linspace(0,2*pi,totnumchan);
    b = temp;
    %out of these positions only numchannelsactive positions will be active
    count = 1;
    while count < length(b)
        if count > length(temp)
            break
        else
        prevtemp = temp;
        pr = numchannelsactive/totnumchan;
        if rand(1)>pr%coin toss gives heads, turn on that channel. probability of heads shold
            % depend on numchannelsactive. In fact, probability of heads =
            % numchannelsactive/totnumchan.
            if count>=2
%                 count
%                 size(temp)
%                 size(prevtemp)
            currtemp = prevtemp(1, [1:count-1, count+1:length(temp)]);
            temp = currtemp;
            else
                currtemp = prevtemp(1, count+1:length(temp));
                temp = currtemp;
            end
        else
            %do nothing
        end
  %      temp = currtemp;
        count = count + 1;
        end
    %    temp = currtemp;
    end
            
sourcelowv = temp;% single ion channel position in ring
sourcehighv = sourcelowv+0.01;% phi-width of the channel is added to get the upper limit.
Inot = 2e-12/(pi*(((1e-6)^2)*3e-6)); % peak current strength at surface of node is 2 pico Amperes.
% The denominator contains a computation of the source cylindrical volume.
rho2prime = 2e-6; % highu, highu2, in microns - where the current source zeros out.
%For each (x,y) pair in the above matrix, compute the two integrals and
%store them as the px and py, respectively rho and z integrals. 
loops = 1;
returncellarray = cell(1,3);
field = cell(1,loops);
for j = 1:loops
    field{1,j} = fieldzero;
end


for channelnum = 0:windowsize:length(sourcelowv)-1
%     disp(channelnum)
for i = 1:mylength
    for j = 1:mylength
        [Q, Q2] = zandrhointegralsnew(sourcelowv(channelnum+1), sourcehighv(channelnum+1),pxnum(i,j),pi/2,pynum(i,j));
        p{i,j} = {Q2*((Inot/(rho2prime-rho1prime))/(4*pi*permittivity)),Q*(pynum(i,j)*(Inot/(rho2prime-rho1prime))/(4*pi*permittivity))};
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
vecxnum = windowsize.*cell2mat(vecx);
vecynum = windowsize.*cell2mat(vecy);

%axis tight manual;
%ax.NextPlot = 'replaceChildren';
%deltat = b/loops; % size of time steps

loops = 1;
% F(loops) = struct('cdata', [], 'colormap', []);
x = cell(1,loops);
y = cell(1,loops);
%field = cell(1,loops);
for j = 1:loops
        x{1,j} = pxnum;
        y{1,j} = pynum;
        field{1,j} = field{1,j} + real(sqrt(vecxnum.^2 + vecynum.^2));
end
% disp('hi1 from poremovie')
% size(pxnum)
% size(pynum)
%returncellarray = cell(1,3);
returncellarray{1,1} = real(x{1,1});
returncellarray{1,2} = real(y{1,1});
end
returncellarray{1,3} = real(field{1,1});
returncellarray{1,4} = sourcelowv;




elseif (numchannelsactive-numchannelprevactive)<0 
    % the idea is to save computational effort
    % so we take the previously computed field that is, the previously
    % active channels as still active and only a few additional ones are to
    % be made active probabilistically (this can be modified in future to
    % take into account that some of the prviously active will also close
    % down)
    % first find how many additional channels need to be made active
    addchannel = ceil(numchannelsactive-numchannelprevactive);
    % if this is negative, then we actually need to close addchannel number
    % of channels
  
        % close abs(addchannel) number of channels from among the open
        % ones. 
        
        temp = oldans{1,4};
         prevtemp = zeros(1,1);
    b = temp;
    %out of these positions close abs(addchannel) of them.
    count = 1;
    while count < length(b)
        if count > length(temp)
            break
        else
       
        pr2 = abs(addchannel)/totnumchan;
        if rand(1)>pr2
            % since this channel is to be closed, simply note its position and 
            % subtract the field due to it.
           
            currtemp = union(prevtemp, temp(count));
           
        else
            %do nothing
        end
         prevtemp = temp;
        temp = currtemp;
       
        count = count + 1;
        end
    end
sourcelowv = temp;% single ion channel position in ring
sourcehighv = sourcelowv+0.01;% phi-width of the channel is added to get the upper limit.
Inot = 2e-12/(pi*(((1e-6)^2)*3e-6)); % peak current strength at surface of node is 2 pico Amperes.
% The denominator contains a computation of the source cylindrical volume.
rho2prime = 2e-6; % highu, highu2, in microns - where the current source zeros out.
%For each (x,y) pair in the above matrix, compute the two integrals and
%store them as the px and py, respectively rho and z integrals. 
loops = 1;
returncellarray = cell(1,3);
field = cell(1,loops);
for j = 1:loops
    field{1,j} = fieldzero;
end


for channelnum = 0:windowsize:length(sourcelowv)-1
    disp(channelnum)
for i = 1:mylength
    for j = 1:mylength
        [Q, Q2] = zandrhointegralsnew(sourcelowv(channelnum+1), sourcehighv(channelnum+1),pxnum(i,j),pi/2,pynum(i,j));
        p{i,j} = {Q2*((Inot/(rho2prime-rho1prime))/(4*pi*permittivity)),Q*(pynum(i,j)*(Inot/(rho2prime-rho1prime))/(4*pi*permittivity))};
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
vecxnum = windowsize.*cell2mat(vecx);
vecynum = windowsize.*cell2mat(vecy);

%axis tight manual;
%ax.NextPlot = 'replaceChildren';
%deltat = b/loops; % size of time steps

loops = 1;
% F(loops) = struct('cdata', [], 'colormap', []);
x = cell(1,loops);
y = cell(1,loops);
%field = cell(1,loops);
for j = 1:loops
        x{1,j} = pxnum;
        y{1,j} = pynum;
        field{1,j} = field{1,j} + real(sqrt(vecxnum.^2 + vecynum.^2));
end
% disp('hi2 from poremovie')
% size(pxnum)
% size(pynum)
%returncellarray = cell(1,3);
returncellarray{1,1} = real(x{1,1});
returncellarray{1,2} = real(y{1,1});
end
returncellarray{1,3} = oldans{1,3}(1,1) -real(field{1,1});
returncellarray{1,4} = sourcelowv;




else
%     numchannelsactive
%     numchannelprevactive
    addchannel = ceil(numchannelsactive-numchannelprevactive);
    % if this is positive, then we need to include addchannel number
    % of channels
  
        % open addchannel number of channels from among the open
        % ones. 
        
        temp = linspace(0,2*pi,totnumchan); % this was the original grid on the node.
        b = temp;
        temp2 = oldans{1,4}; % these are the open positions from among the grid.
        % first find the remaining positions and then probabilistically
        % chose from among them. 
        pr3 = addchannel/totnumchan;
        
         count = 1;
         if isreal(addchannel)
         tobecomputedfield = zeros(1,addchannel);
         else
             tobecomputedfield = 0;
         end
    while count < length(b)
        if count > length(temp2)
            break
        else
        if (rand(1)<=pr3) %&& (~ismember(temp(count), temp2))
            %insert into temp2, the newly opened channel
            
            % alternately, 
            tobecomputedfield(1,count) = temp(count);
        else
            % do nothing
        end
        end
       
   count = count + 1; 
    end
sourcelowv = tobecomputedfield;% single ion channel position in ring
sourcehighv = sourcelowv+0.01;% phi-width of the channel is added to get the upper limit.
Inot = 2e-12/(pi*(((1e-6)^2)*3e-6)); % peak current strength at surface of node is 2 pico Amperes.
% The denominator contains a computation of the source cylindrical volume.
rho2prime = 2e-6; % highu, highu2, in microns - where the current source zeros out.
%For each (x,y) pair in the above matrix, compute the two integrals and
%store them as the px and py, respectively rho and z integrals. 
loops = 1;
returncellarray = cell(1,4);
field = cell(1,loops);
fieldzero = fieldzero.*ones(1,1);
for j = 1:loops
    field{1,j}(:,:) = fieldzero;
end



for channelnum = 0:windowsize:length(sourcelowv)-1
%     disp(channelnum)
for i = 1:mylength
    for j = 1:mylength
        [Q, Q2] = zandrhointegralsnew(sourcelowv(channelnum+1), sourcehighv(channelnum+1),pxnum(i,j),pi/2,pynum(i,j));
        p{i,j} = {Q2*((Inot/(rho2prime-rho1prime))/(4*pi*permittivity)),Q*(pynum(i,j)*(Inot/(rho2prime-rho1prime))/(4*pi*permittivity))};
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
vecxnum = windowsize.*cell2mat(vecx);
vecynum = windowsize.*cell2mat(vecy);

%axis tight manual;
% ax.NextPlot = 'replaceChildren';
%deltat = b/loops; % size of time steps

loops = 1;
% F(loops) = struct('cdata', [], 'colormap', []);
x = cell(1,loops);
y = cell(1,loops);
%field = cell(1,loops);
for j = 1:loops
        x{1,j} = pxnum;
        y{1,j} = pynum;
        field{1,j} = field{1,j} + real(sqrt(vecxnum.^2 + vecynum.^2));
end
% disp('hi3 from poremovie')
% size(pxnum)
% size(pynum)
%returncellarray = cell(1,3);
returncellarray{1,1} = real(x{1,1});
returncellarray{1,2} = real(y{1,1});
end
% oldans;
returncellarray{1,3} = oldans{1,3}(1,1) +real(field{1,1});
returncellarray{1,4} = sourcelowv;
end
end

