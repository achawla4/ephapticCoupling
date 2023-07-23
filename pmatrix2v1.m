function outstring = pmatrix2v1(inputnum)
%pmatrix finder
% n=3 case
% inputnum is a string

% A. Give the Input

%000 input --> Trivial 000 output always
%001 input 
if inputnum=='001'

injectnumber1 = 1;
injectnumber2 = 1;
injectnumber3 = 1;
elseif inputnum=='010'

%010 input
injectnumber1 = 2;
injectnumber2 = 2;
injectnumber3 = 2;
elseif inputnum=='011'
%011 input
injectnumber1 = 1;
injectnumber2 = 2;
injectnumber3 = 1;
elseif inputnum=='100'

%100 input
injectnumber1 = 3;
injectnumber2 = 3;
injectnumber3 = 3;

elseif inputnum=='101'
%101 input
injectnumber1 = 1;
injectnumber2 = 3;
injectnumber3 = 1;
elseif inputnum=='110'
%110 input
injectnumber1 = 2;
injectnumber2 = 3;
injectnumber3 = 2;
elseif inputnum=='111'
%111 input
injectnumber1 = 1;
injectnumber2 = 2;
injectnumber3 = 3;
elseif inputnum=='000'
    v1 = zeros(100,1000);
    v2 = zeros(100,1000);
    v3 = zeros(100,1000);
end


% B. Compute the Response
if inputnum ~= '000'
[v1, v2, v3] = mymainaxong(injectnumber1, injectnumber2, injectnumber3)
end
% C. Use the response to find the output
% output is taken at node 80

% if node 80 of axon 1 crosses 100 mV, then take the output of axon 1 as 1
% and so on

myindex1 = find(v1>0.1, 1); % will return the first index when the voltage on axon 1 is above 0.1
if myindex1 > 200 % i.e. after some reasonable computation, not right at the beginning
    myoutput(1,1) = 1;
else
    myoutput(1,1) = 0;
end


% if node 80 of axon 2 crosses 100 mV, then take the output of axon 2 as 1
% and so on

myindex2 = find(v2>0.1, 1); % will return the first index when the voltage on axon 1 is above 0.1
if myindex2 > 200 % i.e. after some reasonable computation, not right at the beginning
    myoutput(1,2) = 1;
else
    myoutput(1,2) = 0;
end

% if node 80 of axon 3 crosses 100 mV, then take the output of axon 3 as 1
% and so on

myindex3 = find(v3>0.1, 1); % will return the first index when the voltage on axon 1 is above 0.1
if myindex1 > 200 % i.e. after some reasonable computation, not right at the beginning
    myoutput(1,3) = 1;
else
    myoutput(1,3) = 0;
end


% D. Store the input and output in a tabular form. 
outstring = strcat(num2str(myoutput(1,1)),num2str(myoutput(1,2)), num2str(myoutput(1,3)));








