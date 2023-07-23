function matrix= overallrunner2v1(inputnum)
%overall runner
% inputnum='001'
N = 10; % number of trials
for i = 1:N
%     disp(i)
outstring = pmatrix2v1(inputnum);
if outstring == '000'
outstringmy(i) = 0;
elseif outstring == '001'
    outstringmy(i) = 1;
elseif outstring == '010'
    outstringmy(i) = 2;
elseif outstring == '011'
    outstringmy(i)  = 3;
elseif outstring == '100'
    outstringmy(i) = 4;
elseif outstring == '101'
    outstringmy(i) = 5;
elseif outstring == '110'
    outstringmy(i) = 6;
elseif outstring == '111'
    outstringmy(i) = 7;
end
end

indices1 = find(outstringmy == 0);
indices2 = find(outstringmy== 1);
indices3 = find(outstringmy == 2);
indices4 = find(outstringmy== 3);
indices5 = find(outstringmy == 4);
indices6 = find(outstringmy== 5);
indices7 = find(outstringmy == 6);
indices8 = find(outstringmy== 7);

matrix(1,1) = length(indices1)/N;
matrix(1,2) = length(indices2)/N;
matrix(1,3) = length(indices3)/N;
matrix(1,4) = length(indices4)/N;
matrix(1,5) = length(indices5)/N;
matrix(1,6) = length(indices6)/N;
matrix(1,7) = length(indices7)/N;
matrix(1,8) = length(indices8)/N;


