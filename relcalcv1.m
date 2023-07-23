close all; 
clear all;
clc;
clf;
for i=1:8
    if i == 1
        inputnum = '000'
        matrix(i,:)=overallrunner2v1(inputnum);
    elseif i==2
        inputnum='001'
matrix(i,:)= overallrunner2v1(inputnum);
    elseif i==3
        inputnum='010'
        matrix(i,:) = overallrunner2v1(inputnum);
    elseif i==4
             inputnum='011'
        matrix(i,:) = overallrunner2v1(inputnum);
    elseif i==5
             inputnum='100'
        matrix(i,:) = overallrunner2v1(inputnum);
    elseif  i==6
             inputnum='101'
        matrix(i,:) = overallrunner2v1(inputnum);
    elseif i==7
             inputnum='110'
        matrix(i,:) = overallrunner2v1(inputnum);
    elseif i==8
             inputnum='111'
        matrix(i,:) = overallrunner2v1(inputnum);
    end
end

% call channelcodingcaller modified to use above matrix instead of BSC(q)

successrate = channelcodingcallertractver1(matrix);

% [c r]  = balgoaliter(matrix)