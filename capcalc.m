close all; 
clear all;
clc;
clf;
for i=1:7
    if i==1
        inputnum='001'
matrix(i,:)= overallrunner(inputnum);
    elseif i==2
        inputnum='010'
        matrix(i,:) = overallrunner(inputnum);
    elseif i==3
             inputnum='011'
        matrix(i,:) = overallrunner(inputnum);
    elseif i==4
             inputnum='100'
        matrix(i,:) = overallrunner(inputnum);
    elseif  i==5
             inputnum='101'
        matrix(i,:) = overallrunner(inputnum);
    elseif i==6
             inputnum='110'
        matrix(i,:) = overallrunner(inputnum);
    elseif i==7
             inputnum='111'
        matrix(i,:) = overallrunner(inputnum);
    end
end


[c r]  = balgoaliter(matrix)