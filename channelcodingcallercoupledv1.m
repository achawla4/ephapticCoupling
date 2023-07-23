close all;
clear all;
clc;
clf;
% for p = 0.5;
 N = 200;
 x1 = 60;
 x2 = 60;
 p=0.45;
 ratio = .2;
[v1 min1 max1 v2 max1 max2] = mainaxonerrorexp(x1,x2,ratio)
 
% p = 0.65;
% q = 0.4;%BSC
% What is the capacity of the BSC?
% Cbsc = channelcapacityBSC(q);
k=20;
for n = 1:N
    %obtain the present q1 value
    q1 = v1(1,n);
    Cbsc1 = channelcapacityBSC(q1);
    q2 = v2(1,n);
    Cbsc2 = channelcapacityBSC(q2);


%Encoder
%Codeword Generation
% Suppose each codeword component is {0,1} as per Bernoulli with parameter p. 
% p lies between 0 and 1
% p = 0.9;
% Suppose the codeword is of length n
% n = 120;
% Suppose there are k codewords
% k = 7;
for j = 1:k
for i = 1:n
r1 = rand; % This gives a pseudorandom number between 0 and 1, uniformly distributed. 
% if this random number lies between (0,p) we count it as a 1 and if it
% lies between (p,1) we count it as a 0, thereby generating our Bernoulli
% samples. 
r2 = rand;
if r1<p
    bitvalue1(j,i) = 1;
else
    bitvalue1(j,i) = 0;
end
if r2<p
    bitvalue2(j,i) = 1;
else
    bitvalue2(j,i) = 0;
end
end
end

% % Display the codebook
% disp('codebook:')
% disp(bitvalue);
numtransmissions = 30;

%     for numtrials = 70
for m = 1:numtransmissions
sorf1(m) = channelcodingmainv3(p,q1,k,n,bitvalue1);
sorf2(m) = channelcodingmainv3(p,q2,k,n,bitvalue2);
end
numsuccesses1(n) = sum(sorf1);
numsuccesses2(n) = sum(sorf2);
successrate1(n) = numsuccesses1(n)/numtransmissions;
successrate2(n) = numsuccesses2(n)/numtransmissions;
end

% int1 = 1*p*10;
% % % % plot(1:20,numsuccesses)
% figure(floor(int1))
% for i = 1:floor(N/2)
%   subplot(floor(N/2),1,i)
%     plot([10,20,30,40,50, 60, 70],numsuccesses(i,[10,20,30,40,50, 60, 70], p*10), 'r*-');
% %     axis([0,70,0,1])
% end
% int2 = 2*p*10;
% figure(floor(int2))
% count = 1;
% for j = floor(N/2):N
%     while count <= floor(N/2)
%     subplot(floor(N/2),1,count)
%     plot([10,20,30,40,50, 60, 70],numsuccesses(j,[10,20,30,40,50, 60, 70], p*10),'r*-');
% %         axis([0,70,0,1])
%         count = count+1;
%     end
% end
% int3 = 3*p*10;
% disp(int3)
% figure(floor(int3))
% plot(1:N,numsuccesses(1:N))
% axis([0 100 0 50])
% 
% figure(floor(int3)+1)
% plot(1:N,successrate(1:N))
% 
% axis([0 N 0 1])

errorrate1 = 1-successrate1(1:N);
errorrate2 = 1-successrate2(1:N);
% errorexp = -1*log(errorrate)/n;
% hold on;
figure(16)
plot(1:N,errorrate1, 'r--');
hold on;
plot(1:N,errorrate2, 'b--');
figure(17)
plot(1:N, -log(errorrate1)./[1:N], 'r--');
hold on;
plot(1:N, -log(errorrate2)./[1:N], 'b--');

