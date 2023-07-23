function sorf = channelcodingmainv2(p,q,k,n,bitvalue)
% clf;

%Message


epsilon = 0.02; %0.02;

% what is the rate of the above code? 
R = log2(k)/n


    

%Signal

%Channel

%Let the channel be a binary symmetric channel with cross over probability
%q. If the input is 0, it becomes a 1 with probability q and stays a 0 with
%probability 1-q. If the input is 1, it becomes a 0 with probability q and
%stays a 1 with probability 1-q. 

% Thus Pr(1|0) = q = Pr(0|1) and Pr(1|1) = 1-q = Pr(0|0).
% Thus the transition matrix has 1-q on the diagonal and q on the
% off-diagonal.

% q = 1e-4; 

% If the channel input is 0, the output is 1 with probability q and 0 with
% probability 1-q. Thus the output is Bernoulli(q). 

% If the channel input is 1, the output is 0 with probability q and 1 with
% probability 1-q. Thus the output is Bernoulli(1-q). 

% Generate the channel output bit, based on the channel input bit
% Suppose cin = channel input bit

%Suppose the f-th codeword is sent in where f is uniformly (mass points)
%distributed over {1,2,3,4,5,6,...k}.
% fix k = 7 for now.
rf= rand;

for i = 1:k

if (rf>(i-1)/k)&&(rf<(i/k))
    codewordindex = i;
end
end

% if (rf<1/7)&&(rf>0)
%     %choose the 1st codeword
%     codewordindex = 1;
% elseif (rf<2/7)&&(rf>1/7)
%     codewordindex = 2;
% elseif (rf<3/7)&&(rf>2/7)
%     codewordindex = 3;
% elseif (rf<4/7)&&(rf>3/7)
%     codewordindex = 4;
% elseif (rf<5/7)&&(rf>4/7)
%     codewordindex = 5;
% elseif (rf<6/7)&&(rf>5/7)
%     codewordindex = 6;
% elseif (rf>6/7)
%     codewordindex = 7;
% end    

disp('codeword number transmitted:')
disp(codewordindex)
for coutconstructor = 1:k
    for m = 1:size(bitvalue,2)

cin = bitvalue(coutconstructor,m);
r0 = rand; % This gives a pseudorandom number between 0 and 1, uniformly distributed. 

if cin == 0


% if this random number lies between (0,p) we count it as a 1 and if it
% lies between (p,1) we count it as a 0, thereby generating our
% Bernoulli(p) samples. 
if r0<q
    bitvaluem = 1;
%     disp('Flipped: 0 to 1')
else
    bitvaluem = 0;
%     disp('Error free')
end

% cout(coutconstructor,m)=bitvaluem;

else % cin == 1

% r1 = rand; % This gives a pseudorandom number between 0 and 1, uniformly distributed. 
% if this random number lies between (0,p) we count it as a 1 and if it
% lies between (p,1) we count it as a 0, thereby generating our
% Bernoulli(p) samples. 
if r0<q
    bitvaluem = 0;
%     disp('Flipped 1 to 0')
else
    bitvaluem = 1;
%     disp('Error free')
end



end
cout(coutconstructor,m)=bitvaluem;
    end
end

%Show the channel output
% disp('cout:')
% disp(cout)

%Decoding
%a. Jointly typical decoding

% We must figure out which message cout is jointly typical with. 
% First we must find the actual entropies of the input and output of the
% channel. 
% The channel input entropy is binary with parameter p
HX=-p*log2(p)-(1-p)*log2(1-p);
% The channel output entropy is to be found.
p00=(1-q)*(1-p);
p01=q*(1-p);
p10 = q*p;
p11 = p*(1-q);
py1= p01+p11;
py0= p00+p10;
HY=-py0*log2(py0)-py1*log2(py1);
HXY=-p00*log2(p00)-p01*log2(p01)-p10*log2(p10)-p11*log2(p11);
mibetinandout = HX + HY - HXY;


% This will allow the calculation of the joint entropy of the input and
% output.

% Next we need to find the empirical entropies of the channel input and
% output sequences
% Empirical entropy of channel input sequence
for j = 1:k
    codewordtocheck = bitvalue(j,:);
num1s = sum(codewordtocheck); % number of 1s in the sequence
num0s = n- num1s;
empprobinput = (p^(num1s))*((1-p)^(num0s));
empentropyinput(j,1)= (-1/n)*log2(empprobinput);


% Empirical entropy of channel output sequence
% Count number of 1s in the sequences
num1s_o = sum(cout(codewordindex,:));
num0s_o = n- num1s_o;
empproboutput = ((py1)^(num1s_o))*((py0)^(num0s_o));
empentropyoutput = (-1/n)*log2(empproboutput);


% Empirical joint entropy of channel input and output sequence
% Count number of 1s in the sequences
% num1s_j = num1s + num1s_o;
% probof1s_j = num1s_j/(2*n)
% probof0s_j = (2*n-num1s_j)/(2*n);
% empiricalentropyonput=  -(1/(2*n))*log(probof1s_j);

runningprod = 1;

for i = 1:n
    if (bitvalue(j,i)==0)&&(cout(codewordindex,i)==0)
        runningprod = runningprod*p00;
    elseif (bitvalue(j,i)==1)&&(cout(codewordindex,i)==0)
        runningprod = runningprod*p10;
    elseif (bitvalue(j,i)==0)&&(cout(codewordindex,i)==1)
        runningprod = runningprod*p01;
    else
        runningprod = runningprod*p11;
    end
end
% disp(runningprod)
empjointentropyinputoutput(j,1) = (-1/n)*log2(runningprod);

empmi = empentropyinput + empentropyoutput - empjointentropyinputoutput;


end

a=abs(-HY+empentropyoutput);
b=abs(-HX+empentropyinput);
c=abs(-HXY+empjointentropyinputoutput);
decoding = 0;
if a < epsilon
    for i = 1:k
        if (b(i)<epsilon)&&(c(i)<epsilon)&&(decoding==0)
            decoding = i;
        end
    end
else
    disp('Decoding Failure');
end

if decoding == 0
    disp('Decoding Failure')
end
if decoding ~= codewordindex
    sorf = 0;
else
    sorf = 1;
end

% 
% 
% L = floor(n/2);
% [mina(1,:),mindexa(1,:)] = mink(a(1,:),L);
% 
% for j = 1:k
%    
% [minb(j,:),mindexb(j,:)] = mink(b(j,:),L);
% [minc(j,:),mindexc(j,:)] = mink(c(j,:),L);
% end
% returnedmina = zeros(k,L);
% returnedmindexa = zeros(k,L);
% returnedminb = zeros(k,L);
% returnedmindexb = zeros(k,L);
% returnedminc = zeros(k,L);
% returnedmindexc = zeros(k,L);
% 
% 
% % in mina find the first 'delta' transition
% if ~isempty(find(diff(mina(1,:))>0))
%     indices(1) = find(diff(mina(1,:))>0)
%     
%  returnedmina(1,:) = mina(1,indices(1));
%  returnedmindexa(1,:) = mindexa(1,indices(1));
% else
%     returnedmina(1,:) = zeros(1,L)+a;
%     returnedmindexa(1,:) = 1:L;
% end
% 
% 
% for inputnum = 1:k
% % in mina find the first 'delta' transition
% if ~isempty(find(diff(minb(inputnum,:))>0))
%     indices2(inputnum) = find(diff(minb(inputnum,:))>0)
%     
%  returnedminb(inputnum,:) = minb(inputnum,indices2(inputnum));
%  returnedmindexb(inputnum,:) = mindexb(inputnum,indices2(inputnum));
% else
%     returnedminb(inputnum,:) = zeros(1,L)+b(inputnum);
%     returnedmindexb(inputnum,:) = 1:L;
% end
% % in mina find the first 'delta' transition
% if ~isempty(find(diff(minc(inputnum,:))>0))
%     indices3(inputnum) = find(diff(minc(inputnum,:))>0)
%     
%  returnedminc(inputnum,:) = minc(inputnum,indices3(inputnum));
%  returnedmindexc(inputnum,:) = mindexc(inputnum,indices3(inputnum));
% else
%      returnedminc(inputnum,:) = zeros(1,L)+c(inputnum);
%     returnedmindexc(inputnum,:) = 1:L;
% end
% 
% end
% 
% if (size(returnedmindexa,2)==1)&&(size(returnedmindexb,2)==1)&&(size(returnedmindexc,2)==1)
% 
% if returnedmindexa == returnedmindexb
%     if returnedmindexb == returnedmindexc
%         disp('JT Decoding proceeded - 1')
%         decoding = returnedmindexa(codewordindex)
%     else
%         disp('JT Error - 1')
%         decoding = k+1
%     end
% else
%     disp('JT Error - 2')
%     decoding = k+1
% end
% elseif size(returnedmindexa,2)>1
% %     %call weakdecoding
% %     output = weakdecoding(bitvalue,cout, codewordindex, returnedmindexa);
% %     d1 = sqeucdist(cin(codewordindex,:),output(1));
% %     d2 = sqeucdist(cin(codewordindex,:))
%     if returnedmindexb == returnedmindexc
%         disp('JT Decoding proceeded - 2')
%         decoding = returnedmindexb(codewordindex)
%     else
%         disp('JT Error - 3')
%         decoding = k+1
%     end
% elseif size(returnedmindexb,2)>1
%     %call weakdecoding and find the means of the 2 clusters.
%     output = weakdecoding(bitvalue,cout, codewordindex, returnedmindexb);
%     % for the received word, find which mean it is closer to
% 
%     if hammingdist(output(1),cout(codewordindex,:))>hammingdist(output(2), cout(codewordindex,:))
%         % received word belongs to cluster 2
%         % Decode to that codeword which is closest to output(2)
%         for i = 1:k
%         d(i)=hammingdist(bitvalue(i,:),output(2));
%         end
%         [mind, mindexd] = min(d);
%         decoding = mindexd
%    
%         disp('Weak Decoding proceeded - 3')
%      
%     else
%      for i = 1:k
%         d2(i)=hammingdist(bitvalue(i,:),output(1));
%         end
%         [mind2, mindexd2] = min(d2);
%         decoding = mindexd2
%    
%         disp('Weak Decoding proceeded - 4')
%     end
% else
% %     %size(returnedmindexc)>1
% %     % call weakdecoding
% %     output = weakdecoding(bitvalue,cout, codewordindex);
%     if returnedmindexa == returnedmindexb
%         disp('Decoding proceeded - 5')
%         decoding = returnedmindexa(codewordindex)
%     else
%         disp('Error - 5')
%         decoding = k+1
%     end
% end
% 
% if decoding == codewordindex
%     sorf = 1;
%     disp('Transmission Success');
% else
%     sorf = 0;
%     disp('Transmission Failure');
% end


% sizeoftypsetupperbnd = 2^(n*HXY)



