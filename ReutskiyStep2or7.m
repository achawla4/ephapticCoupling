function axon = ReutskiyStep2or7(leftmatrix, Dd, bundle1, axon, time,toggle)
% If toggle is 0, execute step 2, else execute step 7.
if toggle == 0
 VCN = triblocksolve(leftmatrix,Dd,bundle1.N);
       
       % update the voltages on all the axons.
       myhold = 0;
       advance = 1;
       var1 = 0;
       for count = 1:bundle1.N*bundle1.endlength
           
                if myhold == 1
                    advance = 0;
                    var1 = var1;
                else
                    advance = 1;
                    var1 = var1 + 1;
                end
               axnum = mod(count,bundle1.N);
               if axnum == 0
                   axnum = bundle1.N;
               end
               axon(axnum).V(var1,time+1)=VCN(count);
                 % release hold every time count becomes divisible by N. 
                 if mod(count,bundle1.N) == 0
                     myhold = 0;
                 else
                     myhold = 1;
                 end  
       end
else
     VCNprime = triblocksolve(leftmatrix,Dd,bundle1.N);
       
       % update the voltages on all the axons.
       myhold = 0;
       advance = 1;
       var1 = 0;
       for count = 1:bundle1.N*bundle1.endlength
           
                if myhold == 1
                    advance = 0;
                    var1 = var1;
                else
                    advance = 1;
                    var1 = var1 + 1;
                end
               axnum = mod(count,bundle1.N);
               if axnum == 0
                   axnum = bundle1.N;
               end
               axon(axnum).V(var1,time+1)=VCNprime(count);
                 % release hold every time count becomes divisible by N. 
                 if mod(count,bundle1.N) == 0
                     myhold = 0;
                 else
                     myhold = 1;
                 end  
       end
end