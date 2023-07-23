function elevation = randomtrapezoid(fieldzero, start,width,height,b,xpos)
%%RANDOMTRAPEZOID
%% This function creates a trapezoid of input height, random 
%% top side length and random position along the bottom horizontal axis
%% by randomly modulating the width of the underlying trapezoid and placing
%% it randomly on the horizontal axis.
%%% height is the height of the trapezoid
%%% a is the length of the shorter, top side
%%% b is the length of the longer, bottom side
%%% xpos is the position along the bottom side that one 
%%% desires the height for.


%% Based on the random width, generate the top side length
aprime = width;
%% Generate the bottom side length
bprime  = 2*width;

m1 = height/((bprime-aprime)/2);
m2 = -1*m1;
%% Thus we have the features of the new trapezoid starting at start. They are 
%% : aprime, bprime and height, along with m1 and m2. We will need:
y = height/m1;

if (xpos>=0)&&(xpos<=start)
    elevation = fieldzero; 
elseif (xpos>start)&&(xpos<=(start + y))
    elevation = m1*(xpos-start);
elseif (xpos>(start + y))&&(xpos<=(start+y+aprime))
    elevation = height;
elseif (xpos>(start+y+aprime))&&(xpos<=(start+bprime))
    elevation = m1*((start+bprime)-xpos);
elseif (xpos>(start+bprime)) && (xpos<=b)
    elevation = fieldzero;
end