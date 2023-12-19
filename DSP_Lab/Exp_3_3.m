% Exp 3 _ 3
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% part 3_3)a
%% Initialization
w0 = 0.15*pi ;
n = 1:600 ;%sample
N = numel(n);
x1 = cos(w0*n);
xx = zeros(1,N);
tt =(n-1)/200; %time
n1 = 1:200;
n2 = 201 : 400;
n3 = 401 : 600;
xx(n1) = 2*cos(w0*(n1 - 1)) ;
xx(n2) = 4*cos(w0*(n2 - 1)) ;
xx(n3) = 0.5*cos(w0*(n3 - 1)) ;

%% Plotting
figure(1);
plot(tt,xx,"c");
grid on;
xlabel('time(n)'); 
ylabel('x(n)'); 
title('Input Signal');
axis([0 3 -4.5 4.5])

%% Control
clc;
c = zeros(1,N);
w1 = 0;
landa = 0.9;
for i = 1 : N
     c(i) = landa*w1 + (1 - landa)*abs(xx(i)) ;
     w1 = c(i) ;
end

f = zeros(1,N);
c0 = 0.5 ;
p = 0.2;

for i = 1 : N
   if c(i) >= c0  
      f(i) = (c(i)/c0) ^ (p - 1) ;
   else
      f(i) = 1;
   end
end

%% Plotting
figure(2)
subplot(211)
plot(n,f,'b')
ylabel("Amp")
xlabel("Sample")
title("Gain Signal")
axis([0 600 0 1.25])
grid on;
subplot(212)
plot(n,c,'m')
grid on;
ylabel("Amp")
xlabel("Sample")
title("Controler Signal")
axis([0 600 0 3.25])

%% Output Compressor
y = xx.*f ;
figure(3)
plot(n,y,'k')
grid on;
ylabel("Amp")
xlabel("Sample")
title("Compressor Output Signal")
axis([0 600 -3.25 3.25])

%% Moving Average

