% Exp 3 _ 2
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% part 3_2)a
%% Initialization
f1 = 4;
f2 = 8;
f3 = 12;
fs = 400; %sampling Frequency
t1 = 0 : 1/fs : 2 - 1/fs ; %time1
t2 = 2 : 1/fs : 4 - 1/fs ; %time2
t3 = 4 : 1/fs : 6 - 1/fs ; %time3
f0 = 8;
delta_f = 4;
w0 = 2*pi*f0 / fs;
delta_w = 2*pi*delta_f/fs ;
beta = tan(delta_w/2);
%% Input
x1 = cos(2*pi*f1*t1);
x2 = cos(2*pi*f2*t2);
x3 = cos(2*pi*f3*t3);
n = 1:6*fs ;%sample
N = numel(n);
xx = zeros(1,N);
tt =(n-1)/fs; %time
n1 = 1:2*fs;
n2 = 2*fs + 1 : 4*fs;
n3 = 4*fs + 1 : 6*fs;
% n4 = 3*fs + 1 : 5*fs;
xx(n1) = cos(2*pi*f1*(n1 - 1)/fs) ;
xx(n2) = cos(2*pi*f2*(n2 - 1)/fs) ;
xx(n3) = cos(2*pi*f3*(n3 - 1)/fs) ;
% xx(n4) = xx(n4) + 2*sin(2*pi*f1*(n4 -1)/fs);
%%
figure(1);
plot(tt,xx,"r");
grid on;
xlabel('time(n)'); 
ylabel('x(n)'); 
title('Input Signal');
axis([0 6.5 -1.5 1.5])
%% Filter with delta_f = 4 Hz
clc;
b1 = -2*cos(w0);
b2 = 1;
num = [1 b1 b2];
a1 = -2*cos(w0);
a2 = 1-beta;
den = [1+beta a1 a2];
Gp = tf (num , den);
GpssObs = canon(Gp,'companion') %Canonical Presentation
GpssObsA = GpssObs.A
GpssObsB = GpssObs.B
GpssObsC = GpssObs.C
GpssObsD = GpssObs.D

%% Frequency Responce
figure(2)
freqz(num,den,1024,'whole')
title('Filter Frequnecy Responce')
%% 
w = 0 :0.001: 2*pi;
h = freqz(num,den,w);
figure()
plot(w*fs/(2*pi) ,abs(h),"k")
grid on;

%% Buffer with Delay
clc;
N = 300;
y = zeros(1,N);
a1 = 2*cos(w0)/(1+beta);
a2=(1-beta)/(1+beta);
b1 = (-2*cos(w0))/(1+beta) ;
b2 = 1/(1+beta) ;
for i = 1:N
    if (i ==1)
        y(i) = (1/(1+beta))*xx(i) ;
        w1 = y(1) ;
        m1 = xx(1);
    elseif(i ==2)
        y(i) = -a1*w1 + (1/(1+beta))*xx(i) + b1*m1 ;
        w2 = w1;
        w1 = y(2) ;
        m2 = m1;
        m1 = xx(2);
    else
        y(i) = -a1*w1 -a2*w2 + (1/(1+beta))*xx(i) + b1*m1 + b2*m2 ;
        w2 = w1;
        w1 = y(i);
        m2 = m1;
        m1 = xx(i);
    end
end  
%% Plotting
figure(3);
ty = linspace(0,6,length(y));
plot(ty,y ,'k');
grid on;
xlabel('time(n)'); 
ylabel('y(n)'); 
title('Filtered with H1(z) Signal');

%% 3-2)b
clc;
figure(4)
step(Gp,"b")
grid on;
xlabel('time(n)'); 
ylabel('h(n)'); 

%% 3-2)c
y1 = filter(num,den,xx);
figure(5)
ttt = linspace(0,6,length(y1));
plot(ttt,y1,"k");
xlabel('time(n)'); 
grid on;
ylabel('y(n)'); 
title('x(n) filtered by H1(z)');