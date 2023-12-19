%% Exp 2 _ 2
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% Initialization
w1 = 0.05 * pi ;
w2 = 0.2 * pi ;
w3 = 0.35 * pi ;
M = 100 ; %degree of filter
n = 0 : 1 : M-1 ; %time
s = sin(w2.*n); %original signal
v = sin(w1.*n) + sin(w3.*n); %Interference Signal
x = s + v ; %recived signal
wa = 0.15*pi ;
wb = 0.25*pi ;
t1 = n - M/2 ;
h1 = sinc(wa.*t1) ; %first half
h2 = sinc(wb.*t1) ; %second half
hp = h2 - h1 ; %BandPass filter

W = 0.54 - 0.46.*(cos(2*pi.*n/M)) ; %hamming window
h = W.*hp ; %final windowed bandpass filter

%% part 2-2)a
figure(1)
stem(n,x,'*')
hold on;
stem(n,s,'k o')
ylabel("Amp")
xlabel("time")
title("x[n] & s[n] Signal")
legend('Recieved signal','Noisy signal')
grid on;
axis([0 105 -5 5])

%% part 2-2)b
y1  = filter(h,1,x) ; %filter x using Window bandpass filter
figure(2)
subplot(211)
stem(n,s,'k o')
ylabel("Amp")
xlabel("time")
title("s[n] Signal")
legend('Interference signal')
grid on;
axis([0 110 -1.5 1.5])

subplot(212)
stem(n,y1,'b o')
ylabel("Amp")
xlabel("time")
title("Filtered Signal")
legend('Filtered signal')
grid on;
axis([0 110 -2.5 2.5])

%% part 2-2)c



