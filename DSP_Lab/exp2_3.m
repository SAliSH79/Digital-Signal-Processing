%% Exp 2 _ the rest
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%% part1
close all ; clear ; clc;
%% 2 - 3) a Voice Scrambling
clc;
clear;close all;
recObj = audiorecorder; %create an audiorecorder object with the desired properties
fs = 44100 ; %sampling frequency
disp('Start speaking.')
recordblocking(recObj,10);
disp('End of Recording.');
y = getaudiodata(recObj);
%% plotting
y1 = y' ;
figure(1)
plot(y1,'r');
title('Audio Signal (double)');
grid on;
ylabel("Amp")
xlabel("time")
%% Playing
clc;
play(recObj);
%% Saving
filename = 'Audio01.wav';
audiowrite(filename, y, fs); %saving file
%%
clc;
fs = 44100;
[y,Fs] = audioread(filename);
sound(y,fs);

%%
figure(2)
y2 = y' ;
plot(y2,'b');
title('Audio Signal (double)');
grid on;
ylabel("Amp")
xlabel("time")

%% 2 - 3 b ) %low pass filter
clc;
filter1 = load('filter.mat') ;
%% 2 - 3 c ) filttering 
clc;
%filtering signal
% sig1 = filter(filter1.filter,1,y2) ; % low pass filter
sig1 = filter(coef,1,y2);
%% plotting
clc;
figure(3)
plot(sig1,'b');
title('Audio Filtered (sig1) Signal');
grid on;
ylabel("Amp")
xlabel("time")

%% carrier Freq
f0 = fs /4 ; %carrier Freq
f0 = f0 - 10 ;
n = 1 : 1 : numel(sig1) ;
s = 2*cos(2*pi*(f0/fs).*n); %carrier signal
sig2 = sig1.*s ;%signal * carrier
%% plotting
clc;
figure(4)
plot(sig2,'b');
title('Audio Filtered (sig1)*carrier Signal');
grid on;
ylabel("Amp")
xlabel("time")

%% filterring 
sig3 = conv(coef,sig2) ;% low pass filter again
%% plotting
clc;
figure(4)
plot(sig3,'b');
title('Audio Again Filtered (sig3) Signal');
grid on;
ylabel("Amp")
xlabel("time")
%% playing
sound(10*sig3,fs);

%%  2- 3 d )
sigg = conv(coef,sig3) ;%filtering signal
clc;
f0 = fs /4 ; %carrier Freq
f0 = f0 - 10 ;
n = 1 : 1 : numel(sigg) ;
s = 2*cos(2*pi*(f0/fs).*n); %carrier signal
sigg2 = sigg.*s ; %signal * carrier
sigg3 = conv(sigg2,coef); %filtering signal
%% plotting
clc;
figure(5)
plot(sigg3,'b');
title('Audio Returned Filtered (sigg3) Signal');
grid on;
ylabel("Amp")
xlabel("time")
%% playing
sound(10*sigg3,fs);