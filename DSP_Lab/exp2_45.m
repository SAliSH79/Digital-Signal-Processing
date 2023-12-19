%% Exp 2 _ the rest
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%% Clear Recent Data
close all ; clear ; clc;
%% Initialization
clc;
fs = 1000; %Sampling Frequency
t = 0 : 1/fs : 2 - 1/fs ; %time
x1 = sin(2*pi*100*t) ;
x2 = chirp(t,200,2,400) ;
sig = x1 + x2 ;
sig(250) = sig(250) + 50 ; %Impulse in 250th sample
figure(1)
plot(t,sig)
ylabel("Amp")
xlabel("time")
title("Original Signal")
legend('Mixed Sig signal')
grid on;
axis([0 2 -3.0 55])

%% part 2-4)a
clc;
NFFT = 1024*3 ;
f = normalFreqArr(NFFT);
f = fs*f;
spec1 = abs(fftshift(fft(sig,NFFT))/NFFT) ;

figure(2)
plot(f,spec1,"r")
ylabel("Amp")
xlabel("Freq")
title("FFT of Mixed Signal")
legend('Mixed FFT signal')
grid on;
axis([-500 500 0 0.45])

spec2 = abs(fftshift(fft(x2,NFFT))/NFFT) ;
spec3 = abs(fftshift(fft(x1,NFFT))/NFFT) ;
figure(3)
subplot(211)
plot(f,spec2,"g")
ylabel("Amp")
xlabel("Freq")
title("FFT of X2 Chirp Signal")
legend('FFT of X2 chirp signal')
grid on;
axis([-500 500 0 0.1])

figure(3)
subplot(212)
plot(f,spec3,"b")
ylabel("Amp")
xlabel("Freq")
title("FFT of X1 Sin Signal")
legend('FFT of X1 sin signal')
grid on;
axis([-500 500 0 0.5])

%% Adding a Pulse
n = 1.2*fs : 1 : 1.6*fs -1 ;
x3 = 10*sin(2*pi*50*(n/fs));
x33 = zeros(1,numel(x1));
x33(n) = x3 ;
sig1 = sig + x33 ;

figure(3)
plot(t,sig1)
ylabel("Amp")
xlabel("time")
title("Sig Signal")
legend('Mixed Sig signal')
grid on;
% axis([0 2 -35 55])

spec1 = abs(fftshift(fft(sig1,NFFT))/NFFT) ;
figure(4)
plot(f,spec1,"g")
ylabel("Amp")
xlabel("Freq")
title("FFT of Mixed Signal")
legend('FFT Mixed signal')
grid on;
% axis([-500 500 -.010 0.45])
%% part 2-4)b
clc;
N = 100;
win = hamming(N) ;
NOverLap = N - 10;
[y,fArr,timeArr,p] = spectrogram(sig1,win,NOverLap,NFFT,fs);
figure(5);spectrogram(sig1,win,NOverLap,NFFT,fs,'xaxis')
title("Spectrogram of Mixed Signal")
grid on;

figure(6)
surf(timeArr,fArr,10*log10(abs(p)),'Edgecolor','none')
axis xy;
axis tight;
colormap(jet);
view(0,90)
ylabel("Freq")
xlabel("time")
title("WideBand Spectrogram")
%% part 2-5)
clc;
N = 20 ;
index = linspace(0,1,2^N);
[X,XN] = wnoise('quadchirp',N,3.5);
figure(7)
subplot(211)
plot(index,X);
subplot(212)
plot(index,XN);

[cA1,cD1] = dwt(XN,'db1');
[r1,c1] = size(cA1) ;
ind1 = linspace(0,1,c1);
figure(8)
subplot(211)
plot(ind1,cA1,'k')

[cA2,cD2] = dwt(cA1,'db1');
[r2,c2] = size(cA2) ;
ind2 = linspace(0,1,c2);
figure(8)
subplot(212)
plot(ind2,cA2,'k')

Xden = wdenoise(XN);
figure(9)
plot(index,Xden);

%% normalFreqArr
function a = normalFreqArr(m)
     fnor = (0 : m - 1 ) / m ; 
     a = fftshift(fnor) ;
     a(a>=0.5) = a(a>=0.5) - 1 ; 
end