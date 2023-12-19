%% Exp 1 _ the rest
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%% part1
close all ; clear ; clc;
le = 0.01 ; %step size
t = 0 : le : 2; %time domain
A = 5; %AMp
f = 1; %carrier frequency
y = A * sin(2*pi*f*t);  %Original signal 
figure(1)
stem(t,y)
plot(t,y,"r")
grid on;
ylabel("Amp")
xlabel("time")
title("Discrete Original Signal")
axis([0 2.5 -6 6])

figure(2)
stem(t,y,"o")
grid on;
ylabel("Amp")
xlabel("time")
title("Discrete Original Signal")
axis([0 2.5 -6 6])
%% part2
m = -.5 + (1) * rand(1,length(t)); %Random noise : [ a + (b - a)*rand(1,N) ] 
y_n = m + y ; %Noisy signal
% N1 = 21; 
% y_n = smooth(y_n, N1);

figure(3)
subplot(221)
stem(t,y)
grid on;
ylabel("Amp")
xlabel("time")
title("Discrete Original Signal")
axis([0 3 -6 6])

subplot(222)
stem(t,y_n)
grid on;
ylabel("Amp")
xlabel("time")
title("Discrete Signal with noise")
axis([0 2.5 -6 6])

subplot(223)
plot(t,y,"k")
grid on;
ylabel("Amp")
xlabel("time")
title("Discrete Original Signal")
axis([0 2.5 -6 6])

subplot(224)
plot(t,y_n,"k")
grid on;
ylabel("Amp")
xlabel("time")
title("Discrete Signal with noise ")
axis([0 2.5 -6 6])
%% part3
M1 = 0 ; 
M2 = 20 ;
MA_Len = M2 + M1 + 1 ; %length of MA
MA = ones(1,MA_Len)/ (MA_Len); % Moving Average Signal 
Wind_Sig = conv(MA,y_n); %Windowed Noisy signal
t2 = 0 : le : le*(numel(Wind_Sig) - 1) ; %Time domain for conv 

figure(4)
subplot(211)
plot(y_n,"k");
hold on;
plot(Wind_Sig,"b")
grid on;
ylabel("Amp")
xlabel("time")
title("Noisy SIgnal vs Windowed Signal")
axis([0 225 -6 6])
legend('Noisy' , 'By matlab smooth' , 'By conv')

subplot(212)
stem(y_n,"k");
hold on;
stem(Wind_Sig,"b")
grid on;
ylabel("Amp")
xlabel("time")
title("Noisy SIgnal vs Windowed Signal")
axis([0 225 -6 6])
legend('Noisy' , 'By conv')

%% part 4
a = 1 ;
b = ones(1,21) / 21;
yy = filter(b,a,y_n) ;  %%creating y using filter

figure(5)
plot(yy,'*')
hold on
plot(Wind_Sig)
grid on;
ylabel("Amp")
xlabel("time")
title("Filter SIgnal vs Windowed Signal")
axis([0 225 -6 6])
legend('MA' , 'Filter')

%% part5
% We created function at the end of code
y2 = singen(pi/6,300); %Output using Manual Function
t_5 = 1:1:300 ;
figure(6)
plot(t_5,y2)
grid on;
ylabel("Amp")
xlabel("time")
title("Singen Signal")


%% part 6
res = 0.01 ;
t6 = 0 : res : 4 ;
xt = cos(2*pi*t6) + cos(8*pi*t6) + cos(12*pi*t6) ; %Original signal
figure(7)
plot(t6,xt)
fs = 5 ; %sampling Frequency
T = (1/fs) ;
ts = 0 : T : 4; %sampling time
xs = cos(2*pi*ts) + cos(8*pi*ts) + cos(12*pi*ts) ; %Sampled signal
hold on;
plot(ts,xs,'r *')
grid on;
ylabel("Amp")
xlabel("time")
title("Original Signal vs Sampled Signal")
axis([0 5 -3 4])
% legend('Original Signal' , 'Sampled Signal')

wc = 2*pi*2 ;
recunsrtructedsig = zeros(1,numel(t6)); %PreAllocating
for n = 1 : numel(ts)
    recunsrtructedsig = recunsrtructedsig + xs(n)*(wc*T/pi)*(sinc(wc*(t6 - ts(n))/pi));
end
hold on;
plot(t6,recunsrtructedsig,'k')
grid on;
ylabel("Amp")
xlabel("time")
title("Original Signal vs Sampled Signal")
axis([0 5 -3 4])
legend('Original Signal' , 'Sampled Signal','Recunsrtructed Signal')

%% part 7
signalfunc = @(t) sinc(5*t).^2 ;
res = 0.01 ;
Fs = 1/res ;
t7 = -5 : res : 5 ;
x7 = signalfunc(t7) ;
figure(8)
plot(t7,x7,'k')
ylabel("Amp")
xlabel("time")
title("Original Signal")
grid on;
axis([-5 5 -0.2 1.22])
hold on;

L = length(x7) ;
spec_main = abs(fftshift(fft(x7))/L) ;
figure(9)
subplot(311)
plot(normalFreqArr(L) * Fs , spec_main)
ylabel("Amp")
xlabel("freq(hz)")
title("FFT Signal")
grid on;

fsampling7 = 4 ;
T7 = 1/fsampling7 ;
ts7 = -5 : T7 : 5 ;
xs7 = signalfunc(ts7); %sampling signal with fs = 4 Hz
figure(8)
plot(ts7,xs7,'o')
legend('Original' , ' sampled with fs')
hold on ;

Ls = length(xs7) ;
spec_sampled = abs(fftshift(fft(xs7))/Ls) ;
figure(9)
subplot(312)
plot(normalFreqArr(Ls) * fsampling7 , spec_sampled)
ylabel("Amp")
xlabel("freq(hz)")
title("FFT Sampled Signal")
grid on;

N = Ls ;
w = linspace(-4*pi,4*pi,4*N) ;
X = freqz(xs7,1,w)/N ;
subplot(313)
plot(w * fsampling7 /2*pi , abs(X) , 'r')
ylabel("Amp")
xlabel("freq(hz)")
title("FFT Sampled Signal")
grid on;

%% part 8
signalfunc = @(t) sinc(2*t) ;
res = 0.01 ;
Fs = 1/res ;
t8 = -5 : res : 5 ;
x8 = signalfunc(t8) ;
figure(10)
subplot(211)
plot(t8,x8,'k')
ylabel("Amp")
xlabel("time")
title("Sinc(2t)")
grid on;
axis([-5 5 -0.7 1.22])
hold on;

Fsamp = 25.5 *2 ; %sampling freq which we can change
Ts = 1 / Fs ;
ts = -5 : Ts : 5 ;
xs = signalfunc(ts) ;
figure(10)
subplot(212)
plot(ts,xs,'*')
ylabel("Amp")
xlabel("time")
title("Sampled Sinc(2t)")
grid on;
axis([-5 5 -0.7 1.22])
hold on;

L = length(x8) ;
spec_main = abs(fftshift(fft(x8))/L) ;
figure(11)
subplot(211)
plot(normalFreqArr(L) * Fs , spec_main)
ylabel("Amp")
xlabel("freq(hz)")
title("FFT Signal")
grid on;

Ls = length(xs) ;
spec_sampled = abs(fftshift(fft(xs))/Ls) ;
figure(11)
subplot(212)
plot(normalFreqArr(Ls) * Fsamp , spec_sampled , 'r')
ylabel("Amp")
xlabel("freq(hz)")
title("FFT Sampled Signal")
grid on;

%% part 9 a
res = 1/(2*pi) ;
Fs = 1/res ;
signalFunc = @(t,Fs) cos(pi*(1/ 16)*Fs*t ) + cos(pi*(5/ 16)*Fs*t ) + cos(pi*(9/ 16)*Fs*t ) + cos(pi*(13/ 16)*Fs*t );
t = -5 : res : 5 ;
xt = signalFunc(t,Fs) ;
figure(12);subplot(211) ; plot(t,xt,'r');grid on;ylabel("Amp");xlabel("time");
title("Original Signal");axis([-5 5 -4 4.22]);hold on;

Lt = length(xt) ;
spec_sampled = abs(fftshift(fft(xt))/Lt) ;
figure(12) ; subplot(212);hold on;
plot(normalFreqArr(Lt) * Fsamp , spec_sampled , 'k');ylabel("Amp");xlabel("freq(hz)");
title("FFT of Signal");grid on;axis([-25 25 -0.1 0.6])


h = xlsread('filters.xls',1) ; %analysis filter
f = xlsread('filters.xls',2) ; %synthesis filter
N = 256 ;
w = linspace(-pi,pi,N);
figure(13);plot(w/(2*pi),pwlog(freqz(h(1,:),1,w)/N) , 'r');ylabel("Mag");xlabel("freq(hz)");
title("Analysis filters");grid on;
hold on;
plot(w/(2*pi),pwlog(freqz(h(2,:),1,w)/N) , 'k');
plot(w/(2*pi),pwlog(freqz(h(3,:),1,w)/N) , 'b');
plot(w/(2*pi),pwlog(freqz(h(4,:),1,w)/N) , 'y');
legend('analysis1','analysis2','analysis3','analysis4')

figure(14);plot(w/(2*pi),pwlog(freqz(f(1,:),1,w)/N) , 'r');ylabel("Mag");xlabel("freq(hz)");
title("Synthesis filters");grid on;
hold on;
plot(w/(2*pi),pwlog(freqz(f(2,:),1,w)/N) , 'k');
plot(w/(2*pi),pwlog(freqz(f(3,:),1,w)/N) , 'b');
plot(w/(2*pi),pwlog(freqz(f(4,:),1,w)/N) , 'y');
legend('Synthesis1','Synthesis2','Synthesis3','Synthesis4')


coeff = [2 0 1 0.5] ; %% coefficients for frequencies
sub = zeros(4,numel(t));
for i = 1 : 1
   temp1 = filter(h(i,:) , 1, xt); %filter input with analythic filter
   figure(12);subplot(211);plot(t,temp1,'g'); 
   temp2 = downsample(temp1,4);
   L = length(temp2) ;
   spec_out = abs(fftshift(fft(temp2,L))/L);
   figure(12);subplot(212);plot(normalFreqArr(L) * Fs , spec_out , 'y')
   temp3 = coeff(i)*temp2 ;
   temp4 = upsample(temp3,4);
   temp5 = filter(f(i,:) , 1, temp4);
   temp5 = temp5(1 : length(t));
   figure(12);subplot(211);plot(t,temp5,'b');legend('original','analysis flitered','Synthesis Filtered')
   sub(i,:) = temp5;
end

xout = sum(sub) ;
spec_out = abs(fftshift(fft(xout,L))/L);
figure(12) ; subplot(212);hold on;
plot(normalFreqArr(L), spec_out , 'b');legend('original','analysis flitered','Synthesis Filtered')



%% singen
function y = singen(w,n) 
    x1 = zeros(1,n);
    x1(1) = 1;
    b=[0,sin(w)];
    a=[1,-2*cos(w),1];
    y = filter(b,a,x1);
end
%% normalFreqArr
function a = normalFreqArr(m)
     fnor = (0 : m - 1 ) / m ; 
     a = fftshift(fnor) ;
     a(a>=0.5) = a(a>=0.5) - 1 ; 
end
%% Power db
function a = pwlog(x) 
        a = 20*log(abs(x)) ;
end