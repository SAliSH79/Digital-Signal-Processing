%% Exp 1
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
yy = filter(b,a,y_n) ;

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

%% part 6


%% part 7
t3 = -5 : 0.01 : 5 ;
x = sinc(5.*t3).*sinc(5.*t3);
% plot(t3,x)
% grid on;
% axis([-5 5 -0.2 1.22])





%% FUnctions
function y = singen(w,n) 

    x1 = zeros(1,n);
    x(1) = 1;
    b=[0,sin(w)];
    a=[1,-2*cos(w),1];
    y = filter(b,a,x1);

end
