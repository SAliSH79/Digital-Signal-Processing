%% Exp 3 _ 1
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% part 3-1)a
R = 0.8 ;
f0 = 500 ;
fs = 10000;
w = linspace(0,pi,f0);
w0 = 2*pi*(f0/fs) ;
a1 = -2*R*cos(w0);
a2 = R^2 ;
G = (1 - R)*sqrt(1 - 2*R*cos(2*w0) + R^2);
num = [0 0 G];
den = [1 a1 a2];
[h,f] = freqz(num,den,1024,'whole',fs);
%% Plotting
clc;
figure(1)
plot((f/pi),(abs(h)).^2,"r")
grid on;
ylabel("Mag")
xlabel("Freq (Hz)")
title("| H(ejw) |")
axis([0 3150 -0.15 1.25])

%% part 3-1)b
clc;
n = 0 : 1 : 300 - 1 ;
h1 = G /sin(w0) ;
h2 = R.^n ;
h3 = sin(w0.*n + w0);
h = h1.*h2.*h3 ; %Impulse Responce

N = 300;
x = zeros(N,1);
x(1) = 1;
h_f = filter(num,den,x);
figure(2) ; plot(h_f,"b")
hold on; plot(h,"r")

y = zeros(1,N);
for i = 1:N
    if (i ==1)
        y(i) = G*1 ;
    elseif(i ==2)
        y(i) = -a1*y(i - 1) ;
    else
        y(i) = -a1*y(i - 1) -a2*y(i - 2) ;
    end
end
figure(2)
stem(y,"k*");
grid on;
ylabel("Amp")
xlabel("time")
title("Impulse Responce")
legend("filter syntax","Conv","Recursive")

%% 3 -1)c
s = cos(w0.*n);
v = randn(1,N);
x = s + v ;

y = zeros(1,N);
for i = 1:N
    if (i ==1)
        y(i) = G*x(i) ;
        w1 = y(1) ;
    elseif(i ==2)
        y(i) = -a1*w1 + G*x(i) ;
        w2 = w1;
        w1 = y(2) ;
    else
        y(i) = -a1*w1 -a2*w2 + G*x(i) ;
        w2 = w1;
        w1 = y(i);
    end
end  

figure(3)
plot(y,"g");
hold on;
plot(s,"r")
grid on;
ylabel("Amp")
xlabel("time")
legend("filter output","s[n]")
title("Plots")

Y = conv(y,x);
figure(4)
plot(Y,"k");
hold on;
plot(s,"r")
legend("Y","input")

%% 3-1)d
s = cos(w0.*n);
v = randn(1,N);
x = s + v ;

y1 = zeros(1,N);
for i = 1:N
    if (i ==1)
        y1(i) = G*v(i) ;
        w1 = y1(1) ;
    elseif(i ==2)
        y1(i) = -a1*w1 + G*v(i) ;
        w2 = w1;
        w1 = y1(2) ;
    else
        y1(i) = -a1*w1 -a2*w2 + G*v(i) ;
        w2 = w1;
        w1 = y1(i);
    end
end  
figure(5)
plot(v,"g");
hold on;
plot(y1,"r")
grid on;
ylabel("Amp")
xlabel("time")
legend("v[n]","filter output noise")
title("Noise Plots")

%% 3-1)e

s = std(y1) ;
s1 = (1 + R^2)/((1+R)*(1 + 2*R*cos(w0)+ R^2)) ;

display("a little Different is acceptable")



