%% Exp 2 _ 1
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% part 2-1)a
clc;
a = [ 1 1 1 ];
h = [2 2 2 ];
l = 0:1:2;
figure(1)
subplot(311)
stem(l,a,'k-')
ylabel("Amp")
xlabel("time")
title("Original (a) Signal")
legend('signal')
grid on;
axis([0 5 0 2])

figure(1)
subplot(312)
stem(l,h,'r *')
ylabel("Amp")
xlabel("time")
title("h[n] Signal")
legend('filter')
grid on;
axis([0 5 0 2.5])

x = conv(a,h) ; %Using Matlab Function
b = myconv(a,h) ; %Using Manual Function
M = numel(a) + numel(h) - 1 ; %length of Conv
t = 0 : 1 : M - 1 ;
figure(1)
subplot(313)
stem(t,x,'*')
hold on
stem(t,b,'r')
ylabel("Amp")
xlabel("time")
title("Conv & MyConv Signal")
legend('conv' , 'myconv')
grid on;
axis([0 6 -1 7])
display("As You see they are equal in part a")
%% part 2-1)b
clc;
res = 1 ;
t = 0 : res : 199 ; %time
T = 50; %duty cycle
f = 1/T ;
x = square(2*pi*f*t) ; %generate signal
figure(2)
subplot(311)
stem(t,x,'b')
ylabel("Amp")
xlabel("time")
title("Square Wave Signal")
legend('square wave')
grid on;
axis([0 200 -1 1.3])

h = 0.1*ones(1,10) ; %Moving Average
l = 0 : 1 : 9; %time
figure(2)
subplot(312)
stem(l,h,'k *')
ylabel("Amp")
xlabel("time")
title("MA Signal")
legend('MA')
grid on;
axis([0 11 0 0.2])

h(end:numel(x)) = 0 ; %zero padding
y1 = conv(x,h) ;
y2 = myconv(x,h) ;
figure(2)
subplot(313)
stem(y1,'*')
hold on
stem(y2,'r')
ylabel("Amp")
xlabel("time")
title("Conv & MyConv Signal")
legend('conv' , 'myconv')
grid on;
axis([0 220 -1.2 1.2])
display("As You see they are equal in part b")

%% part 2-1)c
clc;
res = 1 ;
t = 0 : res : 199 ; %time
T = 50; %duty cycle
f = 1/T ;
x = square(2*pi*f*t) ; %generate signal
figure(3)
subplot(311)
stem(t,x,'b')
ylabel("Amp")
xlabel("time")
title("Square Wave Signal")
legend('square wave')
grid on;
axis([0 200 -1 1.3])

n = 0:1:14 ; %discrete time
h = 0.25*(0.75).^n ; %New Filter
l = 0 : 1 : 14; %time
figure(3)
subplot(312)
stem(l,h,'k -')
ylabel("Amp")
xlabel("time")
title("part c Filter Signal")
legend('Filter')
grid on;
axis([0 16 0 0.3])

h(end:numel(x)) = 0 ; %zero padding
y1 = conv(x,h) ;
y2 = myconv(x,h) ;
figure(3)
subplot(313)
stem(y1,'*')
hold on
stem(y2,'r')
ylabel("Amp")
xlabel("time")
title("Conv & MyConv Signal")
legend('conv' , 'myconv')
grid on;
axis([0 220 -1.2 1.2])
display("As You see they are equal in part c")

%% part 2-1)d
clc;
res = 1 ;
t = 0 : res : 199 ; %time
T = 50; %duty cycle
f = 1/T ;
x = square(2*pi*f*t) ; %generate signal
figure(4)
subplot(311)
stem(t,x,'b')
ylabel("Amp")
xlabel("time")
title("Square Wave Signal")
legend('square wave')
grid on;
axis([0 200 -1 1.3])

h1 = [1 -1]; %basic h[n]
h2 = conv(h1,h1) ; %1st conv
h3 = conv(h2,h1) ; %2nd conv
h4 = conv(h3,h1) ; %3rd conv
h5 = conv(h4,h1) ; %4th conv
hf = 0.2*conv(h5,h1) ; %5th conv = final impulse response
l = 0 : 1 : 6; %time
figure(4)
subplot(312)
stem(l,hf,'k o')
ylabel("Amp")
xlabel("time")
title("differentiator filter")
legend('differentiator')
grid on;
axis([0 7 -5 5])

hf(end:numel(x)) = 0 ; %zero padding
y1 = conv(x,hf) ;
y2 = myconv(x,hf) ;
figure(4)
subplot(313)
stem(y1,'*')
hold on
stem(y2,'r')
ylabel("Amp")
xlabel("time")
title("Conv & MyConv Signal")
legend('conv' , 'myconv')
grid on;
axis([0 220 -5 5])
display("As You see they are equal in part d")

%% MyconV
function y = myconv(h,x)
    lx = length(x) ;
    lh = length(h) ;
    x = [x,zeros(1,lx)]; %zero padding
    h = [h,zeros(1,lh)]; %zero padding
    y = zeros(1,lx + lh - 1);
    for i = 1 : lx + lh - 1 
        y(i) = 0 ;
        for j = 1 : lx + lh - 1 
            if(j < i + 1)
                y(i) = y(i) + x(j) * h(i - j + 1 );
            else
                break;
            end
        end
    end
end



