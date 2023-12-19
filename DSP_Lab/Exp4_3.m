% Exp 4 _ 3
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%% clear recent data
clc;
close all;
clear;
%% 4-3)a


%% 4-4)a 
clc;
theta = 20;
len = 10;
h = fspecial('motion', len, theta) ;
f=imread('flower.jpeg');
f= im2double(f);
G = imfilter(f,h,'conv','circular');

figure(1);
subplot(2,2,[1, 3])
imshow(f);
grid on; 
title('Flower');

subplot(2,2,[2,4])
imshow(G);
grid on; 
title('Flower Motion Blured filter ');

%% 4-4)b
clc;
NSR = zeros(1,6);
wnr1 = deconvwnr(G,h) ;
figure(2)
subplot(231)
imshow(wnr1)
title('Restored Flower Motion Blured with NSR = 0')

% Noise Parameters
noise_mean = 0;
noise_var = 0.0001;
% Estimated NSR
signal_var = var(f(:));
NSR(2) = noise_var / signal_var;

wnr3 = deconvwnr(G,h,NSR(2));
figure(2)
subplot(232)
imshow(wnr3)
title('Restoration of Flower Motion Blured Image (Estimated NSR)')

% NSR random
NSR(3 : 6) = [0.5 , 1, 2 , 5 ];
%
clc;
for i = 2 : 6
    wnr = deconvwnr(G,h,NSR(i));
    figure(2)
    subplot(2,3,i)
    imshow(wnr)
    title(['Restoration of Flower Motion Blured  Image NSR = ', ...
        num2str(NSR(i))]) 
end

%% 4-4)c
% Noise Parameters
noise_mean = 0;
noise_var = 10;
blurred_noisy = imnoise(G,'gaussian',noise_mean,noise_var);
figure(3)
imshow(blurred_noisy)
title('Blurred and Noisy Flower Image')

%% 4-4)d
clc;
NSR = zeros(1,6);
wnr1 = deconvwnr(blurred_noisy,h) ;
figure(4)
subplot(231)
imshow(wnr1)
title('Restored Flower Motion Noisy Blured with NSR = 0')

% Estimated NSR
signal_var = var(f(:));
NSR(2) = noise_var / signal_var;

wnr3 = deconvwnr(blurred_noisy,h,NSR(2));
figure(4)
subplot(232)
imshow(wnr3)
title('Restoration of Flower Motion Noisy Blured Image (Estimated NSR)')

% NSR random
NSR(3 : 6) = [0.5 , 1, 2 , 5 ];
%
clc;
for i = 2 : 6
    wnr = deconvwnr(blurred_noisy,h,NSR(i));
    figure(4)
    subplot(2,3,i)
    imshow(wnr)
    title(['Restoration of Flower Motion Noisy Blured  Image NSR = ' ...
        ,num2str(NSR(i))]) 
end

%% 4-5)a
clc;
I = imread('glass.tif');
I1 = im2double(I);

figure(5);
subplot(311)
imshow(I1);
grid on; 
title('Glass');
%% 4-5)b
clc;
IF = fft2(I1);
IF = fftshift(I1);

figure(6);
mesh(10*log10(abs(IF)));
grid on; 
title('DFT Mag');

% subplot(313)
figure(7);
mesh(angle(IF));
grid on; 
title('DFT phase');

%% 4-5)c
clc;
%% 4-5)d
clc;
H = FFT_LP_2D(I1,0.1*pi);
figure(8);
subplot(121)
imshow(I1);
grid on; 
title('Glass');

subplot(122)
imshow(H);
grid on; 
title('Glass Low pass Filtered');

%% 4-5 ) h
B = imresize(I,1/4,"nearest"); %Gray Scale and 2-D are DownSampled
figure(9)
subplot(131);
imshow(I1);
grid on; 
title('Glass');

subplot(132)
imshow(B);
grid on; 
title('DownSampled Glass with coeff = 1/4');

B2 = FFT_LP_2D(B,0.11*pi);

subplot(133)
imshow(B2);
grid on; 
title('DownSampled 1/4 Glass LowPassed 0.5pi');

%% Function
function Output_image = FFT_LP_2D(input_image, cutoff_frequency)
% fc is the circular cutoff frequency which is normalized to [0 1], that is,
% the highest radian frequency \pi of digital signals is mapped to 1.
 if cutoff_frequency < 0 || cutoff_frequency > pi
     display("False Frequency You Entered");
 else
     [ir,ic,iz] = size(input_image);
     hr = (ir-1)/2;
     hc = (ic-1)/2;
     [x, y] = meshgrid(-hc:hc, -hr:hr); %Mesh which Has negatives and positives

     mg = sqrt((x/hc).^2 + (y/hr).^2); %CIrcular Image
     lp = double(mg <= cutoff_frequency);%if mg < Fc => save it

     IM = fftshift(fft2(double(input_image))); % DFT of Original Image
     IP = zeros(size(IM)); % Preallocating
    for z = 1:iz
        IP(:,:,z) = IM(:,:,z) .* lp;
    end
     Output_image = abs(ifft2(ifftshift(IP),'symmetric'));
 end
end