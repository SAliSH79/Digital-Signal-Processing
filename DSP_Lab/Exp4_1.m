% Exp 4 _ 1
% Teacher : Dr. Abidi
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% part 4_1)a
f=imread('flower.jpeg');
figure(1);
imshow(f);
grid on; 
title('Flower');
%% part 4_1)b
f_d = im2double(f); %transfer uint8 to double
figure(2);
imshow(f_d)
grid on; 
title('Flower(Double)');
%% part 4_1)c
figure(3);
subplot(211)
imhist(f_d)
grid on; 
title('Flower Histogram Double');
axis([-0.1 1.1 -0.5 12000])

subplot(212)
imhist(f)
grid on; 
title('Flower Histogram Uint8');
axis([-25 275 -0.5 12000])

%% Again
figure(4)
subplot(311)
imhist(f(:,:,1));
grid on; 
title('Flower Histogram f(:,:,1) Uint8');
xlim([-10,260]);

subplot(312)
imhist(f(:,:,2));
grid on; 
title('Flower Histogram f(:,:,2) Uint8');
xlim([-10,260]);

subplot(313)
imhist(f(:,:,3));
grid on; 
title('Flower Histogram f(:,:,3) Uint8');
xlim([-10,260]);

%% part 4-1)d

d = histeq(f(:,:,1)) ;
figure(5)
subplot(3,2,1)
imshow(d);
grid on; 
title('Flower Histeq 1*Uint8');
     
subplot(3,2,2)
imshow(f(:,:,1));
grid on; 
title('Flower Photo 1*Uint8');

d = histeq(f(:,:,2)) ;
figure(5)
subplot(3,2,3)
imshow(d);
grid on; 
title('Flower Histeq 2*Uint8');
     
subplot(3,2,4)
imshow(f(:,:,2));
grid on; 
title('Flower Photo 2*Uint8');

d = histeq(f(:,:,3)) ;
figure(5)
subplot(3,2,5)
imshow(d);
grid on; 
title('Flower Histeq 3*Uint8');
     
subplot(3,2,6)
imshow(f(:,:,3));
grid on; 
title('Flower Photo 3*Uint8');

%% part 4-2)a
f1=imread('face.jpeg');
figure(6);
imshow(f1);
grid on; 
title('Face');
%% part 4-2)b
m = 0;
var = 0.2 ^ 2; 
J = imnoise(f1,'gaussian',m,var); % Noisy Picture
figure(7)
subplot(142)
imshow(J)
grid on; 
title('NOisy photo');

subplot(141)
imshow(f1)
grid on; 
title('Original photo');

%% Part 4-2)c , d
H1 = fspecial('average',3) ;
H2 = fspecial('average',5) ;
P1 = imfilter(J,H1);
P2 = imfilter(J,H2);

figure(7)
subplot(143)
imshow(P1)
grid on; 
title('DeNOised photo by 1/9');

subplot(144)
imshow(P2)
grid on; 
title('Denoised Photo by 1/25');

%% part 4-2 ) h
clc;
x = 0.1 ;
J2 = imnoise(f1,'salt & pepper',x); % Noisy Picture
figure(8)
subplot(142)
imshow(J2)
grid on; 
title('NOisy salt&pepper photo');

subplot(141)
imshow(f1)
grid on; 
title('Original photo');

%% part 4 - 2) v
H11 = fspecial('average',3) ;
H22 = fspecial('average',5) ;
P11 = imfilter(J2,H11);
P22 = imfilter(J2,H22);

figure(8)
subplot(143)
imshow(P11)
grid on; 
title('DeNOised photo salt & pepper by 1/9');

subplot(144)
imshow(P22)
grid on; 
title('Denoised Photo salt & pepper by 1/25');

%% part 4-2) z
clc;
load('LPF1.mat')
h1 = LPF1;
h2 = ftrans2(h1);
figure(66)
freqz(h1)
grid on; 
title('Freqz 1D');

figure(67)
freqz2(h2)
grid on; 
title('Freqz 2D');
%% part 4-2 ) w
clc;
def1 = imfilter(J2,h1);
def12 = imfilter(J2,h2);
figure(9)
subplot(221)
imshow(def1)
grid on; 
title('DeNOised photo Salt&Pepper by 1-D filter');

subplot(222)
imshow(def12)
grid on; 
title('DeNOised photo Salt&Pepper by 2-D filter');

def2 = imfilter(J,h1);
def22 = imfilter(J,h2);
figure(9)
subplot(223)
imshow(def2)
grid on; 
title('DeNOised photo Gaussian by 1-D filter');

subplot(224)
imshow(def22)
grid on; 
title('DeNOised photo Gaussian by 2-D filter');
%% part 4-2) T,K
NoisyMedIm = MedFilt(J2,3,3); % Noisy Image by Median Filter
figure(55)
subplot(131)
imshow(f1)
title("Man Face")
subplot(132)
imshow(J2)
title("Man Face with Salt&pepper")
subplot(133)
imshow(NoisyMedIm)
title("Man Face salt&pepper --> Denoised by Median")

%% 4 - 2 ) L
figure(99)
subplot(121)
imshow(NoisyMedIm)
title("Man Face salt&pepper --> Denoised by Median")

subplot(122)
imshow(P11)
title("Man Face salt&pepper --> Denoised by Mean 3x3")
%%  Functions
function denoiseIm = MedFilt(img,k1,k2)
    [m,n,l] = size(img) ;
    k11 = floor(k1/2) ;
    k22 = floor(k2/2);
    denoiseIm = img;
    for i = 1 : m
        if(mod(k1,2) == 0 ) || (mod(k2,2) == 0)
            disp('even filter size input in my Median Filter design')
            break;
        end
        for j = 1 : n
            for p = 1 : l 
                temp = img(max(i - k11,1):min(i+k11,m) , max(j - k22,1):min(j + k22,n) , p);
                temp = reshape(temp,[],1) ;
                denoiseIm(i,j,p) = median(temp);
            end
        end
    end
end