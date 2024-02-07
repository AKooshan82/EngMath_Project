%% Constructing the MRI Image
clc,close all , clear all
load('rawkneedata.mat')
subplot(1,2,1)
imshow(dat);
title('Frequencies in the K-Space')
inverse = fftshift(abs(ifft2(dat)));
inverse = mat2gray(inverse);
subplot(1,2,2)
imshow(inverse)
imwrite(inverse , 'reconstructedMRI.jpg')
title('Final MRI Image')
%% 4.1 Histograms
clc , close all , clear all
load('rawkneedata.mat')
inverse = fftshift(abs(ifft2(dat)));
inverse = mat2gray(inverse);
kneeMRI = mat2gray(imread('kneeMRI.jpg'));
subplot(1,2,1)
imhist(inverse)
title('histogram of the noisy picture');
subplot(1,2,2)
imhist(kneeMRI)
title('histogram of the noiseless Image')
%% 4.3
subplot(1,2,1)
imshow(inverse);
title('MRI image with noise')
subplot(1,2,2)
imshow(avgFilter(inverse))
title('MRI image filtered with a Mean filter')
[peaksnr1, snr1] = calculateSNRPSNR(avgFilter(inverse), kneeMRI);
fprintf('\n The Peak-SNR value is %0.4f', peaksnr1);
fprintf('\n The SNR value is %0.4f \n', snr1);
%% 4.4
subplot(1,2,1)
imshow(inverse);
title('MRI image with noise')
subplot(1,2,2)
imshow(medfilt2(inverse,[3 3]))
title('MRI image filtered with a Median filter')
[peaksnr2, snr2] = calculateSNRPSNR(medfilt2(inverse,[3 3]), kneeMRI);
fprintf('\n The Peak-SNR value is %0.4f', peaksnr2);
fprintf('\n The SNR value is %0.4f \n', snr2);
%% 4.5
img = im2double(imread('reconstructedMRI.jpg'))
sigma = 1
kernel = zeros(5,5)
W = 0
for i=1:5
    for j=1:5
        sq_dist = (i-3)^2 + (j-3)^2;
        kernel(i,j) = exp(-1*(sq_dist)/(2*sigma*sigma))
        W = W + kernel(i,j)
    end
end
kernel = kernel / W;
[m,n] = size(img);
output = zeros(m,n);
Im = padarray(img,[2 2]);
for i=1:m
    for j=1:n
        temp = Im(i:i+4 , j : j+4);
        temp = double(temp);
        conv = temp.*kernel;
        output(i,j)= sum(conv(:));

    end
end
output = im2double(output)
subplot(1,2,1)
imshow(inverse)
title('Original')
subplot(1,2,2)
imshow(output)
title('Gaussian Filtered Image')
[peaksnr3, snr3] = calculateSNRPSNR(output, kneeMRI);
fprintf('\n The Peak-SNR value is %0.4f', peaksnr3);
fprintf('\n The SNR value is %0.4f \n', snr3);
%% 4.6
subplot(1,2,1)
imshow(inverse);
title('MRI image with noise')
subplot(1,2,2)
imshow(imnlmfilt(inverse))
title('MRI image filtered with a nlm filter')
[peaksnr4, snr4] = calculateSNRPSNR(imnlmfilt(inverse), kneeMRI);
fprintf('\n The Peak-SNR value is %0.4f', peaksnr4);
fprintf('\n The SNR value is %0.4f \n', snr4);
vars = ["Mean" ; "Median" ; "Gaussian" ; "Non-Local-Means"];
SNR = [snr1;snr2;snr3;snr4];
PSNR = [peaksnr1;peaksnr2;peaksnr3;peaksnr4];
evaluation = table(vars,SNR , PSNR)
function filteredImage = avgFilter(input)
    [m, n] = size(input);
    filteredImage = zeros(m, n);
    for i=2:m-1
        for j=2:n-1
            neighbors = input(i-1:i+1, j-1:j+1);
            filteredImage(i, j) = mean(neighbors,'all');
        end
    end
end
function [PSNR, SNR] = calculateSNRPSNR(image1, image2)
% Calculate the mean of the images
[row, col] = size(image2);
MSE = 0;
sum1 = 0;sum2 = 0;
    for i=1:row
        for j=1:col
            sum1 = sum1 + image2(i, j)*image2(i, j);
            sum2 = sum2 + (image1(i, j) - image2(i, j))^2;
        end
    end
    SNR = 10*log10(sum1/sum2); %snr
    for i=1:row
        for j=1:col
            MSE = (image1(i, j) - image2(i, j))^2 + MSE;
        end
    end
    maximum = max(image1, [], 'all');
    PSNR = 10 * log10(maximum^2/(MSE/(row * col))); %psnr
end
