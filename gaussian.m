clc , close all, clear all;
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
imshow(output)
title('Gaussian Filtered Image')