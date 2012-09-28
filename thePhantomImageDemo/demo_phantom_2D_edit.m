close all; clear all;
% This script has to be developed by the students
% addpath(strcat(pwd,'/utils'));
% 
% WavePath;
% 
% if exist('FWT2_PO') <2
% 	error('must have Wavelab installed and in the path');
% end

load brain512
data = phantom('Modified Shepp-Logan', 512);
data = fftshift(fft2(fftshift(data)));
% imshow(real(data), [ ]);

% sampler=mask./pdf;
sampler = mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(data); 	% image Size
DN = size(data); 	% Fourier data Size
param.TVWeight = .1; 	% Weight for TV penalty

% scale data
im_dc = ifftshift(ifft2(ifftshift(data.*sampler))); % matrix E has been defined here
data = data/max(abs(im_dc(:)));

im_dc = im_dc/max(abs(im_dc(:)));

res = im_dc;  %Initial degraded image supplied to fnlcg function

% do iterations
tic
for n=1:5
	res = fnlCgphantom(res,sampler,data, param);  %initialize fnlcg
	im_res = res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc





