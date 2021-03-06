close all; clear all;

load brain512 

sampler=mask./pdf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(data); 	% image Size
DN = size(data); 	% Fourier data Size
param.TVWeight = 0.77; 	% Weight for TV penalty
param.FOVWeight = 1;

% scale data
im_dc = ifftshift(ifft2(ifftshift(data.*sampler))); % matrix E has been defined here
data = data/max(abs(im_dc(:)));

im_dc = im_dc/max(abs(im_dc(:)));

res = im_dc;  %Initial degraded image supplied to fnlcg function
figure(300), imshow(abs(res), []);

% do iterations
tic
for n=1:5
	res = fnlCg(res,sampler,data, param);  %initialize fnlcg
	im_res = res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc