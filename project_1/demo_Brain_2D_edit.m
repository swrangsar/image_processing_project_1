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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(data); 	% image Size
DN = size(data); 	% Fourier data Size
TVWeight = 0.0; 	% Weight for TV penalty
Itnlim = 8;		    % Number of iterations

%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

% scale data
im_dc = ifftshift(ifft2(ifftshift(data.*mask./pdf)));
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

res = im_dc;  %Initial degraded image supplied to fnlcg function

% do iterations
tic
for n=1:5
	res = fnlCg(res,param);  %contains all the function etc for finding optimal solution given the constraints
	im_res = XFM'*res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc





