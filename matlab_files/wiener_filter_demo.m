close all;clear all;
N=256;

f = checkerboard(64); % Create “checkerboard” image
PSF = fspecial('motion', 7, 45); % Move 7 pixels along the 45 degrees line

gb = imfilter(f, PSF, 'circular'); % Degradation ('circular' reduces border effect)
noise = imnoise(zeros(size(f)), 'gaussian', 0, 0.01); % Gassian noise
g = gb + noise; % Add noise

PSF = fspecial('motion', 7, 45);
fr1 = deconvwnr(g, PSF);

figure;imshow(abs(g),[])
%figure;imshow(abs(fr1),[])
%pause

Sn = abs(fft2(noise)) .^2; 

nA = sum(Sn(:))/prod(size(noise));
Sf = abs(fft2(f)) .^2;
fA = sum(Sf(:))/prod(size(f));
R = nA/fA;
fr2 = deconvwnr(g, PSF, R);
figure;imshow(abs(fr2),[])