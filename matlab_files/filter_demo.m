close all;clear all;
N=256;

%For uniform noise distribution
I=0.5*ones(N);

%For salt and pepper
I1=phantom(N);

%For gaussian
I2=phantom(N);

J1 = imnoise(I1,'salt & pepper', 0.02);
J2 = imnoise(I2,'gaussian',0,0.01);

J3=I+0.1*rand(N);

m=3;
n=3;

%geometric mean filter
f = exp(imfilter(log(J2), ones(m, n), 'replicate')) .^ (1/(3*3));
%figure(100);imshow(abs(J2),[])
%figure(200);imshow(abs(f),[])
%figure(100);imagesc(abs(J2))
%figure(200);imagesc(abs(f))

%Arithmetic mean filter
f = imfilter(J2, fspecial('average', [m n]));
%figure(300);imshow(abs(f),[])
%figure(300);imagesc(abs(f))

%max and min filters
f1 = ordfilt2(J1, m*n, ones(m, n)), f2 = ordfilt2(J1, 1, ones(m, n))
figure(200);imshow(abs(f1),[])
figure(300);imshow(abs(f2),[])

pause
f=0.5*(f1+f2);
figure(400);imshow(abs(f),[])


%Adaptive filter
sn=0.00001; %Global variance of noise- assumed to be known
mean = imfilter(J2, ones(m, n), 'replicate')/(m*n); % mean: local mean
var = nlfilter(J2, [m n], 'std2(x).^2'); % var: local variance
f = J2 - (sn ./ (var+eps)) .* (J2 - mean); % sn: variance of noise
figure;imshow(abs(f),[])




