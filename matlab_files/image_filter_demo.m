close all;clear all;
N=256;

I=phantom(N);
%f=lowpassfilter([256 256],0.5,50);
%f=bandpassfilter([256 256],0.001,0.2,5);
%f=highpassfilter([256 256],0.001,50);

f=Efilter([256 256],0.1,0.3,5);

figure;imshow(abs(f),[])

F=fftshift(fft2(fftshift(I)));

Fout=F.*f;

Iout=ifftshift(ifft2(ifftshift(Fout)));

figure;imagesc(abs(Iout))

If=ifftshift(ifft2(ifftshift(f)));

%figure;imagesc((abs(If)))