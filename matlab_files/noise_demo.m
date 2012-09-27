close all;clear all;
N=256;

%For uniform noise distribution
I=ones(N);

%For salt and pepper
I1=0.5*ones(N);

%For gaussian
I2=zeros(N);
I2(1:128,:)=0.75;
I2(129:256,:)=0.25;



J1 = imnoise(I1,'salt & pepper', 0.2);
J2 = imnoise(I2,'gaussian', 0.1);

J3=I2+0.1*rand(N);

bins=100;
% figure(100);imshow(abs(J1),[])
% figure(200);hist(J1,bins)
% 
% figure(100);imshow(abs(J2),[])
% figure(200);hist(J2,bins)

figure(100);imshow(abs(J3),[])
figure(200);hist(J3,bins)