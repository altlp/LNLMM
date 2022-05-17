% clc;
% clear;
close all;

% filename='cleandepth10.png';
% filename='barbara.tif';
% filename='ckb.jpg';

 filename='lena256.tif';
I  =  double(imread(filename));

[m,n,d]=size(I);

sigma=0.02*255; %‘Î…˘«ø∂»
S=8;K=3;
pcadim=6;
fast_flag=1;
Iact=I;
I=I+sigma*randn(m,n,d);

I2=I./255;


error2 = reshape(Iact-I, [d*m*n,1]);
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
max(I(:));
peaksnrvalue=20*log10(maxvalue)-10*log10(mse);
fprintf('\n The Peak-SNR value for noisy image is %0.4f (db) \n', peaksnrvalue);

[h1 ,w1 ,nc1] = size(I2);

radius=3;
Apca=compute_pca(I2, K, pcadim);


Cluster=20;


[m,n,d]=size(I);
pcadim=size(Apca,3);
Ares=reshape(Apca,m*n,pcadim);
[Centre,mincentre]=kmeans_recursive(Ares,Cluster);
Clusterindex=reshape(mincentre,m,n);


stestpca=sum(Apca,3)./pcadim;
figure,imshow(stestpca,[])
% Filtering
spatialtype='box';
convmethod='O1'; % 'matlab' for matlab convolutions and 'O1' for O(1) convolutions

sigma2=3.5*sigma/256;
[Ikmean2 ,WG,WG2]=fastKmeansfiltapproxinter2(I2,S,3.5*sigma/256,Centre,Clusterindex,spatialtype,convmethod,fast_flag,Apca);%denpth
Ikmean2=Ikmean2.*255;

% 
% Tkmeans=toc;
% fprintf('non local means by Kmeans complete with %d clusters \n',Cluster);
% fprintf('time for non local means (ms)=%3.0f \n',Tkmeans*1000);

%% Displaying noisy and filtered image
figure;
imshow(uint8(Iact)),%title('Original image');
figure;
imshow(uint8(I)),%title('input image');
figure;
imshow(uint8(Ikmean2));%title(['fast non local means with ',num2str(Cluster),' clusters']);

% %% PSNR calculation
[m,n,d]=size(I);
error2 = reshape(Iact-Ikmean2, [d*m*n,1]);
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
peaksnrvalue=20*log10(maxvalue)-10*log10(mse)
% ssimval=ssimcalculate(Iact,Ikmean);
% fprintf('\n mean sq error=%f \t The Peak-SNR value for filtered image is %0.4f (db) \t SSIM value is %0.4f \n',sqrt(mse), peaksnrvalue, ssimval);
