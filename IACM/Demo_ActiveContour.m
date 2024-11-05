
clc;
clear;
close all;
acc = 0; fpr = 0; tpr = 0; sn=0; sp=0;  pr=0;  npv=0;  auc=0;  f1=0;
%% 读取图片
[filename, pathname] = uigetfile({('*.jpg;*.png;*.gif;*.tif;*.bmp')},'Select a pic');

I = imread([pathname filename]);

figure; 

subplot(2,2,1), imshow(I(:,:,1)),title('R'); 
subplot(2,2,2), imshow(I(:,:,2)),title('G'); 
subplot(2,2,3), imshow(I(:,:,3)),title('B');
subplot(2,2,4), imshow(I(:,:,1)+I(:,:,2)+I(:,:,3)),title('R+G+B/3');

figure; imshow(I(:,:,2)),title('G'); 

GrayPic=rgb2gray(I); 
figure;
imshow(GrayPic);
title('RGB图像转化为灰度图像')

thresh=graythresh(I);
Pic2=im2bw(I,thresh);  
figure;
imshow(Pic2);
title('RGB图像转化为二值化图像')

thresh=graythresh(GrayPic);
Pic2_=im2bw(GrayPic,thresh);
figure;
imshow(Pic2_);
title('灰度图像转化为二值化图像')

tic;
I = I(:,:,2);  

imsize = 1;
I = imresize(I,imsize); 
[row col]= size(I);

ginv = imcomplement(I);  
adahist = adapthisteq(ginv); 

figure;
imshow(adahist);
title('自适应直方图均衡')

se = strel('ball',10,10); 
gopen = imopen(adahist,se); 
godisk = adahist - gopen;   

medfilt = medfilt2(godisk);                    
background = imopen(medfilt,strel('disk',15)); 
I2 = medfilt - background;  

%%
I = adapthisteq(I); 
I = double(I);
I = medfilt2(I, [3, 3]);
Img = I;

name = filename(1:2);
imgmsk = imread(strcat('DRIVE\mask\',name,'test_mask.gif'));

Omask = imgmsk;

[imgMF mfbw] = MatchFilter_function(I,imgmsk);
figure; 
subplot(1,3,1), imshow(I,[]); 
subplot(1,3,2), imshow(imgmsk); 
subplot(1,3,3), imshow(imgMF);

imgMF        = normalize01(imgMF);   
mf_bw        = ~mfbw;
ot = im2uint8(mat2gray(mf_bw));    
figure; 
subplot(1,3,1), imshow(mfbw,[]); 
subplot(1,3,2), imshow(mf_bw,[]); 
subplot(1,3,3), imshow(ot);
mf_bw = bwareaopen(ot,20);         

Ivessel = FrangiFilter2D(I);         
Ivessel = Ivessel.*double(imgmsk);   
Ivessel = normalize01(Ivessel);
Iv = otsu(Ivessel,2);                
for i = 1:size(Iv,1)
    for j = 1:size(Iv,2)
        if (Iv(i,j) == 1)
            Iv(i,j) = 0;
        end
    end
end

for i = 1:size(Iv,1)
    for j = 1:size(Iv,2)
        if (Iv(i,j) == 2)
            Iv(i,j) = 1;
        end
    end
end

figure;
imshow(Iv);
title('Hessian Frangi')

figure; subplot(1,3,1), imshow(Img,[]);
Img = acode_main_retin_seg(Img, imgmsk);
subplot(1,3,2), imshow(Img,[]); 
Img = im2uint8(mat2gray(Img));
subplot(1,3,3), imshow(Img);
Img = double(Img(:, :, 1));
Img = medfilt2(Img, [3, 3]);    

figure;
imshow(Img);
title('图像增强')

Img = Img+double(I2);             
u = initialcurve(Img,'gradient');
ui = u;

figure;
imshow(u);
title('图像增强')

mu = 1; lambda1 = 0.7; lambda2 = 0.6; lambda3 = 5; lambda4 = lambda3; lambda5 = 4;
timestep = 0.001; v = 1; epsilon = 0.6; iterNum = 100;
figure;
for n = 1:iterNum
    [u,u1,u2,u3] = acwe_dwt(u, Img, timestep, mu, v, lambda1, lambda2, 1, epsilon, 1, lambda3, lambda4, lambda5, Iv, mf_bw, row, col);
    if mod(n,10)==0
        pause(0.1);
        imshow(Img, []);
        hold on; axis off,axis equal
        [c,h] = contour(u,[0 0],'r'); 
        set(gcf,'color','white');
        whitebg('white')
        iterNum = [num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
end

figure;
[c,h] = contour(u,[0 0],'k');

figure;
imshow(Img, []); hold on; axis off,axis equal
[c,h] = contour(u,[0 0],'w');
totalIterNum = [num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);

