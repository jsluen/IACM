
function [u bw] = MatchFilter_function(Img,mask)

[s1 s2 s3] = size(Img);
free = 2;
M = 0;
count = 0;
alpha = 1;
BETHA = 1;
mu = 1;
Beta = 0.5;
Sigma = 1;    
Length = 8;   
Size = 25; 
Bound = 9;  

NF = 36;
middle = round(Size/2);
F = zeros(Size,Size,NF);

GaussFilter = fspecial('Gaussian', 13, 1);

for A = 1:NF
    Ang = (A)*(pi/NF);   
    for x = -fix(Size/2):fix(Size/2)
        for y = -fix(Size/2):fix(Size/2)
            u = ((x)*cos(Ang)+(y)*sin(Ang));
            v = ((y)*cos(Ang)-(x)*sin(Ang));
            F(x+middle,y+middle,A) = 0;
            if (u>=-Bound && u<=Bound)&&(v>-Length/2 && v<Length/2)
                count = count+1;
                F(x+middle,y+middle,A) = -evpdf(u,mu,Beta); 
                M = M+F(x+middle,y+middle,A);
            end
        end
    end
    m = M/count;
    
    for x = -fix(Size/2):fix(Size/2)
        for y = -fix(Size/2):fix(Size/2)
            u = ((x)*cos(Ang)+(y)*sin(Ang));
            v = ((y)*cos(Ang)-(x)*sin(Ang));
            if (u>=-Bound && u<=Bound)&&(v>-Length/2 && v<Length/2)
                F(x+middle,y+middle,A)=(10*(F(x+middle,y+middle,A)-m));
            end
        end
    end
end

I = Img;
for i=1:NF
    Filtered_image(:,:,i)=(conv2(I,F(:,:,i),'same'));
end
Filtered_image_Reshaped=zeros(NF,s1*s2);
A1=zeros(1,s1*s2);

for i=1:NF
    Filtered_image_Reshaped(i,:)=reshape(Filtered_image(:,:,i),1,s1*s2);
end
A1=max(Filtered_image_Reshaped);
Max=max(A1);
Min=min(A1);
IG=reshape(A1,s1,s2);

IG(:,:)=((IG(:,:)-Min)/(Max-Min))* 255;
for i=1:s1
    for j=1:s2
        IG(i,j)=(IG(i,j)+2*IG(i,j)*log(IG(i,j)));
    end
end

sa = 1.1;
rt = mim(IG,sa);

Max=max(max(rt));
Min=min(min(rt));
rt(:,:)=round(((rt(:,:)-Min)/(Max-Min))* 255);
rt=otsu(rt,255);

[tt1,e1,cmtx] = myThreshold(rt);

ms = 27; 
mk = msk(rt,ms);

rt2 = 255*ones(s1,s2);
for i=1:s1
    for j=1:s2
        if rt(i,j)>=tt1 && mk(i,j)==255
            rt2(i,j)=0;
        end
    end
end
J = im2bw(rt2);

J= ~J;
[Label,Num] = bwlabel(J);
Lmtx = zeros(Num+1,1);
for i=1:s1
    for j=1:s2
        Lmtx(double(Label(i,j))+1) = Lmtx(double(Label(i,j))+1) + 1;
    end
end
sLmtx = sort(Lmtx);
cp = 0; 
for i=1:s1
    for j=1:s2
        if (Lmtx(double(Label(i,j)+1)) > cp) && (Lmtx(double(Label(i,j)+1)) ~= sLmtx(Num+1,1))
            J(i,j) = 0;
        else
            J(i,j) = 1;
        end
    end
end
for i=1:s1
    for j=1:s2
        if mk(i,j)==0
            J(i,j)=1;
        end
    end
end

Mres=imresize(mask,.995);  
dim1=round((size(mask,1)-size(Mres,1))/2);
dim2=round((size(mask,2)-size(Mres,2))/2);

for i=1:size(Mres,1)
    for j=1:size(Mres,2)
        mask(i+dim1,j+dim2)=Mres(i,j);
    end
end

for i=1:s1
    for j=1:s2
        if mask(i,j)==0
            im3(i,j)=1;
        end
    end
end

u = im2uint8(mat2gray(IG));
bw = im3;

end



