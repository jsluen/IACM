function [u,u1,u2,u3] = acwe_dwt(u0, Img,  timestep, mu, v, lambda1, lambda2, pc, ...
                    epsilon, numIter, lambda3, lambda4, lambda5, Iv, mf_bw, row, col)
 

[cA4,cH4,cV4,cD4] = dwt2(u0,'db15');             %
UU      = idwt2(cA4,[],[],[],'db15',[row col]);  % 
[cA3,cH3,cV3,cD3] = dwt2(UU,'db10');
cA     = idwt2(cA3,[],[],[],'db10',[row col]);
[cA,cH,cV,cD] = dwt2(cA,'db9');
cod_X = wcodemat(u0);     % 
cod_cA = wcodemat(cA); cod_cH = wcodemat(cH); 
cod_cV = wcodemat(cV); cod_cD = wcodemat(cD); 
cod_cH = imresize(cod_cH,[row col]);
cod_cV = imresize(cod_cV,[row col]);
cod_cD = imresize(cod_cD,[row col]);
Dch = cod_cH;
Dcv = cod_cV;
Dcd = cod_cD;

%%
u = u0;
for k1 = 1:numIter
    u = NeumannBoundCond(u);
    K = curvature_central(u);

    e = exp(1);
    Hu = 0.75*(0.5*(1 + (u./(sqrt(1+(u.^2)))) + 0.5*(1./(1+exp(-u)))));
    DrcU = 0.2.* (1./sqrt(u.^2+1)-u.^2./((u.^2+1).^(3/2)) + (e.^(-u))./((1+(e.^(-u))).^2));

    th = mean(Hu(:));
    inside_idx = find(Hu(:) < th);
    outside_idx = find(Hu(:) >= th);

    c1 = mean(Img(inside_idx));
    c2 = mean(Img(outside_idx));
    
    data_force = -DrcU.*(mu*K - v - lambda1*(Img-c1).^2 + lambda2*(Img-c2).^2 + ... % Eq.(1) or Eq.(2)?
                  lambda3*real(sqrt(Dch.*Img)) + lambda4*real(sqrt(Dcv.*Img)) + lambda5*real(sqrt(Dcd.*Img)));  % 小波项

    u = u + timestep*(data_force);   %

    alpha1 = 0.8;
    alpha2 = 0.07;
    alpha3 = 0.08;
    u = alpha1*u + alpha2*Iv.*u + alpha3*mf_bw.*u;   % Eq.(6)  Iv来着Hessian矩阵  mf_bw来自匹配滤波

    u1 = 0; u2 = 0; u3 = 0;
    
    
end                 %

function g = NeumannBoundCond(f)
[nrow, ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);

function k = curvature_central(u)
[ux, uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./normDu; Ny = uy./normDu;
[nxx, junk] = gradient(Nx); 
[junk, nyy] = gradient(Ny);
k = nxx+nyy;






