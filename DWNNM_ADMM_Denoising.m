function [E_Img, Par]   =  MCWNNM_ADMM2_Denoising( N_Img, O_Img,W_Img, Par )
E_Img           = N_Img;   % Estimated Image
[h, w, ch]  = size(E_Img);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par ); %寻找最近邻索引，包括所有关键块的自身索引，可搜索的所有图像块的数量，可搜索的所有图像块的索引序列
% noisy image to patch
NoiPat =	Image2Patch( N_Img, Par );  %划分图像块，每块包含三个通道
WeiPat=Image2Patch(W_Img,Par);  % add对于权重图像划分图像块，也对应其通道数
Par.TolN = size(NoiPat, 2);   %列数，即总的图像块数量
Sigma_arrCh = zeros(Par.ch*Par.ps2, Par.TolN);
for iter = 1 : Par.Iter
    Par.iter = iter;
    % iterative regularization
    E_Img =	E_Img + Par.delta * (N_Img - E_Img);
    % image to patch
    CurPat =	Image2Patch( E_Img, Par );
    % estimate local noise variance
    for c = 1:Par.ch
        if (iter == 1) && (Par.Iter > 1)
            TempSigma_arrCh = sqrt(max(0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
%                         TempSigma_arrCh = sqrt(abs(repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        else
            TempSigma_arrCh = Par.lambda*sqrt(max(0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
%                         TempSigma_arrCh = Par.lambda*sqrt(abs(repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        end
        Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
    end
    if (mod(iter-1, Par.Innerloop) == 0)
        Par.nlsp = Par.nlsp - 10;  % Lower Noise level, less NL patches
        NL_mat  =  Block_Matching(CurPat, Par);% Caculate Non-local similar patches for each
    end
    % Denoising by DWNNM
    [Y_hat, W_hat]  =  DWNNM_ADMM_Estimation( NL_mat, Sigma_arrCh, CurPat,WeiPat, Par );   % Estimate all the patches
    E_Img = PGs2Image(Y_hat, W_hat, Par);
    PSNR  = csnr( O_Img, E_Img, 0, 0 );
    SSIM      =  cal_ssim( O_Img, E_Img, 0, 0 );
    fprintf( 'Iter = %2.3f, PSNR = %2.2f, SSIM = %2.2f \n', iter, PSNR, SSIM );
    Par.PSNR(iter, Par.image)  =   PSNR;
    Par.SSIM(iter, Par.image)      =  SSIM;
end
return;





