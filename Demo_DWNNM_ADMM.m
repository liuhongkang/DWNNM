%-------------------------------------------------------------------------------------------------------------
% This is an implementation of the MCWNNM algorithm for real color image denoising.
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk
%              The Hong Kong Polytechnic University
%
% Please refer to the following paper if you use this code:
%
% @article{MCWNNM,
% 	author = {Jun Xu and Lei Zhang and David Zhang and Xiangchu Feng},
% 	title = {Multi-channel Weighted Nuclear Norm Minimization for Real Color Image Denoising},
% 	journal = {ICCV},
% 	year = {2017}
% }
%
% This is not the final version of MCWNNM, please do not distribute.
% Please see the file License.txt for the license governing this code.
%-------------------------------------------------------------------------------------------------------------

clear;
Original_image_dir  =    'kodak_color/';
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

nSig = [5 30 15];

Par.nSig      =   nSig;                                 % Variance of the noise image
Par.win =   20;                                   % Non-local patch searching window
Par.Constant         =  2 * sqrt(2);                              % Constant num for the weight vector
Par.Innerloop =   2;                                    % InnerLoop Num of between re-blockmatching
Par.ps       =   6;                            % Patch size
Par.step        =   5;
Par.Iter          =   6;                            % total iter numbers
Par.display = true;


% Par.method = 'WNNM_ADMM'
Par.method = 'MCWNNM_ADMM'
Par.maxIter = 10;

Par.model = '2';

Par.delta     =   0.1;                                  % Parameter between each iter
Par.lambda = 0.75;
Par.mu = 1.001;
Par.rho = 0.05;
% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
Par.SSIM = zeros(Par.Iter, im_num, 'single');
for i = 1:im_num
    Par.image = i;
    Par.nSig0 = nSig;
    Par.nlsp        =   70;   % Initial Non-local Patch number
    Par.I =  double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
    S = regexp(im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);
    Par.nim = zeros(size(Par.I));
    for c = 1:ch
        randn('seed',0);
        Par.nim(:, :, c) = Par.I(:, :, c) + Par.nSig0(c) * randn(size(Par.I(:, :, c)));
    end
    fprintf('%s :\n',im_dir(i).name);
    PSNR =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM      =  cal_ssim( Par.nim, Par.I, 0, 0 );
    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
    %
    time0 = clock;
    if Par.model == '1'
        [im_out, Par] = MCWNNM_ADMM1_Denoising( Par.nim, Par.I, Par );
    elseif Par.model == '2'
        [im_out, Par] = MCWNNM_ADMM2_Denoising( Par.nim, Par.I, Par );
    else
        [im_out, Par] = MCWNNM_ADMM_Denoising( Par.nim, Par.I, Par );
    end
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
    im_out(im_out>255)=255;
    im_out(im_out<0)=0;
    % calculate the PSNR
    Par.PSNR(Par.Iter, Par.image)  =   csnr( im_out, Par.I, 0, 0 );
    Par.SSIM(Par.Iter, Par.image)      =  cal_ssim( im_out, Par.I, 0, 0 );
    imname = sprintf([Par.method '_nSig' num2str(nSig(1)) num2str(nSig(2)) num2str(nSig(3)) '_' Par.model '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(mu) '_lambda' num2str(lambda) '_' im_dir(i).name]);
    imwrite(im_out/255, imname);
    fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.Iter, Par.image),Par.SSIM(Par.Iter, Par.image)     );
end
mPSNR=mean(Par.PSNR,2);
[~, idx] = max(mPSNR);
PSNR =Par.PSNR(idx,:);
SSIM = Par.SSIM(idx,:);
mSSIM=mean(SSIM,2);
fprintf('The best PSNR result is at %d iteration. \n',idx);
fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
name = sprintf([Par.method '_' Par.model '_nSig' num2str(nSig(1)) num2str(nSig(2)) num2str(nSig(3)) '_' Par.model '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(mu) '_lambda' num2str(lambda) '.mat']);
save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM');