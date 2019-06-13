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
Original_image_dir  =    ' infrared_image/';
fpath = fullfile(Original_image_dir, '*clean.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

nSig = [5 30 15];  %三个通道添加的噪声强度不同

Par.nSig      =   nSig;                                 %  噪声图像标准差
Par.win =   20;                                   % 非局部图像块搜索窗口尺寸
Par.Constant         =  2 * sqrt(2);                              % 权重矩阵的固定常数项
Par.Innerloop =   2;                                    % 内部循环次数，重新进行块匹配
Par.ps       =   6;                            % 图像块尺寸
Par.step        =   5;                         % 块匹配时的搜索步长
Par.Iter          =   6;                            %  总的迭代次数
Par.display = true;


% Par.method = 'WNNM_ADMM'
Par.method = ' DWNNM_ADMM'
Par.maxIter = 10;  %ADMM分解迭代次数

%Par.model = '2';      模式选择删除

Par.delta     =   0.1;                                  % Parameter between each iter
%未搞懂   迭代过程中削减噪声强度
Par.lambda = 0.75;         %for different noise levels, this parameter should be tuned to achieve  better performance
Par.mu = 1.001;       %更新参数 rho
Par.rho = 0.05;       %惩罚系数
% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
Par.SSIM = zeros(Par.Iter, im_num, 'single');
for i = 1:im_num
    Par.image = i;   %图像 index
    Par.nSig0 = nSig;   %噪声强度
    Par.nlsp        =   70;   % 初始化非局部图像块数量
    Par.I =  double( imread(fullfile(Original_image_dir, im_dir(i).name)) );    %读取图像
    S = regexp(im_dir(i).name, '\.', 'split');   %暂时无用
    [h, w, ch] = size(Par.I);     %图像尺寸和通道数
    Par.nim = zeros(size(Par.I));
    for c = 1:ch
        randn('seed',0);
        Par.nim(:, :, c) = Par.I(:, :, c) + Par.nSig0(c) * randn(size(Par.I(:, :, c))); %每个通道添加不同的噪声强度
    end
    fprintf('%s :\n',im_dir(i).name);
    PSNR =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM      =  cal_ssim( Par.nim, Par.I, 0, 0 );%分别计算噪声图像与 ground truth 的 PSNR 和 SSIM 值
    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
    %
    time0 = clock;
    [im_out, Par] = DWNNM_ADMM_Denoising( Par.nim, Par.I, Par );%采用唯一的去噪模式
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );   %记录去噪过程消耗时间
    %结果修正
    im_out(im_out>255)=255;
    im_out(im_out<0)=0;

    % calculate the PSNR  SSIM
    Par.PSNR(Par.Iter, Par.image)  =   csnr( im_out, Par.I, 0, 0 );
    Par.SSIM(Par.Iter, Par.image)      =  cal_ssim( im_out, Par.I, 0, 0 );
    imname = sprintf([Par.method '_nSig' num2str(nSig(1)) num2str(nSig(2)) num2str(nSig(3)) '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(mu) '_lambda' num2str(lambda) '_' im_dir(i).name]);
    imwrite(im_out/255, imname);
    fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.Iter, Par.image),Par.SSIM(Par.Iter, Par.image)     );
end
mPSNR=mean(Par.PSNR,2);    %计算所有图像 PSNR 均值，形成列向量
[~, idx] = max(mPSNR);     %寻找最大值索引
PSNR =Par.PSNR(idx,:);
SSIM = Par.SSIM(idx,:);  %两者均取所有迭代次数中的均值最大的代数
mSSIM=mean(SSIM,2);
fprintf('The best PSNR result is at %d iteration. \n',idx);
fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
name = sprintf([Par.method '_' Par.model '_nSig' num2str(nSig(1)) num2str(nSig(2)) num2str(nSig(3)) '_' Par.model '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(mu) '_lambda' num2str(lambda) '.mat']);
save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM');