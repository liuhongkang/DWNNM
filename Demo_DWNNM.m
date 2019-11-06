
clear;
clc;
Original_image_dir  =    'clean_images/';
Denoising_image_dir='denoising_images/';
im_num = 20;

nSig = [5 30 15];  %

Par.nSig      =   nSig;                                 %  ����ǿ�ȣ�����ͨ�����ӵ�����ǿ�Ȳ�ͬ�����㷨Ϊ��ͨ��
Par.win =   20;                                   %   �Ǿֲ�ͼ����������ڳߴ�
Par.Constant         =  2 * sqrt(2);                              % Ȩ�ؾ���Ĺ̶�������
Par.Innerloop =   2;                                    % �ڲ�ѭ�����������½���ͼ���ƥ��
Par.ps       =   6;                            % ͼ���ߴ�
Par.step        =   5;                         % ��������
Par.Iter          =   6;                            %  �ܵĵ����������ظ�ȥ��
Par.display = true;


Par.method = ' DWNNM_ADMM';
Par.maxIter = 10;  %ADMM �ֽ��������


Par.delta     =   0.1;                                  % Parameter between each iter
%δ�㶮���������������ӵ�����
Par.lambda = 0.75;         %Ȩ�ز�����Ӧ���ݲ�ͬ�����ȼ����е���
Par.mu = 1.001;       %ADMM���²��� mu
Par.rho = 0.05;       % ADMM�ͷ�ϵ��
% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
Par.SSIM = zeros(Par.Iter, im_num, 'single');
for i = 1:im_num
    Par.imgidx = i;   %ͼ������
    Par.nSig0 = nSig;   %�����ȼ�
    Par.nlsp        =   70;   % ��ʼ���Ǿֲ�ͼ�������
    imgName=[num2str(i,'%02d'),'clean.png'];
    imgPath=[Original_image_dir,imgName];
    Par.I =  double( imread(imgPath));    
    [h, w, ch] = size(Par.I);     
    Par.nim = zeros(size(Par.I));    %����ͼ��
    Par.weight=zeros(size(Par.I));   %Ȩ�ؾ���
    for c = 1:ch
        randn('seed',0);
        Par.nim(:, :, c) = Par.I(:, :, c) + Par.nSig0(c) * randn(size(Par.I(:, :, c))); %��������
        Par.weight(:,:,c)=local_structure_weight(Par.nim(:, :, c),'Gaussian',Par.nSig0(c));
    end
    fprintf('%s denoising process start:\n',imgName);
    PSNR =   csnr( Par.nim, Par.I, 0, 0 );
    SSIM      =  cal_ssim( Par.nim, Par.I, 0, 0 );%calculate PSNR and SSIM
    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
    %
    time0 = clock;
    [im_out, Par] = DWNNM_Denoising( Par.nim, Par.I,Par.weight, Par );%DWNNM model
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );   %time cost
    %correct the result
    im_out(im_out>255)=255;
    im_out(im_out<0)=0;

    % calculate the PSNR  SSIM
    Par.PSNR(Par.Iter, Par.imgidx)  =   csnr( im_out, Par.I, 0, 0 );
    Par.SSIM(Par.Iter, Par.imgidx)      =  cal_ssim( im_out, Par.I, 0, 0 );
    imname = sprintf([Par.method '_nSig' num2str(nSig(1))  '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(mu) '_lambda' num2str(lambda) '_' imgName]);
    imwrite(im_out/255, [Denoising_image_dir,imname]);
    fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.Iter, Par.imgidx),Par.SSIM(Par.Iter, Par.imgidx)     );
end
mPSNR=mean(Par.PSNR,2);    %mean value of rows
[~, idx] = max(mPSNR);     %the index of max mean value
PSNR =Par.PSNR(idx,:);
SSIM = Par.SSIM(idx,:);  %SSIM value of max index
mSSIM=mean(SSIM,2);
fprintf('The best PSNR result is at %d iteration. \n',idx);
fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
name = sprintf([Par.method '_nSig' num2str(nSig(1)) '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(mu) '_lambda' num2str(lambda) '.mat']);
save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM');