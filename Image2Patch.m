function       X = Image2Patch( im_out, par )
% record the non-local patch set and the index of each patch in
% of seed patches in image
im_out         =  single(im_out); %%单浮点数
X          =  zeros(par.ps2, par.maxrc, 'double'); %%maxrc 为整幅图像划分图像块的个数
k    =  0;
for l = 1:par.ch  %%patch 为通道数
    for i = 1:par.ps
        for j = 1:par.ps
            k    =  k+1;
            blk  = im_out(i:end-par.ps+i,j:end-par.ps+j, l);
            X(k,:) = blk(:)';
        end
    end
end
%%将每个图像块化为一列
