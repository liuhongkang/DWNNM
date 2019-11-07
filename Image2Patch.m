function       X = Image2Patch( im_out, par )
% record the non-local patch set and the index of each patch in
% of seed patches in image
im_out         =  single(im_out); %%??æµ??¹æ??
X          =  zeros(par.ps2, par.maxrc, 'double'); %%maxrc ä¸ºæ?´å??¾å???????¾å??????ä¸???
k    =  0;
for i = 1:par.ps
    for j = 1:par.ps
        k    =  k+1;
        blk  = im_out(i:end-par.ps+i,j:end-par.ps+j);
        X(k,:) = blk(:)';
    end
end
