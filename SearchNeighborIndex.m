function  par  =  SearchNeighborIndex(par)
% This Function Precompute the all the patch indexes in the Searching window
% -NeighborIndex is the array of neighbor patch indexes for each keypatch
%NeighborIndex为每一个关键块有效搜索块索引
% -NumIndex is array of the effective neighbor patch numbers for each keypatch
%NumIndex每一个关键块的有效搜索块数量
% -SelfIndex is the index of keypatches in the total patch index array
%关键块在所有图像块索引阵列中的索引
par.maxr = par.h - par.ps + 1;
par.maxc = par.w - par.ps + 1;%%划分图像块按照 step 为 1 全部划分
r          =  1:par.step:par.maxr;
par.r          =  [r r(end) + 1:par.maxr];
c          =  1:par.step:par.maxc;
par.c          =  [c c(end) + 1:par.maxc];%% 去噪过程中仅采用关键图像块，而不是遍历所有划分的图像块
par.lenr = length(par.r);
par.lenc = length(par.c);
par.ps2 = par.ps^2;  %图像块中像素个数
par.ps2ch = par.ps2 * par.ch;   %三通道图像块总的像素个数 3p^2
% Total number of patches in the test image
par.maxrc = par.maxr * par.maxc;%单通道，三通道要乘以 3
% Total number of seed patches being processed   关键图像块个数
par.lenrc = par.lenr * par.lenc;
% index of each patch in image
par.Index     =   (1:par.maxrc);
par.Index    =   reshape(par.Index, par.maxr, par.maxc);  %所有划分图像块的索引，从 1 到 maxrc
% “搜索”窗口中所有补丁索引的预设变量
par.NeighborIndex    =   int32(zeros(4 * par.win^2, par.lenrc));   
par.NumIndex        =   int32(zeros(1, par.lenrc));
par.SelfIndex   =   int32(zeros(1, par.lenrc));

for  i  =  1 : par.lenr
    for  j  =  1 : par.lenc
        row = par.r(i);
        col = par.c(j);
        off = (col-1) * par.maxr + row;   %关键块在所有图像块中的索引
        off1 = (j-1) * par.lenr + i;    %关键块索引
        
        % the range indexes of the window for searching the similar patches
        % 关键块在所有图像块中的搜索范围
        rmin    =   max( row - par.win, 1 );
        rmax    =   min( row + par.win, par.maxr );
        cmin    =   max( col - par.win, 1 );
        cmax    =   min( col + par.win, par.maxc );
        
        idx     =   par.Index(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        
        par.NumIndex(off1)  =  length(idx);
        par.NeighborIndex(1:par.NumIndex(off1), off1)  =  idx;
        par.SelfIndex(off1) = off;
    end
end
