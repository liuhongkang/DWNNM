function [struct_weight] = local_structure_weight(img, filter_type, sigma)
%计算�?部结构先验，进�?�得出与图像尺寸相同的权重矩�?
%filter_type 表示滤波种类，pat_wid表示模糊核尺�?
    [img_hei, img_wid] = size(img);
    pat_wid=3;%% 卷积核尺�?
    G = fspecial('gaussian', [pat_wid pat_wid], sigma); % Gaussian kernel
    u = imfilter(img, G, 'symmetric');
    % u = img;
    [Gx, Gy] = gradient(u);
    J_rho = zeros(img_hei, img_wid, 2, 2);    
    
    if isequal(filter_type, 'Gaussian')
        K = fspecial('gaussian', [pat_wid pat_wid], 2); % Gaussian kernel
        J_rho(:,:,1,1) = imfilter(Gx.^2, K, 'symmetric'); 
        J_rho(:,:,1,2) =  imfilter(Gx.*Gy, K, 'symmetric');
        J_rho(:,:,2,1) = J_rho(:,:,1,2);
        J_rho(:,:,2,2) = imfilter(Gy.*Gy, K, 'symmetric');
    elseif isequal(filter_type, 'Wiener')
        J_rho(:,:,1,1) = wiener2(Gx.^2, [pat_wid pat_wid]);
        J_rho(:,:,1,2) =  wiener2(Gx.*Gy, [pat_wid pat_wid]);
        J_rho(:,:,2,1) = J_rho(:,:,1,2);
        J_rho(:,:,2,2) = wiener2(Gy.*Gy, [pat_wid pat_wid]);     
    end
    
    sqrt_delta = sqrt((J_rho(:,:,1,1) - J_rho(:,:,2,2)).^2 + 4*J_rho(:,:,1,2).^2);
    lambda1 = 0.5*(J_rho(:,:,1,1) + J_rho(:,:,2,2) + sqrt_delta);
    lambda2 = 0.5*(J_rho(:,:,1,1) + J_rho(:,:,2,2) - sqrt_delta);
    struct_weight = (mat2gray(lambda1 - lambda2)).^3;%%此处�?要调整，得到归一化的权重