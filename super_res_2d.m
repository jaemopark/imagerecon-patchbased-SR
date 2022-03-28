function [c_super_res,e_final,iteration_final] = super_res_2d(c_low,h_high,p,s_patch,target_size,...
                                   target_th_type,target_th)
% Iuputs: 
%      c_low: acquired low-res 13C image;
%      h_high: acquired high-res 1H image, same prescription as c_low;
%      p: possibility maps from high-res 1H image, # segments * # matrix size * # matrix size
%         for brain imaging, we use # segments = 3 (GM, WM and CSF));
%      s_patch: size of patch (odd number,recommended s_patch = 5); 
%      target_size: size of targeted high-res 13C image;
%      target_th_type: type of threshold that determines when the algorithm
%                      ends, either 0 or 1:
%                      0: threshold is the number of iterations needs to be run;
%                      1: threshold is the error allowed;
%      target_th: shreshold that determines when the algorithm ends
%
% Outputs:
%      c_super_res: reconstructed super-res 13C image
%      e_final: reconstructed error when the iterated reconstruction ends
%      iteration_final: number of iterations it takes to end the super-res reconstruction

tic
%% parameters
n_seg = size(p,1);
%% Preprocessing
% dowmsample the high-res 1H image if its matrix size is higher than target_size
h_low = imresize(double(h_high),target_size/size(h_high,1),'nearest');

% dowmsample the possibility map if its matrix size is higher than target_size
p_low = zeros(size(p,1),target_size,target_size);
for i = 1:size(p,1)
    p_low(i,:,:) = imresize(squeeze(p(i,:,:)),target_size/size(p,2),'nearest');
end
% Upsample low-res 13C data for initialization
[up_x,up_y] = meshgrid(1:(size(c_low,1)-1)/(target_size-1):size(c_low,1));
c_init = interp2(c_low,up_x,up_y,'linear');

%% Iterative super-res algorithm
p_1 = zeros(n_seg,s_patch-1+size(p_low,2),s_patch-1+size(p_low,3));
p_1(:,1+(s_patch-1)/2:size(p_1,2)-(s_patch-1)/2,1+(s_patch-1)/2:size(p_1,3)-(s_patch-1)/2) = p_low;

h_low_1 = zeros(s_patch-1+size(h_low,1));
h_low_1(1+(s_patch-1)/2:size(h_low_1,1)-(s_patch-1)/2,1+(s_patch-1)/2:size(h_low_1,2)-(s_patch-1)/2) = h_low;

input_1 = zeros(size(c_init)+s_patch-1);
input_1(1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2,1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2) = c_init;
input_2 = zeros(size(input_1));

% reconstruction 
for i = 1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2
    for j = 1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2
        [w] = Cal_weight_H([i,j],n_seg,p_1,s_patch,h_low_1);
        input_2(i,j) = sum(sum(w.*input_1(i-(s_patch-1)/2:i+(s_patch-1)/2,j-(s_patch-1)/2:j+(s_patch-1)/2)));
    end
end

% mean correction
input_2_temp = input_2(1+(s_patch-1)/2:size(input_2,1)-(s_patch-1)/2,1+(s_patch-1)/2:size(input_2,2)-(s_patch-1)/2);
input_2_downsample = down_sampling(input_2_temp,size(c_low,1));
error = interp2(c_low - input_2_downsample,up_x,up_y,'spline');
for i = 1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2
    for j = 1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2
        input_2(i,j) = input_2(i,j) + error(i-(s_patch-1)/2,j-(s_patch-1)/2);
    end
end

e(1) = max(max(abs(input_2 - input_1)));

iteration = 1;
if target_th_type == 0
    while iteration < target_th
         input_1 = input_2;
         for i = 1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2
            for j = 1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2
                [w] = Cal_weight_H([i,j],n_seg,p_1,s_patch,h_low_1);
                input_2(i,j) = sum(sum(w.*input_1(i-(s_patch-1)/2:i+(s_patch-1)/2,j-(s_patch-1)/2:j+(s_patch-1)/2)));
            end
         end
         
         input_2_temp = input_2(1+(s_patch-1)/2:size(input_2,1)-(s_patch-1)/2,1+(s_patch-1)/2:size(input_2,2)-(s_patch-1)/2);
         input_2_downsample = down_sampling(input_2_temp,size(c_low,1));
         error = interp2(c_low - input_2_downsample,up_x,up_y,'spline');
         for i = 1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2
            for j = 1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2
                input_2(i,j) = input_2(i,j) + error(i-(s_patch-1)/2,j-(s_patch-1)/2);
            end
         end
         iteration = iteration + 1;
         e(iteration) = max(max(abs(input_2 - input_1)));
   end
elseif target_th_type == 1
   while max(max(abs(input_2 - input_1))) > target_th
         input_1 = input_2;
         for i = 1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2
            for j = 1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2
                [w] = Cal_weight_H([i,j],n_seg,p_1,s_patch,h_low_1);
                input_2(i,j) = sum(sum(w.*input_1(i-(s_patch-1)/2:i+(s_patch-1)/2,j-(s_patch-1)/2:j+(s_patch-1)/2)));
            end
         end
         
         input_2_temp = input_2(1+(s_patch-1)/2:size(input_2,1)-(s_patch-1)/2,1+(s_patch-1)/2:size(input_2,2)-(s_patch-1)/2);
         input_2_downsample = down_sampling(input_2_temp,size(c_low,1));
         error = interp2(c_low - input_2_downsample,up_x,up_y,'spline');
         for i = 1+(s_patch-1)/2:size(input_1,1)-(s_patch-1)/2
            for j = 1+(s_patch-1)/2:size(input_1,2)-(s_patch-1)/2
                input_2(i,j) = input_2(i,j) + error(i-(s_patch-1)/2,j-(s_patch-1)/2);
            end
         end
         iteration = iteration + 1;
         e(iteration) = max(max(abs(input_2 - input_1)));
   end
end

c_super_res = input_2(1+(s_patch-1)/2:size(input_2,1)-(s_patch-1)/2,1+(s_patch-1)/2:size(input_2,2)-(s_patch-1)/2);

%% iteration information
e_final = max(max(abs(input_2 - input_1)));
iteration_final = iteration;
figure;plot([1:iteration],e);xlabel('Number of Iterations');ylabel('max(x(i+1)-x(i))');
toc
%% Comparison
c_near = interp2(c_low,up_x,up_y,'nearest');  % nearest-neighbor interpolation
c_bli = interp2(c_low,up_x,up_y,'linear');   % bilinear interpolation
c_spl = interp2(c_low,up_x,up_y,'spline');   % spline interpolation
c_sinc = Sinc_interpolation(c_low,target_size);  % sinc interpolation

figure;colormap(jet);
subplot(2,3,1);imagesc(c_low);title('low-res');axis image;
subplot(2,3,2);imagesc(c_near);title('nearest-neighbor interpolaiton');axis image;
subplot(2,3,3);imagesc(c_bli);title('bilinear interpolation');axis image;
subplot(2,3,4);imagesc(c_spl);title('spline interpolation');axis image;
subplot(2,3,5);imagesc(c_sinc);title('sinc interpolation');axis image;
subplot(2,3,6);imagesc(c_super_res);title('patch-based super-res algorithm');axis image;

end

