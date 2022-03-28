function output = down_sampling(data,target_size)
% Downsmple the high-res image to low-res image via mean calculation
% Inputs:
%       data: high-res data that needs to be downsampled;
%       target_size: targeted matrix size for low-res image.
%  
% Outputs:
%       output: downsampled low-res image

output = zeros(target_size);
factor = length(data)/target_size;
for i = 1:length(output)
    for j = 1:length(output)
        output(i,j) = mean(mean(data(factor*(i-1)+1:factor*i,factor*(j-1)+1:factor*j)));
    end
end
end

