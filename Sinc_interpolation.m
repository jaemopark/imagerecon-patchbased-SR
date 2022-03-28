function up_result = Sinc_interpolation(low_data,up_size)
% Sinc interpolation for upsampling
% Inputs: 
%       low_data: low-res data that needs to be upsampled;
%       up_size: targeted matrix size for high-res image.
%
% Outputs:
%       up_result:upsampled high-res image after sinc interpolation

up_result = zeros(up_size);
up_factor = up_size/length(low_data);

for i = 1:up_size
    for j = 1:up_size
        for m = 1:length(low_data)
            for n = 1:length(low_data)
                up_result(i,j) = up_result(i,j) + low_data(m,n)...
                                 *sinc((mod(i-1,up_factor)+1)/up_factor-1/up_factor*(up_factor/2+0.5)-(m-1-floor((i-1)/up_factor)))...
                                 *sinc((mod(j-1,up_factor)+1)/up_factor-1/up_factor*(up_factor/2+0.5)-(n-1-floor((j-1)/up_factor)));
    
            end
        end    
    end
end
end

