function [w] = Cal_weight_H(loc,K,p,N,data)
% Calculate weights for super-res algorithm
% Inputs:
%       loc: [x,y], location of the point that the weight is calculated for;
%       K: number of segments;
%       p: possibility maps;
%       N: size of the patch;
%       data: processed 1H image used for weight calculation.
%
% Outputs:
%       w: weight for location [x,y]

temp_data = zeros(N-1+length(data));
temp_data(1+(N-1)/2:length(temp_data)-(N-1)/2,1+(N-1)/2:length(temp_data)-(N-1)/2) = data;

N_c = reshape(data(loc(1)-(N-1)/2:loc(1)+(N-1)/2,loc(2)-(N-1)/2:loc(2)+(N-1)/2),[N*N,1]);
N_c_neighbor = N_c([1:(N*N+1)/2-1,(N*N+1)/2+1:N*N]);

w = zeros(N,N);
for i = 1:N
    for j = 1:N
        N_n = reshape(temp_data(loc(1)-(N-1)/2+i-1:loc(1)+(N-1)/2+i-1,loc(2)-(N-1)/2+j-1:loc(2)+(N-1)/2+j-1),[N*N,1]);
        N_n_neighbor = N_n([1:(N*N+1)/2-1,(N*N+1)/2+1:N*N]);
        if i ~= (N+1)/2 && j ~= (N+1)/2
            for k = 1:K
                if std(N_c_neighbor) == 0
                    w(i,j) = 0;
                else
                    w(i,j) =w(i,j) + p(k,loc(1),loc(2))...
                        *p(k,loc(1)-(N-1)/2+i-1,loc(2)-(N-1)/2+j-1)...
                        *exp(-(sqrt(sum((N_c_neighbor-N_n_neighbor).^2))/std(N_c_neighbor)).^2/2/(N*N-1))...
                        /K;
                end
            end
        end
    end
end

if sum(sum(w)) ~= 0
    w = w/sum(sum(w));
end