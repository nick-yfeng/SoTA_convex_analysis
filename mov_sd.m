function s = mov_sd(x, k)

x = x(:);
N = length(x);
W = 2*k + 1;
s = zeros(size(x));
w = ones(W,1);

temp = sqrt((conv(x.^2, w) - conv(x, w).^2/W)/W); 

s(k+1:N-k) = temp(W:end-2*k); 

end