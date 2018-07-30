function [y, wx, wy, T_level, y_low] = watar(x, J, T, flag_low)
% y = watar(x, J, T)
% Wavelet trasnsient artifact reduction
%
% Input:
%   x : input signal (1D array)
%   J : number of levels
%   T : threshold
%
% Output:
%   y : estimated artifact signal
%
% Use watar(x, J, T, 'nolow') to exclude low-pass band from y.
%
% Addtional outputs:
%   [y, wx, wy, T_level, y_low] = watar(x, J, T)
%   wx : wavelet transform of x
%   wy : wavelet transform of y
%   T : threshold at each level
%   y_low : low-pass component

h0 = [1 1]/sqrt(2);
h1 = [1 -1]/sqrt(2);

N = length(x);

wx = udwt(x, J, h0, h1);

T_level = T * sqrt(2).^(0:J-1);

shrink = @(x, T) sign(x) .* max(abs(x) - T^2./abs(x), 0);

wy = wx;
for j = 1:J
    wy{j} = shrink(wy{j}, T_level(j));
end

% Remove low-pass component
if nargin > 3
    if strcmp(flag_low, 'nolow')
        wy{J+1} = 0 * wy{J+1};
    end
end

y = iudwt(wy, J, h0, -h1);
y = y(1:N);
y = reshape(y, size(x));

if nargout > 4    
    wy_low = wx;
    for j = 1:J
        wy_low{j} = 0 * wy_low{j};
    end
    y_low = iudwt(wy_low, J, h0, -h1);
    y_low = y_low(1:N);
    y_low = reshape(y_low, size(x));    
end


if 0
    
    % Verify perfect reconstruction
    z = iudwt(wx, J, h0, -h1);
    z = z(1:N);
    z = reshape(z, size(x));
    err = z - z;
    max(abs(err))
    
end
