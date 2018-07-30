function y = hard(x, T)
% y = hard(x, T)
%
% HARD THRESHOLDING
% for real or complex data.
%
% INPUT
%   x : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%
% OUTPUT
%   y : output of hard thresholding
%
% If x and T are both multidimensional, then they must be of the same size.

y = x;

y( abs(x) < T ) = 0;

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu
