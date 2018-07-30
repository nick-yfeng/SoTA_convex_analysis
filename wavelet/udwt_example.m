%% EXAMPLE for "udwt"

% TEST PERFECT RECONSTRUCTION PROPERTY

K = 3;
[h0, h1, g0, g1] = daubf(K);    % Get filters

K = 1;
[h0, h1, g0, g1] = daubf(K);    % Get filters

% Haar:
h0 = [1 1]/sqrt(2);
h1 = [1 -1]/sqrt(2);
g0 = h0;
g1 = -h1;

N = 300;
x = rand(1,N);
J = 4;
w = udwt(x, J, h0, h1);
y = iudwt(w, J, g0, g1);
y = y(1:N);
err = x - y;
max(abs(err))

