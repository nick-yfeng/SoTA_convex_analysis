%% EXAMPLE for "udwt"

% TEST PERFECT RECONSTRUCTION PROPERTY

for K = 1:4
    [h0, h1, g0, g1] = daubf(K);    % Get filters

    for N = 500:510
        x = rand(1,N);
        for J = 1:5
            w = udwt(x, J, h0, h1);
            y = iudwt(w, J, g0, g1);
            y = y(1:N);
            err = x - y;
            max(abs(err))
        end
    end
end
