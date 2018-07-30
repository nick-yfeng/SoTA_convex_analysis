%% Example_LPFCSD_NIRS_2
%
%  Nick Feng,
% NYU Tandon,
% July 2018

%% Start

clear
% close all

printme_pdf = @(filename) print('-dpdf', sprintf('figures/Example_LPFCSD_NIRS_%s', filename));

printme_eps = @(filename) print('-depsc', sprintf('figures/Example_LPFCSD_NIRS_%s', filename));

%% Load data

y = load('wl1_col02.txt');

N = length(y);
n = 1:N;

figure(1)
clf

subplot(2, 1, 1)
plot(n, y);
xlim([0 N])
ax1 = axis;
axis(ax1);
title('Data')
xlabel('Time (n)')



%% Set filter parameters

fc = 0.05;       % fc : cut-off frequency (cycles/sample)
d = 2;           

alpha = 1/((4*(sin(pi*fc))^2)^d);

e = ones(N, 1);

D = spdiags([-e 2*e -e], 0:2, N-2, N);

DTD = D'*D;
eye = speye(N);



HTH = @(x) (eye/alpha + DTD)\(DTD*x);
LPF = @(x) (eye + alpha*DTD)\x;

% D2 = spdiags([-e 2*e -e], 0:2, N-2, N);
D1 = spdiags([e -e], 0:1, N-1, N);

%% Regarding lambdas
imp = zeros(N, 1);
imp( round(N/2) ) = 1;

hh = HTH(imp);
hh_norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter H


beta1 = 0.1;
beta2 = 0.04;
lam0A = beta1 * hh_norm;
lam1A = beta2 * hh_norm * sqrt(N);


rho = 1/(1/(alpha*4^d) + 1);

%% Generalized fused lasso



lam0 = lam0A;
lam1 = lam1A;


Nit = 300;
Innit = 500;


gamma = 0.8;

[x_gm, u, f, cost] = TAS_GM(y, lam0, lam1, HTH, LPF,gamma, rho, Nit, Innit);




%% Wavelet method
% Undecimated discrete wavelet transform (UDWT) - 
% (Translation-invariant) (Stationary WT, SWT)

addpath wavelet

% compute wavelet transform of y

h0 = [1 1]/sqrt(2);
h1 = [1 -1]/sqrt(2);
J = 5;

% Verify perfect reconstruction
w = udwt(y, J, h0, h1);
z = iudwt(w, J, h0, -h1);
z = z(1:N)';
err = y - z;
max(abs(err))

% figure(1)
% clf
% plotW(w, 2)
% title('Wavelet subbands')

%% Set thresholds

T0 = 0.5;
T = T0 * sqrt(2).^(0:J-1);

% figure(2)
% clf
% plotW(w, 2, T)
% title('Wavelet subbands')


%% Apply nonlinear shrinkage (nn garrote)

shrink = @(x, T) sign(x) .* max(abs(x) - T^2./abs(x), 0);

w2 = w;
for j = 1:J
    w2{j} = shrink(w2{j}, T(j));
end

% set low-pass to zero
w2{J+1} = 0*w2{J+1};

% figure(1)
% clf
% plotW(w2, 2)
% title('Wavelet subbands')


%% Invert the wavelet transform

x_w = iudwt(w2, J, h0, -h1);
x_w = x_w(1:N)';

% figure(2)
% clf
% plot(n, x_w, n, x_gm)

%% spline interpolation

% Calculate moving standard deviation and threshold
s = mov_sd(y,5);
s1 = hard(s,2.8);

figure(1)
clf
plot(y)
hold
plot(s)
plot(s1)
xlim([1150,1450])


% Segment and cubic interpolation

[y_seg, N_seg, IsMA] = segmentation(y,s1);

x_sp = zeros(N,1);

for i = 2:2:length(N_seg)
    temp = csaps(1:N_seg(i),y_seg{i},0.8,1:N_seg(i));
    x_sp((1:N_seg(i))+sum(N_seg(1:i-1))) = temp;
    y_seg{i} = y_seg{i} - temp';

end


% Reconstruction

alpha = floor(6.25/3);
beta = floor(6.25*2);

f_sp = reconstruction(y_seg, N_seg, alpha, beta);

figure(2)
clf
plot(f_sp)
hold
plot(y)


%%

GRAY = [1 1 1]*0.7;


figure(1)
clf

ylim1 = [-25 50];

K = 7;
subplot(K, 1, 1)
line(n, y, 'color', 'k')
title('Data');
ylim(ylim1)


subplot(K, 1, 2)
line(n, x_w, 'color', 'k')
title('Artifacts Extraction Wavelet thresholding');
ylim(ylim1)

subplot(K, 1, 3)
line(n, y, 'color', GRAY, 'lineWidth', 1)
line(n, y - x_w, 'color', 'k')
xlabel('Time (n)')
title('Corrected Data from Wavelet thresholding');
ylim(ylim1)

subplot(K, 1, 4)
line(n, x_sp, 'color', 'k')
title('Artifacts Extraction Spline Interpolation');
ylim(ylim1)

subplot(K, 1, 5)
line(n, y, 'color', GRAY, 'lineWidth', 1)
line(n, f_sp, 'color', 'k')
xlabel('Time (n)')
title('Corrected Data from Spline Interpolation');
ylim(ylim1)

subplot(K, 1, 6)
line(n, x_gm, 'color', 'k')
title('Artifacts Extraction SoTA Generalized fused lasso');
ylim(ylim1)

subplot(K, 1, 7)
line(n, y, 'color', GRAY, 'lineWidth', 1)
line(n, y - x_gm, 'color', 'k')
xlabel('Time (n)')
title('Corrected Data from SoTA Generalized fused lasso');
ylim(ylim1)

xlim1 = [1150 1450];
for i = 1:K
    subplot(K, 1, i)
    xlim(xlim1)
end

orient portrait
set(gcf, 'PaperPosition', [1.1 1.1 4.5 9] )

printme_eps('summary_wav')


