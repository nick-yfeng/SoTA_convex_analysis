%% Example_LPFCSD_toy
%
% Nick Feng,
% NYU Tandon SOE,
% July 2018

%% Start

clear
% close all

printme_pdf = @(filename) print('-dpdf', sprintf('figures/Example_LPFCSD_toy_%s', filename));

printme_eps = @(filename) print('-depsc', sprintf('figures/Example_LPFCSD_toy_%s', filename));

%% Make data
% 'toy' signal

N = 300;
n = (0:N-1)';

freq1 = 1/N;
freq2 = 5/N;
f_ = 1*sin(freq1*2*pi*n) + 0.5*sin(freq2*2*pi*n + 0.3*pi);

x_ = zeros(N, 1);
x_(40:45) = 1.5;
x_(46:48) = 3.5;
x_(78:83) = 2;
x_(84:89) = -2;


x_(190:195) = 4;
x_(196:206) = -1.5;
x_(207:212) = 5;
x_(212:222) = -1;


rng(1);                 

sigma = 0.5;

w = sigma * randn(N, 1);

data = f_ + x_ + w;

figure(1)
clf

subplot(2, 1, 1)
plot(n, f_ + x_);
xlim([0 N])
title('Noise-free test signal')
box off

subplot(2, 1, 2)
plot(n, data);
xlim([0 N])
title('Noisy test signal')
xlabel('Time (n)')
box off


%% Set filters for HTH and LPF
% Set filter parameters

% fc : cut-off frequency (cycles/sample)

fc = 2 * max([freq1 freq2]);

d = 4;           

alpha = 1/((4*(sin(pi*fc))^2)^d);

e = ones(N, 1);
% D = spdiags([e -3*e 3*e -e], 0:3, N-3, N);
D = spdiags([e -4*e 6*e -4*e e], 0:4, N-4, N);
DTD = D'*D;
eye = speye(N);

% figure(2)
% clf
% plot(f_)
% hold
% plot((eye + alpha*DTD)\f_)
% title('Verification of LPF')

HTH = @(x) (eye/alpha + DTD)\(DTD*x);
LPF = @(x) (eye + alpha*DTD)\x;

D3 = spdiags([e -3*e 3*e -e], 0:3, N-3, N);

rho = 1/(1/(alpha*4^d) + 1);


%% L1 Solution

imp = zeros(N, 1);
imp( round(N/2) ) = 1;

hh = HTH(imp);
hh_norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter H


beta1 = 0.1;
beta2 = 0.2;
lam0 = beta1 * sigma * hh_norm;
lam1 = beta2 * sigma * hh_norm * sqrt(N);

Nit = 200;


[x_L1, f_L1, cost_L1] = TAS_L1(data, lam0, lam1, HTH, LPF, rho, Nit);       % Run LPF/CSD algorithm

figure(3)
clf
plot(x_L1)
hold
plot(x_)
title(['L1 Solution, RMSE:', num2str(sqrt(sum((x_L1 - x_).^2)/N))])
legend('L1 solution','Original')


%% Generalized Fused Lasso Solution
imp = zeros(N, 1);
imp( round(N/2) ) = 1;

hh = HTH(imp);
hh_norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter H


beta1 = 0.1;
beta2 = 0.2;
lam0 = beta1 * sigma * hh_norm;
lam1 = beta2 * sigma * hh_norm * sqrt(N);

gamma = 0.6;


Nit = 200;
Innit = 200;

[x_gm, u_gm, f_gm, cost_gm] = TAS_GM(data, lam0, lam1, HTH, LPF, gamma, rho, Nit, Innit);



figure(5)
clf
plot(x_gm)
hold
plot(x_)
title(['GFL Solution, RMSE:', num2str(sqrt(sum((x_gm - x_).^2)/N))])
legend('GFL solution','Original')

%% Figure for paper : L1 and atan penalties

figure(2)
clf

ylim1 = [-3.5 5];
GRAY = [1 1 1] * 0.8;

gray_line = @(n, x) line(n, x, 'color', GRAY, 'lineWidth', 1);

n = 0:N-1;


subplot(3, 1, 1)
line(n, data,'color','k')
title('(a) Noisy observation');
xlim([0 N])
ylim(ylim1)


subplot(3, 1, 2)
gray_line(n, x_);
line(n, x_L1,'color','k');
xlim([0 N])
title(['(b) TAS with L1 penalty, RMSE:', num2str(round(sqrt(sum((x_L1 - x_).^2)/N),3))]);
ylim(ylim1)


subplot(3, 1, 3)
gray_line(n, x_)
line(n, x_gm,'color','k');
title(['(c) TAS with GF penalty, RMSE:', num2str(round(sqrt(sum((x_gm - x_).^2)/N),3))]);
xlabel('Time (n)')
xlim([0 N])
ylim(ylim1)


orient portrait

set(gcf, 'PaperPosition', [1 1 4.5 4.5] )

printme_eps('output')


