function [b,r] = sfact(h)
% [b,r] = sfact(h)
% spectral factorization of a polynomial h.
% b: new polynomial
% r: roots of new polynomial
%
% % example:
%    g = rand(1,10);
%    h = conv(g,g(10:-1:1));
%    b = sfact(h);
%    h - conv(b,b(10:-1:1)) % should be 0

% required subprograms: seprts.m, leja.m
%
% Ivan Selesnick, Polytechnic University, Brooklyn, NY
% selesi@poly.edu
%
% leja.m is by Markus Lang, and is available from the 
% Rice University DSP webpage.

if length(h) == 1
	b = sqrt(h);
	r = [];
	return
end

% Get the appropriate roots.
r = seprts(h);	

% Form the polynomial from the roots
r = leja(r);
b = poly(r);	

if isreal(h)
	b = real(b);
end

% normalize
b = b*sqrt(max(h)/sum(b.*b));

% ------------------------------------------------------------


function r = seprts(p)
% r = seprts(p)
% This program is for spectral factorization.
% The roots on the unit circle must have even degree.
% Roots with high multiplicity will cause problems,
% they should be handled by extracting them prior to
% using this program.

% Ivan Selesnick, Polytechnic University, Brooklyn, NY
% selesi@poly.edu

SN = 0.00001;	% Small Number (criterion for deciding if a 
		% root is on the unit circle).

rts = roots(p);  

% The roots INSIDE the unit circle
irts = rts(abs(rts)<(1-SN));    

% The roots ON the unit circle            
orts = rts((abs(rts)>=(1-SN)) & (abs(rts)<=(1+SN)));
N = length(orts);
if rem(N,2) == 1
        disp('Sorry, but there is a problem (1) in seprts.m')
        r = [];
	return
end

% The roots at -1
f = find(abs(orts+1)<SN);
K = length(f);
if K > 0
	orts(f) = [];
	N = N - K;    % number of roots on unit circle, (besides -1)
end
if rem(K,2) == 1
        disp('Sorry, but there is a problem (2) in seprts.m')
	r = [];
	return
end
% The roots at -1 must be separated from the other roots on
% the unit circle, otherwise the following sort by angle will 
% put half the roots at -1 first, and the other half last.

% Sort roots on the unit circle by angle
[tmp,k] = sort(angle(orts));
orts = orts(k);

% Make final list of roots
a = (angle(orts(1:2:N)) + angle(orts(2:2:N)))/2;
orts = exp(sqrt(-1)*a);
r = [irts; orts; -ones(K/2,1)];


% ------------------------------------------------------------


function [x_out] = leja(x_in)
% function [x_out] = leja(x_in)
%
%    Input:   x_in
%
%    Output:  x_out
%
%    Program orders the values x_in (supposed to be the roots of a 
%    polynomial) in this way that computing the polynomial coefficients
%    by using the m-file poly yields numerically accurate results.
%    Try, e.g., 
%               z=exp(j*(1:100)*2*pi/100);
%    and compute  
%               p1 = poly(z);
%               p2 = poly(leja(z));
%    which both should lead to the polynomial x^100-1. You will be
%    surprised!
%
%

%File Name: leja.m
%Last Modification Date: %G%	%U%
%Current Version: %M%	%I%
%File Creation Date: Mon Nov  8 09:53:56 1993
%Author: Markus Lang  <lang@dsp.rice.edu>
%
%Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 1993  Rice University
%
%Permission is granted for use and non-profit distribution providing that this
%notice be clearly maintained. The right to distribute any portion for profit
%or as part of any commercial product is specifically reserved for the author.
%
%Change History:
%

x = x_in(:).'; n = length(x);

a = x(ones(1,n+1),:);
a(1,:) = abs(a(1,:));

[dum1,ind] = max(a(1,1:n));  
if ind~=1
  dum2 = a(:,1);  a(:,1) = a(:,ind);  a(:,ind) = dum2;
end
x_out(1) = a(n,1);
a(2,2:n) = abs(a(2,2:n)-x_out(1));

for l=2:n-1
  [dum1,ind] = max(prod(a(1:l,l:n)));  ind = ind+l-1;
  if l~=ind
    dum2 = a(:,l);  a(:,l) = a(:,ind);  a(:,ind) = dum2;
  end
  x_out(l) = a(n,l);
  a(l+1,(l+1):n) = abs(a(l+1,(l+1):n)-x_out(l));
end
x_out = a(n+1,:);




