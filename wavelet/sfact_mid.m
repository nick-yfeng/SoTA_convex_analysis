function [b,r] = sfact(h)
% [b,r] = sfact(h)
% mid-phase spectral factorization of a polynomial h.
% this program assumes that h does not have roots on the unit circle.
% for h with zeros on the unit circle, they should be factored out
% before using this program.
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


if max(b)+min(b) < 0
   b = -b;
end


% ------------------------------------------------------------


function r = seprts(p)
% r = seprts(p)
% This program is for spectral factorization.
% This program assumes that p does NOT have roots on the unit circle.
% Roots with high multiplicity will cause problems,
% they should be handled by extracting them prior to
% using this program.

% Ivan Selesnick, Polytechnic University, Brooklyn, NY
% selesi@poly.edu

rts = roots(p);  

% The roots INSIDE the unit circle
irts = rts(abs(rts)<1);    

k = imag(irts)==0;
Crts = irts(~k);		% complex roots
Rrts = sort(irts(k));		% real roots
Crts = cplxpair(Crts);
Crts = Crts(end:-1:1);

Prts = Rrts(Rrts > 0);		% positive real roots
Nrts = Rrts(Rrts < 0);		% negative real roots

A = [Crts(1:2:end), Crts(2:2:end)];
A1 = A(1:2:end,:);
A2 = A(2:2:end,:);

Prts = [Prts(1:2:end); 1./Prts(2:2:end)];
if rem(length(Prts),2) == 0
	if length(A2) > 0
		Crts = [A1(:); 1./A2(:)];
	else
		Crts = [A1(:)];
	end
else
	Crts = [A2(:); 1./A1(:)];
end

r = [Prts(:); Crts];

if rem(length(Prts)+length(Crts)/2,2) == 1
	Nrts = [Nrts(2:2:end); 1./Nrts(1:2:end)];
else
	Nrts = [Nrts(1:2:end); 1./Nrts(2:2:end)];
end

r = [r; Nrts(:)];



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




