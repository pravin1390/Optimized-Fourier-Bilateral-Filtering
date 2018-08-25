function g = compress(src, coeff, sigmaS, K, omega)
% create spatial filter
w  = round(3*sigmaS); if (mod(w,2) == 0); w  = w+1; end
filt     = fspecial('gaussian', [w w], sigmaS);
numer=coeff(1).*imfilter(src, filt, 'symmetric');
denom=coeff(1).*ones(size(src));
for k = 2:K
	Iphase = omega*(k-1)*src;
	Ic = cos(Iphase);
	Is = sin(Iphase);
	CIc  = imfilter(Ic, filt, 'symmetric');
	CIs  = imfilter(Is, filt, 'symmetric');
	CIcp = imfilter(Ic.*src, filt, 'symmetric');
	CIsp = imfilter(Is.*src, filt, 'symmetric'); 
    numer = numer +coeff(k)*(Ic.*CIcp+Is.*CIsp);
 	denom = denom +coeff(k)*(Ic.*CIc +Is.*CIs );
end
g = numer./denom;
