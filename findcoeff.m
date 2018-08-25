function [c,reconproposed,erropt]=findcoeff(T,K,sigma)
samples = (-255: 255)';
N = 2*T +1;
w0 = 2*pi/N;
b =  (1/(sqrt(2*pi)*sigma)).*exp(-0.5*samples.^2/sigma^2);
A=ones(size(samples));
for i = 2: K
    A = [A cos((i-1)*w0*samples)];
end
c=pinv(A)*b;
reconproposed=A*c;

erropt=norm((b-reconproposed),2);

