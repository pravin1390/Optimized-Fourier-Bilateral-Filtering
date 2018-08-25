function  hatx  =  directImplementationfbf (g, sigmas, sigmar)
w  = round(3*sigmas);
if (mod(w,2) == 0)
    w  = w+1;
end
c = (w-1)/2;
d=size(g,3);
for i=1:d
    f(:,:,i)=padarray(g(:,:,i),[c c],'symmetric');
end
% % Pre-compute Gaussian domain weights.
[X,Y] = meshgrid(-c:c,-c:c);
G = exp(-(X.^2+Y.^2)/(2*sigmas^2));
[m,n,d] = size(f);
for i=1:m-w+1
    for j=1:n-w+1
         I = f(i:i+w-1,j:j+w-1,:);
         H=(I(:,:,1)-f(i+c,j+c,1)).^2;
         for k=2:d
             H=H+(I(:,:,k)-f(i+c,j+c,k)).^2;
         end
         H=exp(-H./((2*sigmar^2)));
         % Calculate bilateral filter response.
         F = H.*G;
         norm_F = sum(F(:));
         for k=1:d
         hatx(i,j,k) = sum(sum(F.*I(:,:,k)))/norm_F;
         end     
    end
end
