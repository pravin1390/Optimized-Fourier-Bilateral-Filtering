clear, close all ; %clc;
f  =  double( imread('barbara512.png') );
[m,n,d]=size(f);


%% filter parameters
sigmas = 5;
sigmar = 20; %% Should be given in the range (10,150) and should be an integer

% Algorithm (proposed) specific error tolerance
eps = 1e-1; %% Value should be among {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7}


%% Proposed method
tic,
% Get optimum T, K and coefficients and plot the approximation error
filename = 'LUT.mat';
mfile = matfile(filename);
Kours = mfile.Kstar(ceil(log10(1/eps)),sigmar-9);
Tours = mfile.Tstar(ceil(log10(1/eps)),sigmar-9);

[coeffproposed ,reconproposed, errorours]=findcoeff(Tours,Kours,sigmar);

% Filtering
omegaours=(2*pi)/(2*Tours+1);
g_opt = compress(f, coeffproposed, sigmas, Kours, omegaours); 
Timeproposed=toc;

%% Direct implementation
img1 = directImplementationfbf(f,sigmas, sigmar); 

%% Plotting reconstructed kernel and error
samples = (-255: 255)';
b =  (1/(sqrt(2*pi)*sigmar)).*exp(-0.5*samples.^2/sigmar^2);
figure;
plot(b,'color','k','LineStyle','--'); axis tight; hold on; plot(reconproposed,'r'); axis tight; grid on; hold off;
legend('Target kernel','Proposed')
t1 = title('Reconstructed Approximation');
t1.Visible = 'on'; 

% Measuring error
error2prop = reshape(img1-g_opt, [d*m*n,1]);
MSE_mcbf2prop = sqrt(sum(error2prop.^2)/(d*m*n));
PSNR2prop=20*log10(255/(MSE_mcbf2prop));

% Displaying parameters
fprintf('Spatial deviation sigmas = %d and Range deviation sigmar = %d \n \n',sigmas,sigmar);
fprintf('Proposed:');
fprintf('\n');
fprintf('K = %d  \t \t \t \t',Kours);
fprintf('\n');
fprintf('T = %d  \t \t \t \t',Tours);
fprintf('\n');
fprintf('Kernel error = %e  \t \t',errorours);
fprintf('\n');
fprintf('PSNR = %f db  \t \t',PSNR2prop);
fprintf('\n');
fprintf('Time taken is %f sec \t \t',Timeproposed);
fprintf('\n');

%% Plotting images and measuring error

figure, 
subplot(1,3,1); imshow(uint8(f)), title('Input');
subplot(1,3,2); imshow(uint8(img1)), title('Brute-force');
subplot(1,3,3); imshow(uint8(g_opt)), title('Proposed');
