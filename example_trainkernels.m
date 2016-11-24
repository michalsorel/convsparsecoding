% MULTISCALE EXAMPLE OF TRAINING CONVOLUTION KERNELS
%
disp('Warning: for large numbers of input images algorithm needs a lot of memory, be sure you have at least 16GB of RAM.');
disp('It also takes time, 1000 iterations can take a day or two, depending on your machine.');
load  out_mesto_4

nexp = 1; % experiment number
method = -4;

L=50;  % number of input images
%L = [1 10 50];  % multiple values can be used to compute several experiments at once
K = {[8 32 64]}; % number of kernels for each scale (for one scale simply K = {8})
% again cell-array can have more entries to compute several experiments at once

% Kernels are trained on 2D data, for this reason color images are first converted to grayscale
% Data are locally normalized, it is necessary to get beautiful kernels! 
j = 1;
clear In
In = zeros(128,128,size(im,3));
for i=1:size(im,3)    
    I = double(rgb2gray(squeeze(im(:,:,i,:))))/65535;
    I = I(101:228,101:228); % smaller crop 128x128    
    [In(:,:,j)] = localnormalize(I,4,4); % for nature sigma = 5, otherwise sigma = 4    
    j = j+1;
end  

iters.maxiter_main = 1000; % max number of iterations in the main loop
iters.maxiter_A = 10; % max number of iterations in the minAstep
iters.maxiter_H = 10; % max number of iterations in the minHstep
iters.showims = false;
iters.beta = 1e3; 
iters.xi = 1e3; 

% Consensus algorithm works better with different parameters
%itersConsensus = iters;
%itersConsensus.beta = 1e2;
%itersConsensus.xi = 1e2;

rep = cell(length(L),length(K)); % reports 
t = zeros(length(L),numel(K)); % time per experiment
for  l = 1:length(L)
    if L(l) == 1
        warning('For one image, all inverse formula based methods are equivalent');
    end
    input = In(:,:,1:L(l));
    for k = 1:numel(K)        
        iH = cell(length(K{k}),1);
        for i = 1:length(K{k})
            iH{i} = randn(2^(i+2),2^(i+2),K{k}(i)); 
        end
        tic;
        [Un, A, H, rep{l,k}] = convsparseF(input,iH,method,iters);         
        t(l,k)=toc;
        disp(['L=' num2str(L(l)) ', K=' num2str(K{k}) ', time: ' num2str(toc) 's']);
    end
end
save(['example_trainingresults_' num2str(nexp)],'t','rep','H','A'); %tiling only            

% Energy vs. time for the last experiment
figure(49); 
nK = length(K); 
nL = length(L); 
plot(rep{nL,nK}.timeit_global,rep{nL,nK}.E_global,'k-.','linewidth',2);
xlabel('Time [s]');
ylabel(['Energy']);
title(['Convergence, K=' num2str(K{nK}) ', L=' num2str(L(nL))  ', P=' num2str(iters.maxiter_A)])

% kernels
ind = 1;
for iii = 1:length(H)
    indend = ind + size(H{iii},3) - 1;
    [v i] = sort(squeeze(sum(sum(sum(abs(A(:,:,ind:indend)).^2,1),2),4)),'descend'); % only for one size - sort by A energy
    aux = reshape(mat2cell(squeeze(H{iii}),size(H{iii},1),size(H{iii},2),ones(1,size(H{iii},3))),[size(squeeze(H{iii}),3) 1]);
    figure(49+iii);w=showmask(tileims(aux(i),2,8,1,-5),-3);
    imwrite(w,['learnedkernels_L' num2str(size(H{iii},1)) '_' num2str(nexp) '.png']);
    ind = indend+1;
end
