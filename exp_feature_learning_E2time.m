%im = read_out_mesto_4(true);
%save  out_mesto_4 im
% Read input images stored in a .mat file
load  out_mesto_4

j = 1; 
clear In
In = zeros(128,128,size(im,3));
for i=1:size(im,3)    
    I = double(rgb2gray(squeeze(im(:,:,i,:))))/255;
    I = I(101:228,101:228); % smaller crop 128x128    
    [In(:,:,j) k q] = localnormalize(I,4,4); % for nature sigma=5 otherwise sigma=4
    j = j+1;
end

% Set numbers of iterations
iters.maxiter_main = 100; %50; % max number of iterations in the main loop
iters.maxiter_A = 5; % max number of iterations in the minAstep
iters.maxiter_H = 5; % max number of iterations in the minHstep
iters.showims = false;
iters.beta = 1e3;
iters.xi = 1e3;

itersConsensus = iters;
itersConsensus.beta = 1e2;
itersConsensus.xi = 1e2;

nexp = 50; % experiment number
withBristow =false;
if ~withBristow, disp('For speed the Bristow''s algorithm is skipped'); end
L=10; %L = [1 10 50];
K=50; %K = [1:10 15:5:50 60:10:100]; Max. 59
M=9; % kernel size is MxM
rep1 = cell(length(L),length(K)); % reports 
rep2 = cell(size(rep1));
rep3 = cell(size(rep1));
rep4 = cell(size(rep1));
rep5 = cell(size(rep1));
kernels = randn(M,M,K(end)); 
t = zeros(length(L),length(K),4);    
for  l = 1:length(L)
    input = In(:,:,1:L(l));
    for k = 1:length(K)        
        iH{1} = kernels(:,:,1:K(k));
        if withBristow
            disp('---------------------------------------------------');
            disp('Bristow');
            disp('---------------------------------------------------');
            tic;[Un, A1, H1, rep1{l,k}] = convsparseF(input,iH,2,iters);t(l,k,1)=toc; %Bristow
        end
        disp('---------------------------------------------------');
        disp('Proposed 3D');
        disp('---------------------------------------------------');
        tic;[Un, A3, H3, rep3{l,k}] = convsparseF(input,iH,-1,iters);t(l,k,3)=toc; % 3D Woodbury
        if L == 1
            tic;[Un, A4, H4, rep4{l,k}] = convsparseF(input,iH,0,itersConsensus);t(l,k,4)=toc; % single image alg.
            %for L>1 this version uses recursion, which  is  extremely slow for more than just a few images
        else
            disp('---------------------------------------------------');
            disp('Proposed consensus');
            disp('---------------------------------------------------');
            tic;[Un, A4, H4, rep4{l,k}] = convsparseF(input,iH,-3,itersConsensus);t(l,k,4)=toc; % consensus
            disp('---------------------------------------------------');
            disp('Proposed consensus approximative');
            disp('---------------------------------------------------');
            tic;[Un, A4, H4, rep5{l,k}] = convsparseF(input,iH,-2,itersConsensus);t(l,k,4)=toc; % consensus approximative
            disp('---------------------------------------------------');
            disp('Proposed tiling');
            disp('---------------------------------------------------');
            tic;[Un, A2, H2, rep2{l,k}] = ...
               convsparseF(reshape(input,[size(In,1) size(In,2)*size(input,3)]),...
                                   iH,0,iters);t(l,k,2)=toc; % Tiling - not computed for L=1                            
        end
        disp(['L=' num2str(L(l)) ', K=' num2str(K(k)) ', time: ' num2str(toc) 's']);
    end        
    % with huge A's stored
     if withBristow
        save(['timegraph_' num2str(nexp)],'t','rep1','rep2','rep3','rep4',...
                 'H1','H2','H3','H4','A1','A2','A3','A4');
     else
         save(['timegraph_' num2str(nexp)],'t','rep1','rep2','rep3','rep4',...
         'H2','H3','H4','A2','A3','A4');
     end
end

if withBristow
    % Energy vs. time
    q = 200;
    figure(q+10); 
    %K = [1:10 15:5:50 60:10:100]; 
    nK = length(K); % max. number of kernels
    nL = length(L); % max. number of input images
    plot( rep1{nL,nK}.timeit_global,rep1{nL,nK}.E_global,'b--',...       % Bristow
            rep2{nL,nK}.timeit_global,rep2{nL,nK}.E_global,'k-.',...       % Proposed tiling
            rep3{nL,nK}.timeit_global,rep3{nL,nK}.E_global,'r:',...        % Proposed 3D convolution
            rep4{nL,nK}.timeit_global,rep4{nL,nK}.E_global,'g-',...        % Proposed consensus
            'linewidth',2);
    xlabel('Time');
    ylabel(['Energy']);
    title(['Convergence, K=' num2str(K(nK)) ', L=' num2str(L(nL))  ', P=' num2str(iters.maxiter_A)])
    legend('Bristow et al.','Proposed - tiling','Proposed - 3D','Proposed - consensus','Location','NorthEast');
end

q = 100;   % comparison of our methods (without approximative consensus)
figure(q+11); 
nK = length(K); % max. number of kernels
nL = length(L); % max. number of input images
plot( ... %rep1{nL,nK}.timeit_global,rep1{nL,nK}.E_global,'b--',...       % Bristow
        rep2{nL,nK}.timeit_global,rep2{nL,nK}.E_global,'k-.',...       % Proposed tiling
        rep3{nL,nK}.timeit_global,rep3{nL,nK}.E_global,'r:',...         % Proposed 3D convolution
        rep4{nL,nK}.timeit_global,rep4{nL,nK}.E_global,'g-',...        % Proposed consensus
        'linewidth',2);
xlabel('Time [s]');
ylabel(['Energy']);
title(['Convergence of proposed methods, K=' num2str(K(nK)) ', L=' num2str(L(nL)), ', P=' num2str(iters.maxiter_A)])
%legend('Bristow et al.','Proposed tiling','Proposed 3D','Proposed consensus','Location','NorthWest');
legend('Proposed - tiling','Proposed - 3D','Proposed - consensus','Location','NorthEast');
print( q+11, ['E2time_L10_' num2str(nexp) '.eps'], '-deps2c', '-tiff' );
%eps2pdf(['E2time_L10_' num2str(nexp) '.eps'],'c:\program files\gs\gs9.10\bin\gswin64c.exe');

q = 300; % comparison of our methods (with approximative consensus)
figure(q+11); 
%K = [1:10 15:5:50 60:10:100]; 
nK = length(K); % max. number of kernels
nL = length(L); % max. number of input images
plot( ... %rep1{nL,nK}.timeit_global,rep1{nL,nK}.E_global,'b--',...       % Bristow
        rep2{nL,nK}.timeit_global,rep2{nL,nK}.E_global,'k-.',...       % Proposed tiling
        rep3{nL,nK}.timeit_global,rep3{nL,nK}.E_global,'r:',...        % Proposed 3D convolution
        rep4{nL,nK}.timeit_global,rep4{nL,nK}.E_global,'g-',...        % Proposed consensus
        rep5{nL,nK}.timeit_global,rep5{nL,nK}.E_global,'b--',...        % Proposed consensus approximative
        'linewidth',2);
xlabel('Time [s]');
ylabel(['Energy']);
title(['Convergence of proposed methods, K=' num2str(K(nK)) ', L=' num2str(L(nL)), ', P=' num2str(iters.maxiter_A)])
%legend('Bristow et al.','Proposed tiling','Proposed 3D','Proposed consensus','Location','NorthWest');
legend('Proposed - tiling','Proposed - 3D','Proposed - consensus','Proposed - consensus approx.','Location','NorthEast');
print( q+11, ['E2time_L10_' num2str(nexp) '.eps'], '-deps2c', '-tiff' );
%eps2pdf(['E2time_L10_' num2str(nexp) '.eps'],'c:\program files\gs\gs9.10\bin\gswin64c.exe');

%figure;plot(rep2{1,1}.E_global); - the same as our algorithm

% kernel figure
% Energy distribution
[v i] = sort(squeeze(sum(sum(sum(abs(A3).^2,1),2),4)),'descend');
cs = cumsum(v'); 
figure(q+12);plot(cs./max(cs(:)),'linewidth',2);
title('Cumulative kernel energy');
xlabel('Kernel');
ylabel('Energy');
print( q+12, ['KernelEnergyCum' num2str(nexp) '.eps'], '-deps2c', '-tiff' );
%eps2pdf(['KernelEnergyCum' num2str(nexp) '.eps'],'c:\program files\gs\gs9.10\bin\gswin64c.exe');

% kernels
aux = reshape(mat2cell(squeeze(H3{1}),9,9,ones(1,50)),[size(squeeze(H3{1}),3) 1]);
figure;w=showmask(tileims(aux(i),2,10,1,-5),-3);
%figure;w=showmask(tileims(aux,2,10,1,-5),-3); % bez trideni
imwrite(w,'tiledkernels_H1.png');

%M = mat2cell(reshape(H1{1},9,[]),[9],repmat([9],[1 50]));
%figure;dispIm(cell2mat(M(i)))
