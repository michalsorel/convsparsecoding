% MULTISCALE EXPERIMENT WITH TILING ONLY
%im = read_out_mesto_4(true);
%save  out_mesto_4 im
disp('Warning: this experiment needs a lot of memory - for machines with smaller memory can be extremely slow because of swaping');
load  out_mesto_4

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

%itersConsensus = iters;
%itersConsensus.beta = 1e2;
%itersConsensus.xi = 1e2;

nexp = 35; % experiment number
L=50; %L = [1 10 50];
K=30; %K = [1:10 15:5:50 60:10:100]; Max. 59
Levels = 3;
KK = [8 32 64];

rep1 = cell(length(L),length(K)); % reports 
rep2 = cell(size(rep1));
rep3 = cell(size(rep1));
rep4 = cell(size(rep1));
t = zeros(length(L),length(K),4); 
iH = cell(Levels,1);
for  l = 1:length(L)
    input = In(:,:,1:L(l));
    for k = 1:length(K)        
        for i = 1:Levels
            iH{i} = randn(2^(i+2),2^(i+2),KK(i)); %kernels(:,:,1:K(k));            
        end
        %tic;[Un, A1, H1, rep1{l,k}] = convsparseF(input,iH,2,iters);t(l,k,1)=toc; %Bristow
        %tic;[Un, A3, H3, rep3{l,k}] = convsparseF(input,iH,-1,iters);t(l,k,3)=toc; % 3D Woodbury
        %if L == 1
        %    tic;[Un, A4, H4, rep4{l,k}] = convsparseF(input,iH,0,itersConsensus);t(l,k,4)=toc; % cons. Woodbury
        %    %recursion is  extremely slow for more than just a few images
        %else
            %tic;[Un, A4, H4, rep4{l,k}] = convsparseF(input,iH,-2,itersConsensus);t(l,k,4)=toc; % consensus Woodb.                    
            tic;[Un, A2, H2, rep2{l,k}] = ...
                convsparseF(reshape(input,[size(In,1) size(In,2)*size(input,3)]),...
                                    iH,0,iters);t(l,k,2)=toc; % Tiling - Variant 1
        %end
        disp(['L=' num2str(L(l)) ', K=' num2str(K(k)) ', time: ' num2str(toc) 's']);
    end        
save(['timegraph_' num2str(nexp)],'t','rep2','H2','A2'); %tiling only            
end

% Energy vs. time
q = 200;
figure(q+10); 
%K = [1:10 15:5:50 60:10:100]; 
nK = length(K); % max. number of kernels
nL = length(L); % max. number of input images
plot( ...%rep1{nL,nK}.timeit_global,rep1{nL,nK}.E_global,'b--',...       % Bristow
        rep2{nL,nK}.timeit_global,rep2{nL,nK}.E_global,'k-.',...       % Proposed tiling
        ...%rep3{nL,nK}.timeit_global,rep3{nL,nK}.E_global,'r:',...        % Proposed 3D convolution
        ...%rep4{nL,nK}.timeit_global,rep4{nL,nK}.E_global,'g-',...        % Proposed consensus
        'linewidth',2);
xlabel('Time');
ylabel(['Energy']);
title(['Convergence, K=' num2str(K(nK)) ', L=' num2str(L(nL))  ', P=' num2str(iters.maxiter_A)])
%legend('Bristow et al.','Proposed - tiling','Proposed - 3D','Proposed - consensus','Location','NorthEast');
legend('Proposed - tiling','Location','NorthEast');

%figure;plot(rep2{1,1}.E_global); - the same as our algorithm

% kernel figure
% Energy distribution
[v i] = sort(squeeze(sum(sum(sum(abs(A2).^2,1),2),4)),'descend');
figure;plot(squeeze(sum(sum(abs(A2),1),2)));title('Energy distribution among kernels'); % energy distribution
cs = cumsum(v'); 
figure(q+12);plot(cs./max(cs(:)),'linewidth',2); % cumulative energy distrib.
title('Cumulative kernel energy');
xlabel('Kernel');
ylabel('Energy');
print( q+12, ['KernelEnergyCum' num2str(nexp) '.eps'], '-deps2c', '-tiff' );
%eps2pdf(['KernelEnergyCum' num2str(nexp) '.eps'],'c:\program files\gs\gs9.10\bin\gswin64c.exe');

ind = 1;
for iii = 1:length(H2)
    indend = ind + size(H2{iii},3) - 1;
    [v i] = sort(squeeze(sum(sum(sum(abs(A2(:,:,ind:indend)).^2,1),2),4)),'descend'); % only for one size
    aux = reshape(mat2cell(squeeze(H2{iii}),size(H2{iii},1),size(H2{iii},2),ones(1,size(H2{iii},3))),[size(squeeze(H2{iii}),3) 1]);
    figure(50+iii);w=showmask(tileims(aux(i),2,8,1,-5),-3);
    imwrite(w,['learnedkernels_L' num2str(size(H2{iii},1)) '_' num2str(nexp) '.png']);
    ind = indend+1;
end
