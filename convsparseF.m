function [U, A, H, report] = convsparseF(G,iH,method,iters,lambda)
%CONVSPARSEF Convolution Sparse Coding
%
% [U, A, H, report] = convsparseF(G,iH,method,iters,lambda)
%
% input:
% G ... gray-scale input images MxNxL (L... number of input images)
% iH ... initial filters; cell array { scale1(M1xN1xK1) scale2(M2xN2xK2) ...}
%       default: see code below
% iters ... structure containing numbers of iterations for the main loop (iters.maxiter_main), 
%       optimization of feature maps A (iters.maxiter_A), optimization of kernels (iters.maxiter_H)
% method 
%           3 = in min_A and min_H Conjugate Gradients (similar to Zeiler 2010)
%           2 = in min_A and min_H block inverse (exactly Bristow, CVPR 2013) 
%           1 = in min_A inverse formula and in min_H block inverse (modified, more efficient, Bristow, CVPR 2013)
%           0 = in min_A inverse formula and in min_H iterative inverse formula (for multiple inputs probably corresponds to Wohlberg 2016), 
%           -1 = min_A inverse formula, min_H 3D inverse formula (proposed algorithm - variant 2)
%           -2 = min_A inverse formula, min_H consensus approximative (a bit worse than standard consensus but uses less memory)
%           -3 = min_A inverse formula, min_H consensus inverse formula (proposed algorithm - variant 3)
%           -4 = min_A inverse formula, min_H tiling (proposed algorithm - variant 1)           
%     
% lambda ... ratio <0,1> in the denosing step
%       default: lambda = 0; no denoising 
%
%The expression inverse formula above means using matrix inversion lemma (aka Woodbury
%formula). For one image,  methods 0,-1,-2,-3,-4 are equivalent.
%
% output:
% U ... inferred image conv(A,H) (if lambda > 0 then denoised image)
% A ... sparse representation for each filter (MxNxsize(iH))
% H ... learned filters; same structure as iH
% report ... structure; field E: energy in every iteration
%
%Michal Sorel & Filip Sroubek 2014-16

%% PARAMETERS
if method == -4 % tile and run basic inverse formula algorithm
    method = 0;
    G = reshape(G,[size(G,1) size(G,2)*size(G,3)]);
end
% Number of iterations in each step
if ~exist('iters','var') || isempty(iters)
    iters.maxiter_main = 50; %50; % max number of iterations in the main loop
    iters.maxiter_A = 10; % max number of iterations in the minAstep
    iters.maxiter_H = 1; % max number of iterations in the minHstep
    iters.showims = true;
end
maxiter_main = iters.maxiter_main; % max number of iterations in the main loop
maxiter_A = iters.maxiter_A; % max number of iterations in the minAstep
maxiter_H = iters.maxiter_H; % max number of iterations in the minHstep

timeit_global = zeros(1,maxiter_main);
E_global = zeros(1,maxiter_main);
timeit_A = zeros(maxiter_main,maxiter_A);
%timeit_H = zeros(1,maxiter_H);
t0 = tic;

%% MAIN PARAMETERS
%(alpha, beta, xi, maxbeta, tau, gamma)
% 1e1 1e1 1e1 1e3 1.1 1e0
% 1e1 1e2 1e2 1e4 1.1 1e1
Lp = 1; % type of sparse norm l_p
if isfield(iters,'alpha'), alpha = iters.alpha; 
else
    alpha = 1e1; % sparse weight
end
if isfield(iters,'beta'), beta = iters.beta; 
else
    beta = 1e3; %1e1; %1e2; % beta || a - v ||^2
end
if isfield(iters,'xi'), xi = iters.xi; 
else
    xi = 1e3; %1e1; %1e2; % xi || a - w ||^2
end
if isfield(iters,'maxbeta'), maxbeta = iters.maxbeta; 
else
    maxbeta = 1e1; %1e3; %1e1;
end
if isfield(iters,'tau'), tau = iters.tau; 
else
    tau = 1.1;  
end
if isfield(iters,'gamma'), gamma = iters.gamma; 
else
    gamma = 0.5e0; % data term weight
end

reltol = 1e-9; % CG reltol
%% PARAMETERS END HERE
% default method is our
if ~exist('method','var') || isempty(method)
     method = 0;
end

if ~exist('lambda','var') || isempty(lambda)
    doDenoise = false;
    lambda = 0; % when performing denoising, ratio between G and inferred image U
else
    doDenoise = true;
end

% set convolution kernels at different scales (initialization)
if exist('iH','var') && ~isempty(iH)
    H = iH(:);
    S = length(H);
    hsize = zeros(S,3);
    for s = 1:S
        hsize(s,:) = [ size(H{s},1) size(H{s},2) size(H{s},3) ];
    end
else
    %% default initialization of convolution kernels H
    S = 7; % number of scales
    K = [1 1 1 1 1 1 1]; %1*ones(S,1); % number of kernels in each scale
    H = cell(S,1);
    hsize = zeros(S,3);
    for s = 1:S
        %H{s} = randn((s-1)*2+1,(s-1)*2+1,K(s));
        H{s} = ones((s-1)*2+1,(s-1)*2+1,K(s));
        %H{s} = rand(7,7,K(s));
        hsize(s,:) = [size(H{s},1) size(H{s},2) size(H{s},3)];
    end
end

gsize = [size(G,1),size(G,2),size(G,3)];
L = gsize(3); % number of input images

% variables used in min_A
asize = [gsize(1)*gsize(2),sum(hsize(:,3)),L];
A = zeros(asize);
FA = zeros(asize);

V = zeros(asize); % Auxiliar variable
B = zeros(asize); % Lagrange multiplier
%% additive
Pr = asetupLnormPrior(Lp,alpha,beta);
V = Pr.fh(A-B,abs(A-B));

% variables used in min_H
FH = zeros(asize(1),asize(2));
i = 0;
for s = 1:S
    H{s} = projSphere(H{s}); % assure that H satisfy the same condition on its magnitude as in H minimization
    FH(:,i+[1:hsize(s,3)]) = reshape(fft2(H{s},gsize(1),gsize(2)),[asize(1), hsize(s,3)]);
    i = i+hsize(s,3);
end
if method == -1 % pak zkontrolovat, jestli opravdu musim udrzovat
    FH3 = repmat(FH,[L 1]);
    W = zeros([gsize(1),gsize(2),L,asize(2)]);
    for i = 1:asize(2)
        W(:,:,:,i) = real(ifftn(reshape(FH3(:,i),[gsize(1),gsize(2),L]))); % Auxiliar variable
    end    
    C = zeros(gsize(1),gsize(2),gsize(3),asize(2)); % Lagrange multiplier
    FG3 = repmat(reshape(fftn(G),[asize(1)*L,1]),[1 size(FH,2)]);
elseif method == -3 % consensus (W is L x copy of initial H, C = 0)    
    % FH3 = repmat(FH,[L 1]); ASI BUDU MUSET PRIDAT
    W = repmat(real(ifft2(reshape(FH,[gsize(1),gsize(2),asize(2)]))),[1 1 1 L]); % Auxiliary variable
    C = zeros(size(W)); % Lagrange multiplier    N x K x L
else
    W = real(ifft2(reshape(FH,[gsize(1),gsize(2),asize(2)]))); % Auxiliary variable
    C = zeros(gsize(1),gsize(2),asize(2)); % Lagrange multiplier    
end
FG = repmat(reshape(fft2(G),[asize(1),1,L]),[1 size(FH,2)]);
report.E = [];

for i_main = 1:maxiter_main
    disp([num2str(i_main)]);
    t_g=tic;t_Astep=tic;
    % A = min_A E(A,H)
    minAstep;
    disp(['Time of A step: ' num2str(toc(t_Astep)) ' seconds']);
    % H = min_H E(A,H)
    minHstep;
    [E,E1,l1] = calcE;
    disp(['After minH: E=',num2str(E),', data=',num2str(E1),', l1=',num2str(l1)]);
    % update G (denoising step)
    if doDenoise
        minUstep;
    end
    
    % update beta, xi
    dstr = '';
    if beta < maxbeta
        beta = beta*tau;
        dstr = [dstr,'beta=',num2str(beta),' '];
    end
    if xi < maxbeta 
        xi = xi*tau;
        dstr = [dstr,'xi=',num2str(xi)];
    end
    disp(dstr);
    E_global(i_main) = E;
    timeit_global(i_main) = toc(t0);
    disp(['Time of iteration: ' num2str(toc(t_g)) ' seconds']);
end
report.timeit_global = timeit_global;
report.E_global = E_global;
report.timeit_A = timeit_A;
A = reshape(A,[gsize(1),gsize(2),asize(2),asize(3)]);
% if no denoising is selected return in U the synthesized image
if ~doDenoise
    lambda = 1;
    minUstep;
end

if iters.showims
    figure;
    m = ceil((asize(2)+2)/2);
    subplot(m,2,1); dispIm(G(:,:,1)); title('input');
    subplot(m,2,2); dispIm(U(:,:,1)); title(['synthesized: \alpha=',num2str(alpha),', \beta=',num2str(beta),', \gamma =',num2str(gamma)]);
    for s=1:asize(2)
        subplot(m,2,2+s); dispIm(A(:,:,s,1)); colorbar; %title(num2str(hsize(s,1)));
    end
    %if maxiter_H > 0
        figure;
        for s=1:S
            subplot(S,1,s); 
            z = []; 
            for i=1:hsize(s,3) 
                z = [z,linscale(H{s}(:,:,i)),zeros(hsize(s,1),1)]; 
            end;
            dispIm(z);
        end
    %end
end
%% END of main loop

%% update A
function minAstep
% size asize(1) x asize(2) x L
FHtG = repmat(conj(FH),[1 1 L]).*FG;

% prepare HtH variable for Bristow's method
if method == 2
    FHtH = zeros(asize(2),asize(2),asize(1));
    for i = 1:asize(2)
       FHtH(i,i,:) = abs(FH(:,i)).^2; 
       for j = i+1:asize(2)
           FHtH(i,j,:) = conj(FH(:,i)).*FH(:,j);
           FHtH(j,i,:) = conj(FHtH(i,j,:));
       end
    end  
    if maxiter_A > 1 % precompute inverses 
            FHtHinv = zeros(size(FHtH));tic;
            for n=1:asize(1)
                    FHtHinv(:,:,n) = inv(FHtH(:,:,n) + beta/gamma*eye(asize(2)));
            end;
            disp(['Time of inversion: ' num2str(toc) 's']);    
    end    
end

for i=1:maxiter_A
    % size asize(1) x asize(2) x L
    b = FHtG + beta/gamma*reshape(fft2(reshape(V + B,[gsize(1),gsize(2),asize(2),asize(3)])),asize);    
    switch method
        case {-3,-2,-1,0, 1}
            if (L>1)
                rFH = repmat(FH,[1 1 L]);
                FA = gamma/beta*(b-conj(rFH).*...
                    repmat(sum(rFH.*b,2)./(beta/gamma+sum(rFH.*conj(rFH),2)),[1 asize(2) 1]));
            else
                FA = gamma/beta*(b-conj(FH).*...
                    repmat(sum(FH.*b,2)./(beta/gamma+sum(FH.*conj(FH),2)),[1 asize(2) 1]));
            end
        case 2
            b = reshape(shiftdim(b,1),[asize(2),L,asize(1)]);
            if maxiter_A > 1
                for n=1:asize(1)
                    FA(n,:,:) = FHtHinv(:,:,n)*b(:,:,n);
                end            
            else
                for n=1:asize(1)
                    FA(n,:,:) = (FHtH(:,:,n) + beta/gamma*eye(asize(2)))\b(:,:,n);
                end
            end
        case 3
            [xmin,flag,relres,iter,resvec] = mycg(@gradA,vec(b),reltol,1000,[],vec(FA));
            disp(['flag, iter:',num2str([flag iter])]);
            FA = reshape(xmin,asize);
    end
    
    A = reshape(real(ifft2(reshape(FA,[gsize(1),gsize(2),asize(2),asize(3)]))),asize);
    %% additive
    Pr = asetupLnormPrior(Lp,alpha,beta); 
    %% additive
    V = Pr.fh(A-B,abs(A-B)); % Lp norm 
    %V = A-B; % no contraint
    %V = beta*(A-B)/(alpha+beta); % L2 norm
    %V(V<0) = 0;
    % update Bregman variable (Lagrange variable)
    B = B + V - A;
    timeit_A(i_main,i) = toc(t0);
end
[E,E1,l1] = calcE;
disp(['E=',num2str(E),', data=',num2str(E1),', l1=',num2str(l1)]);
report.E = [report.E, E];
end

function r = gradA(x)
    X = reshape(x,asize);
    R = repmat(sum(repmat(FH,[1 1 L]).*X,2),[1 asize(2)]).*repmat(conj(FH),[1 1 L]) ...
        + beta/gamma*X;
    r = vec(R);
end

%% update H
function minHstep
    if maxiter_H == 0 
        return;
    end
    if method == -1 % A'G for 3d convolution
        FA3 = FG3; %=zeros(size(FG3)); % could be faster like this
        for i = 1:size(FG3,2)
            FA3(:,i) = reshape(fftn(reshape(A(:,i,:),[gsize(1) gsize(2) L])),[size(FG3,1) 1]);
        end
        FAtG = conj(FA3).*FG3;   %ve 3D verzi bych nemel scitat!
    else
        FAtG = sum(conj(FA).*FG,3);   
    end

    % prepare AtA variable for Bristow's method and our method "1"
    if (method == 1) || (method == 2)
         FAtA = zeros(asize(2),asize(2),asize(1));
         for i = 1:asize(2)
             FAtA(i,i,:) = sum(abs(FA(:,i,:)).^2,3);
             for j = i+1:asize(2)
                 FAtA(i,j,:) = sum(conj(FA(:,i,:)).*FA(:,j,:),3);
                 FAtA(j,i,:) = conj(FAtA(i,j,:));
             end
         end
    end

    for i=1:maxiter_H
        if method == -1  % 3D inverse formula
            b = FAtG;
            for ii = 1:asize(2)
                b(:,ii) = b(:,ii) + xi/gamma*reshape(fftn(W(:,:,:,ii) + C(:,:,:,ii)),[asize(1)*asize(3) 1]);
            end
            FH3 = gamma/xi*(b-conj(FA3).*...
                    repmat(sum(FA3.*b,2)./(xi/gamma+sum(FA3.*conj(FA3),2)),[1 asize(2)]));         

            TH = zeros([gsize(1) gsize(2) gsize(3) asize(2)]);
            for ii = 1:asize(2)
                TH(:,:,:,ii) = real(ifftn(reshape(FH3(:,ii),[gsize(1),gsize(2),gsize(3)])));
            end    

            % update auxiliary variable
            W = TH-C;
            % zero W outside H support
            n = 0;
            for s = 1:S
                W(hsize(s,1)+1:end,:,1,n+[1:hsize(s,3)]) = 0; % everything bellow
                W(1:hsize(s,1),hsize(s,2)+1:end,1,n+[1:hsize(s,3)]) = 0; % to the right            
                n = n+hsize(s,3);
            end
            W(:,:,2:end,:) = 0; % zero W outside of the main plane

             % project W into sphere 
            W = projSphere(W); % function is in the same file under convsparseF()            
        else    
            if (method ~= -2 && method ~= -3)
                b = FAtG + xi/gamma*reshape(fft2(W + C),asize(1:2));
            end
            switch method
                case -2  % consensus inverse formula - approximative
                    b = conj(FA).*FG + xi/gamma*repmat(reshape(fft2(W+C),asize(1:2)),[1 1 L]);
                    FH = mean(gamma/xi*(b-conj(FA).* ...
                        repmat(sum(FA.*b,2)./(xi/gamma + sum(FA.*conj(FA),2)),[1 asize(2)])),3);             
                case -3 % consensus inverse formula                  
                    b = conj(FA).*FG + xi/gamma*reshape(fft2(W+C),asize);
                    FH = gamma/xi*(b-conj(FA).* ...
                        repmat(sum(FA.*b,2)./(xi/gamma + sum(FA.*conj(FA),2)),[1 asize(2)])); % jako 0, cili -2 bez prumerovani             
                case 0
                    if (L>1)
                        % iterative inverse formula
                        FH = calcInverse(b);
                    else  % L=1 or tiling
                        FH = gamma/xi*(b-conj(FA).*...
                            repmat(sum(FA.*b,2)./(xi/gamma+sum(FA.*conj(FA),2)),[1 asize(2)])); 
                    end
                case {1, 2}
                    for n=1:asize(1)
                        FH(n,:) = (FAtA(:,:,n) + xi/gamma*eye(asize(2)))\(b(n,:)).';
                    end
                case 3
                    [xmin,flag,relres,iter,resvec] = mycg(@gradH,vec(b),reltol,1000,[],vec(FH));
                    disp(['flag, iter:',num2str([flag iter])]);
                    FH = reshape(xmin,size(FH));
            end
            if method == -3
                TH = real(ifft2(reshape(FH,[gsize(1),gsize(2),asize(2:3)])));
                % update auxiliary variable
                W = mean(TH-C,4);        
            else
                TH = real(ifft2(reshape(FH,[gsize(1),gsize(2),asize(2)])));
                W = TH - C;
            end
             % zero W outside H support, average and project to sphere
            n = 0;
            for s = 1:S
                W(hsize(s,1)+1:end,:,n+(1:hsize(s,3))) = 0;
                W(1:hsize(s,1),hsize(s,2)+1:end,n+(1:hsize(s,3))) = 0;
                n = n+hsize(s,3);
            end
            % project W into sphere 
            W = projSphere(W); % function is in the same file under convsparseF()
        end        
         if method == -3
             W = repmat(W,[1 1 1 L]);
         end
         % update Bregman variable (Lagrange multiplier)    
         C = C + W - TH;         
        %calcE;
    end

    % update H - for now we update by the auxiliary variable
    if method == -1 % 3D
        n = 0;
        for s = 1:S
            %H{s} = TH(1:hsize(s,1),1:hsize(s,2),1,n+[1:hsize(s,3)]);
            H{s} = W(1:hsize(s,1),1:hsize(s,2),1,n+[1:hsize(s,3)]);
            FH(:,n+[1:hsize(s,3)]) = reshape(fft2(H{s},gsize(1),gsize(2)),[asize(1), hsize(s,3)]);
            n = n+hsize(s,3);
        end   
    elseif method == -3 % true consensus
        n = 0;
        FH = zeros(asize(1),asize(2));
        for s = 1:S
            H{s} = W(1:hsize(s,1),1:hsize(s,2),n+(1:hsize(s,3)),1); % take 1st, all are the same
            FH(:,n+[1:hsize(s,3)]) = reshape(fft2(H{s},gsize(1),gsize(2)),[asize(1), hsize(s,3)]);
            n = n+hsize(s,3);
        end       
    else % all others
        n = 0;
        for s = 1:S
            %H{s} = TH(1:hsize(s,1),1:hsize(s,2),n+[1:hsize(s,3)]);
            H{s} =  W(1:hsize(s,1),1:hsize(s,2),n+[1:hsize(s,3)]);
            FH(:,n+(1:hsize(s,3))) = reshape(fft2(H{s},gsize(1),gsize(2)),[asize(1), hsize(s,3)]);
            n = n+hsize(s,3);
        end
    end

end % minHstep

function W = projSphere(W) % projection into sphere in 2D (circle) - must be changed to work for other dims                       
            nW = squeeze(sum(sum(W.^2,1),2));            
            d = hsize(1,1)*hsize(1,2); %d = ones(asize(2),1)*hsize(1,1)*hsize(1,2);
            %            n = 0;
            %         for s = 1:S    %  d(n+[1:hsize(s,3)]) = hsize(s,1)*hsize(s,2); 
            %             d(n+[1:hsize(s,3)]) = (hsize(1,1)*hsize(1,2))^2/(hsize(s,1)*hsize(s,2)); 
            %             n = n+hsize(s,3);
            %         end
            mask = nW>d;
            if sum(mask) ~= 0 
                %W(:,:,mask) = W(:,:,mask)./repmat(reshape(sqrt(nW(mask)./d(mask)),1,1,[]),[gsize(1),gsize(2)]); 
                W(:,:,mask) = W(:,:,mask)./repmat(reshape(sqrt(nW(mask)./d),1,1,[]),[size(W,1),size(W,2)]); 
            end          
end

function r = gradH(x)
        X = reshape(x,size(FH));
        R = sum(repmat(sum(FA.*repmat(X,[1 1 L]),2),[1 asize(2)]).*conj(FA),3) ...
            + xi/gamma*X;
        r = vec(R);
end

function FH = calcInverse(b)
iA = zeros(asize(1),asize(2),asize(2));
% inv(I + H_1^T H_1)
ia = (gamma/xi)^2./(1 + gamma/xi*sum(abs(FA(:,:,1)).^2,2));
T1 = zeros(size(iA)); T2 = zeros(size(iA));
for k=1:asize(2)
    for l=1:asize(2)
        iA(:,k,l) = -FA(:,l,1).*conj(FA(:,k,1)).*ia;
        if k==l 
            iA(:,k,l) = gamma/xi + iA(:,k,l);
        end
    end
end
% inverse of the rest
for p = 2:L
    for k=1:asize(2)
        for l=1:asize(2);
            T1(:,k,l) = FA(:,k,p).*conj(FA(:,l,p));
        end
    end
    ia = 1./(1 + sum(sum(T1.*iA,2),3));
    for k=1:asize(2)
        for l=1:asize(2)
            T2(:,k,l) = T1(:,l,k).*ia;
        end
    end
    for k=1:asize(2)
        for l=1:asize(2)
            T1(:,k,l) = sum(squeeze(T2(:,k,:)).*iA(:,:,l),2);
        end
    end
    for k=1:asize(2)
        for l=1:asize(2)
            T2(:,k,l) = sum(squeeze(iA(:,k,:)).*T1(:,:,l),2);
        end
    end
    iA = iA - T2;
end
for k=1:asize(2)
    FH(:,k) = sum(squeeze(iA(:,k,:)).*b,2);
end

end

%% update G
function minUstep
if lambda == 0
    U = G;
    return;
end
U = real(ifft2(reshape(sum(FA.*repmat(FH,[1 1 asize(3)]),2),gsize(1),gsize(2),asize(3))));  
U = (1-lambda)*G + lambda*U;
FG = repmat(reshape(fft2(U),[asize(1),1,L]),[1 size(FH,2)]);
end

        
function [E,E1,l1] = calcE
   U = real(ifft2(reshape(sum(FA.*repmat(FH,[1 1 asize(3)]),2),gsize(1),gsize(2),asize(3))));
   l1 = sum(abs(A(:))); 
   E1 = sum((U(:)-G(:)).^2);
   E = gamma*E1 + alpha*l1;
   %disp(['E=',num2str(E),', data=',num2str(E1),', l1=',num2str(l1)]);
end

function [E,E1,l1] = calcEcons
   U = real(ifft2(reshape(sum(FA.*repmat(FH,[1 1 asize(3)]),2),gsize(1),gsize(2),asize(3))));
   l1 = sum(abs(A(:))); 
   E1 = sum((U(:)-G(:)).^2);
   E = gamma*E1 + alpha*l1;
   %disp(['E=',num2str(E),', data=',num2str(E1),', l1=',num2str(l1)]);
end

function E1 = calcRes % residual
   U = real(ifft2(reshape(sum(FA.*repmat(FH,[1 1 asize(3)]),2),gsize(1),gsize(2),asize(3))));
   E1 = sum((U(:)-G(:)).^2);
end


end

