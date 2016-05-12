function u = tileims(su,dims,nsz,bd,val)
%TILEIMS tiles images for displaying
%
%   function u = tileims(su,dims,nsz,border)
%
%   su ... cell array of images (1D or 2D)
%   dims ... 1 ... images tiled along x-axis
%        ... 2 ... tilex along both axes
%   nsz ... width of the tiling if su is 1D and dims==2
%   border ... white grid of border between images
%
%
%Michal Sorel (c) 2007

if ~exist('val','var'), val = 0.5; end

[m,n,q] = size(su{1});
[mc, nc] = size(su);
if mc==1 | nc==1 % 1D cell array input
    if dims == 1 %TADY TREBA UPRAVIT podle dims==2
        if q>1, u = zeros([m n*length(su) q]); 
        else u = zeros([m n*length(su)]); end
        for i = 1:length(su)
            u(:,1+(i-1)*n:i*n,:) = su{i};
        end
    else %dims == 2        
        mc = ceil(length(su)/nsz);nc = nsz;
        m = m+bd;n=n+bd;
        if bd>0
            if q>1, u = val*ones([m*mc+1 n*nc+1 q]); 
            else u = val*ones([m*mc+1 n*nc+1]); end    
        else
            if q>1, u = val*ones([m*mc n*nc q]); 
            else u = val*ones([m*mc n*nc]); end            
        end
        rest = mc*nc-length(su);
        for i = 1:mc
            for j = 1:nc
                u(bd+1+(i-1)*m:i*m,bd+1+(j-1)*n:j*n,:) = su{(i-1)*nsz+j};
                if i==mc & j+rest>=nc
                    break; 
                end
            end
            if i==mc, break; end;
        end           
    end
else % 2D array input (ignore nsz)
    if exist('nsz','var')
        if ~isempty(nsz), error('With 2D su don''t use parameter nsz.'); end;end;
    %u = zeros([m*mc n*nc]);
    m = m+bd;n=n+bd;
    if bd>0
        if q>1, u = val*ones([m*mc+1 n*nc+1 q]); 
        else u = val*ones([m*mc+1 n*nc+1]); end    
    else
        if q>1, u = val*ones([m*mc n*nc q]); 
        else u = val*ones([m*mc n*nc]); end            
    end
    for i = 1:mc
        for j = 1:nc
            u(bd+1+(i-1)*m:i*m,bd+1+(j-1)*n:j*n,:) = su{i,j};
        end
    end    
end

