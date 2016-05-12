function r = shrink(a,k)
%SHRINK reduces image a to half of size (it is done k x when k specified) 
%
%   function r = shrink(a,k)
%
%   a ...      image (can be RGB), or cell aray of  images
%   k ... > 0  reduces k-times
%   k ... < 0  doubles image -k-times
%
%   r ... shrinked or enlarged image / cell array of images
%
%Remark: in case of odd size it removes odd columns/rows   
%
%History: 7. 4. 2005 'same'-> 'valid' in shrinking(l.20)
%               28.2.2006 a can be a cell array 
%           23.5.2006 made faster, works with uint variables
%
%Michal Sorel (c) 2004, 2005, 2006

if nargin == 1
    k = 1;
end

if iscell(a)
    r = cell(size(a));
    for j = 1:numel(a)
        r{j} = shraux(a{j},k);
    end
else
    r = shraux(a,k);
end

function r = shraux(a,k)
if k >= 0 %  k x zmensim na pulku
    if ndims(a) < 3
        for i = 1:k      
            %a = conv2(a,[1,1;1,1]/4,'valid');
            %a = a(1:2:end,1:2:end);
            a = (a(1:2:end-1,1:2:end-1)+a(1:2:end-1,2:2:end)...
                +a(2:2:end,1:2:end-1)+a(2:2:end,2:2:end))/4;
        end
    else
        for i = 1:k      
            %for ch = 1:size(a,3)
                %a(:,:,ch) = conv2(a(:,:,ch),[1,1;1,1]/4,'same');            
                a = (a(1:2:end-1,1:2:end-1,:)+a(1:2:end-1,2:2:end,:)...
                +a(2:2:end,1:2:end-1,:)+a(2:2:end,2:2:end,:))/4;     
            %end
            %a = a(1:2:end-1,1:2:end-1,:);
        end
    end
else % k < 0 ... enlarges image
    if ndims(a) < 3
        for i = 1:-k
        r = zeros(2*size(a));
        r(1:2:end,1:2:end) = a;
        r(2:2:end,1:2:end) = a;
        r(1:2:end,2:2:end) = a;
        r(2:2:end,2:2:end) = a;
        a = r;
        end
    else
        for i = 1:-k            
            r = zeros(2*size(a,1),2*size(a,2),size(a,3));
            r(1:2:end,1:2:end,:) = a;
            r(2:2:end,1:2:end,:) = a;
            r(1:2:end,2:2:end,:) = a;
            r(2:2:end,2:2:end,:) = a;
            a = r;            
        end
    end
end
r = a;