function r = showmask(m,k,scale)
%SHOWMASK displays enlarge mask as image
%
%   function r = showmask(m,k,scale)
%
%   m ... mask (convolution kernel)
%   k ... 2^{-k} times zoomed
%   scale ... contrast change (impl. normalized to <0,1>)
%
%Michal Sorel (c) 2006

if ~exist('k','var'), k = -4; end
if ~exist('scale','var'), scale = 1; end

r = normalize(shrink(m,k));
imshow(r*scale);
