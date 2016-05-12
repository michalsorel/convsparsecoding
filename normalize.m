function r = normalize(a)
%NORMALIZE normalizes values of A to be within the interval 0 and 1
%
% function r = normalize(a)
%
% NaNs and Infs are eliminated before normalization 

a(~isfinite(a)) = Inf;
minv = min(a(:));
a(~isfinite(a)) = -Inf;
maxv = max(a(:));
r = (a  - minv) / (maxv - minv);  
