function V = aLn(DU,normDU,u_star,k)
%
% V = aLn(DU,normDU,u_star,k)
%
% half-quadratic algorithm 
% additive version
%
% For prior function phi(s) = alpha*|s-t|^q, |s|>u_star 
% phi(s) = beta/2*s^2, |s|<=u_star 
%
% Shrinkage formula:
% min_v { (u-v)^2 + lambda*|v| } 
% v_min = sign(u)*max(|u|-lambda,0);
%
% This is generalized for problems with |v|^q, 0<=q<=1

V = zeros(size(DU));

DUp = DU(normDU>u_star);
normDUp = normDU(normDU>u_star);
V(normDU>u_star) = DUp.*(normDUp-k)./normDUp; 

