function P = asetupLnormPrior(q,alpha,beta)



switch q
    case 1
        v_star = 0;
        u_star = alpha/beta; 
    case 0
        v_star = sqrt(2*alpha/beta);
        u_star = sqrt(2*alpha/beta);
    otherwise % for 0<q<1
        leftmarker = fzero( @(v) -v+alpha/beta*v^(q-1)*(1-q)*q, [eps 10]);
        v_star = fzero( @(v) -0.5*v^2+alpha/beta*v^q*(1-q), [leftmarker 10]);
        u_star = v_star + alpha/beta*q*v_star^(q-1);
        
end
P.fh = @(x,nx) aLn(x,nx,u_star,u_star-v_star);


        