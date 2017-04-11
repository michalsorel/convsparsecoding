function [x flag relres iter resvec] = mycg(fh,b,defrelres,noiter,ph,x0)
% CG implementation
% same input and output parameters as original MATLAB pcg
flag = 1;
    nb = sqrt(b'*b);
    
    x = x0;
    r=b-fh(x);
    if isempty(ph)
        z = r;
    else
        z = ph(r);
    end;
    p=z;
    rzold=r'*z;
    for i=1:noiter
        Ap=fh(p);
        alpha=rzold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        relres = sqrt(r'*r)/nb; 
        if relres<defrelres
            flag = 0;
            break;
        end
        if isempty(ph)
            z = r;
        else
            z = ph(r);
        end
        rznew=r'*z;
        p=z+rznew/rzold*p;
        rzold=rznew;
    end
   iter = i;
   resvec = r;
end