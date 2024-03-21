using LinearAlgebra
using Printf



#
# Evalue Mubness
#
function Mubness(M)
    d,d,n=size(M);

    val = 0
    for i in 1:n-1;
        for j in i+1:n;
            M1 = M[:,:,i];
            M2 = M[:,:,j];
            val -=  2*sum(sqrt.( abs.(ones(d,d)-abs.(M1'*M2).^2)))/d;
        end
    end
    return val;

end

#
# Evalue all matrices in M are a basis
#
function IsBasis(M);
    d,d,n=size(M);

    t = 0;
    for i in 1:d;
        M1 = M[:,:,i];
        t += sum(abs.(I(d)-abs.(M1'*M1)));
    end

    return t;
end



#
# Generate one configurations of matrices
#
function GeneraOneMatConf(d,n,k,conf)
    #M = zeros(Complex,d,d,n);
    #k = size(phases,1);
    #@show(d,n,k,conf)
    t = digits(conf,base=k,pad=d^2*n);
    p = reshape(t,d,d,n);
    
    M = exp.(im*2*pi*p/k);
    return M
end

#
# Check all the posibles configurations 
# for a discrete  phases
#
function VariaConf(d,n,k)
    
    Mubval = 0;
    Mubval2 = 0;

    Ortval = 0
    Confval = 0;
    Confval2 = 0;

    for i in 0:k^(d^2*n)-1;
        M = GeneraOneMatConf(d,n,k,i)./sqrt(d);
        tv = Mubness(M);
        tb = IsBasis(M);
        if (tv < Mubval) & (tb<1e-13);
            Mubval = tv;
            Confval = i;
        end
        if (tv < Mubval) ;
            Mubval2 = tv;
            Confval2 = i;
            Ortval = tb;
        end
    end

    return [Mubval,Confval, Mubval2, Ortval, Confval2]
end


