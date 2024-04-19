
include("functions.jl")

#
# The idea is check by sequantial* searching modify  the phaes
# by example each entry has c*exp(i 2*pi/k), what happen when k
# is 2,3,4, etc. 
# starting by  d=2 for d basis 

function main()
    
    d = 2;
    n = d;
    kmax = 15;
    M = zeros(Complex,d,d,n);

    #M[:,:,1] = [1 1 ; 1  -1]./sqrt(2);
    #M[:,:,2] = [1 1;im  -im]./sqrt(2);    
    #val = Mubness(M);
    result = zeros(kmax-1,6);
    io = open("d$(d)-n$(n)-kmax$(kmax).dat","w");
    @printf(io,"# k   Mval1 conf1  Mval2 Orto2 conf2\n")

    for k in 2:kmax;
        result[k-1,1] = k;
        result[k-1,2:6] =  VariaConf(d,n,k);
        @printf(io,"%d  %f  %d  %f  %f  %d \n",k,result[k-1,2],result[k-1,3],result[k-1,4],result[k-1,5],result[k-1,2])
        @show(k)
    end
    Mubval = -n*(n-1)*sqrt(d^2-d);
    @show(Mubval)
    
    close(io)
    #return result;

end


t=main()
