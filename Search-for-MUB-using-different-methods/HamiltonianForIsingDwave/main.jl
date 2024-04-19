
include("fun.jl")

#
# Queremos contruir una busqueda discreta de mubs
# donde cada vector formado solo por fases conocidas
# formado por un conjunto de qbits 
#
# Queremos que remos intentar con fases difertentes por 
# ejemplo 4 fases por entrada haciendo 16 confinfiguraciones
# por vector, y necesitaremos 16 qbits, donde un qbit en el 
# estado up representa una de las 16 configurciones 
#


function main()
    d = 2; n = 1; k = 2;
    Hu =  oneUpforVec(d,n,k,0,1);
    Ho  =  HOrth(d,n,k,0,1);
    Hi  = Hising(d,n,k,2);
    if n>1 
        Hm = HMubs(d,n,k,0,1)
    end
    
  # H = [-1  1  1 1 ;
  #      0  -1 1  1;
  #      0  0  -1 1;
  #      0  0  0  -1];
  # E = Energia(H);
  

end

function Nminst()
    d = 2 ; n = 2; k = 2;
    Ho = HOrth(d,n,k,0,1);
    Hi = Hising(d,n,k,2);
    Hu = oneUpforVec(d,n,k,0,1);
    
    B = collect(-7:0.01:7);
    of = open("HoHu-B-d"*string(2)*"n"*string(n)*"k"*string(k)*".dat","w");
   
    for b in B;
        Em,t = findM(Ho+Hu,b);
        nt = size(t,1);
        mt = mean(sum.(VecConf.(t.-1,n*d*d^k)));
        @printf(of,"%f    %f   %f   %f\n",b,nt,Em,mt);
        println(b,"" ,nt," ",Em," " ,mt)
    end
   
    close(of);
end 





function saveH()
    d = 2 ; n = 2; k = 4;
    Ho = HOrth(d,n,k,0,1);
    Hi = Hising(d,n,k,1);
    Hu = oneUpforVec(d,n,k,0,1);
    if n>1 
        Hm = HMubs(d,n,k,0,1)
    end

    # Horht= Ortogonalidad(d,k,0,1);
    # Hmub = Mubsness(d,k,0,1);
    #  open("Horht"*string(d)*"n"*string(n-1)*"k"*string(k)*".dat","w") do io
    #     writedlm(io,Horht)
    # end
    # open("Hmubs"*string(d)*"n"*string(n-1)*"k"*string(k)*".dat","w") do io
    #     writedlm(io,Hmub)
    # end  
 
   
    open("HoHm"*string(d)*"n"*string(n)*"k"*string(k)*".dat","w") do io
        writedlm(io,Ho+Hm)
    end

end
