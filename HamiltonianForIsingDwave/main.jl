
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
    d = 2; n = 2; k = 4;
    Hu =  oneUpforVec(d,n,k,0,2);
    Ho  =  HOrth(d,n,k,-1,0);  
    Hi  = Hising(d,n,k,2);
    if n>1 
        Hm = HMubs(d,n,k,0,1)
    end
    
    cont1 = 0;
    cont2 = 0;
    for ik in 1:10;
    st ,E = findMS(Ho+Hu+Hm,41,1000,3);st'
    M = ConverStToVec(st,d,n,k);

    if  M != [];
        isO = CheckOrthM(M,d,n,k);
        isM = CheckMubsM(M,d,n,k);
        cont1 = cont1 +1;
   # @show(isO,isM)
    if (isO+isM) ==2;
        cont2 = cont2 +1;
    end
    end
    end
  # H = [-1  1  1 1 ;
  #      0  -1 1  1;
  #   11   0  0  -1 1;
  #      0  0  0  -1];
  # E = Energia(H);
  return cont1,cont2;

end

function Nminst()
    d = 2 ; n = 1; k = 2;
    Ho = HOrth(d,n,k,-1,0);
    Hi = Hising(d,n,k,1);
    Hu = oneUpforVec(d,n,k,0,2);
    
    B = collect(2:0.001:6);
    of = open("HooHu2-B2-3-d"*string(2)*"n"*string(n)*"k"*string(k)*".dat","w");
   
    for b in B;
        Em,t = findM(Ho+Hu,b);
        nt = size(t,1);
        mt = mean(sum.(VecConf.(t.-1,n*d*d^k)));
        @printf(of,"%f    %f   %f   %f\n",b,nt,Em,mt);
        println(b," " ,nt," ",Em," " ,mt)
    end
   
    close(of);
end 





function saveH()
    d = 2 ; n = 1; k =2 ;
    Ho = HOrth(d,n,k,-1,0);
    Hi = Hising(d,n,k,1);
    Hu = oneUpforVec(d,n,k,0,2);
    if n>1 
        Hm = HMubs(d,n,k,0,1);
    else
        Hm = zeros(n*k^d,n*k^d);
    end

    # Horht= Ortogonalidad(d,k,-1,0);
    # Hmub = Mubsness(d,k,0,1);
    #open("Hmat/Hoorht"*string(d)*"n"*string(n-1)*"k"*string(k)*".dat","w") do io
    #    writedlm(io,Horht)
    #end
    # open("Hmubs"*string(d)*"n"*string(n-1)*"k"*string(k)*".dat","w") do io
    #     writedlm(io,Hmub)
    # end  
 
   # corregir DhoHmHu esta mal lo sobre escribir por error  con DHooHmHu
    open("Hmat/DHooHmHu"*string(d)*"n"*string(n)*"k"*string(k)*".dat","w") do io
        writedlm(io,Ho+Hu+Hm)
    end

end


function Dsave();
    d = 2 ; n = 1; k = 2 ;

    Ho = HOrth(d,n,k,0,1);
    Hi = Hising(d,n,k,1);
    Hu = oneUpforVec(d,n,k,0,1);
    if n>1 
        Hm = HMubs(d,n,k,0,1);
    end


    open("DHoHu"*string(d)*"n"*string(n)*"k"*string(k)*".dat","w") do io
        writedlm(io,Ho+Hu)
    end
end
