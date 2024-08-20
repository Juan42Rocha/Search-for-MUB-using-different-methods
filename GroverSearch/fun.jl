using LinearAlgebra
using Printf
using Statistics
using Random
using DelimitedFiles



#  en logisim  and es solo eun espacio 
#   a & b ==  a b
#  
#  y or es con una suma  
#  a or b == a + b
#
# la negacion es  este simbolo
# las tablas en logisim es de izquierda a derecha 
# aqui estoy hacienod las cosas de derecha a izquierda 
# en logisim  2 en binario es  10 
# aqui es     2 es             01
function ConverSt2logisimCir(sts);
 num = size(sts,1);
 var = "hgfedcba";
 t= digits.(sts,base=2,pad=8);

 for k in 1:num;
     for i in 1:8;
       if i > 1;
          print("*")
       end
       if t[k][i]==0
          #print("~")
          print("n")
       end
       print(var[i:i]*" ")
     end
     print("+ ")
 end


end


# tabla pra mubsness
function CcircuitMubsness(d,k)
    Out =  [];

    for nv1 in 0:k^d-1;
        v1 = VecConFases(d,k,nv1);
        for nv2 in 0:k^d-1;
            v2 = VecConFases(d,k,nv2);
            o12 = IsMubs(v1,v2);
            if o12 == 1;      # si v1 y v2 no son ortogonales no importa v3 y v4
                 val = FushNum2(d,k,nv1,nv2);
                 push!(Out,val);
             end           
            
        end
    end
    return Out
end

# tabla para ortogonallidad
function  CcircuitOrth(n,k);
    Out =  [];

    for nv1 in 0:k^d-1;
        v1 = VecConFases(d,k,nv1);
        for nv2 in 0:k^d-1;
            v2 = VecConFases(d,k,nv2);
            o12 = IsOrthogonality(v1,v2);
            if o12 == 1;      # si v1 y v2 no son ortogonales no importa v3 y v4
                 val = FushNum2(d,k,nv1,nv2);
                 push!(Out,val);
             end           
            
        end
    end
    return Out
end


# tabla total final
function COracle(d,n,k);
    OrOnes =  [];

    for nv1 in 0:k^d-1;
        v1 = VecConFases(d,k,nv1);
        for nv2 in 0:k^d-1;
            v2 = VecConFases(d,k,nv2);
            o12 = IsOrthogonality(v1,v2);
            if o12 == 1;      # si v1 y v2 no son ortogonales no importa v3 y v4
            for nv3 in 0:k^d-1;
                v3 = VecConFases(d,k,nv3);
                m13 =  IsMubs(v1,v3);
                m23 =  IsMubs(v2,v3);
                if  (m13 & m23) == 1;   # si v3 no es mub con v1 y v2 no importa v4
                for nv4 in 0:k^d-1;
                    v4 = VecConFases(d,k,nv4);
                    o34 =  IsOrthogonality(v3,v4);
                    m14 =  IsMubs(v1,v4);
                    m24 =  IsMubs(v2,v4);   
                    if (o34 & m14 & m24 )==1 ;
                        val = FushNum(d,n,k,nv1,nv2,nv3,nv4);
                        push!(OrOnes,val);
                    end           
                end
                end
            end

            end 
        end
    end
    return OrOnes;
end

#
# fusiona dos vectos 
function  FushNum2(d,k,nv1,nv2);
    st1 = digits(nv1,base=2,pad=4)';
    st2 = digits(nv2,base=2,pad=4)';

    st = [ st2 st1];
    pw = 2 .^collect(0:7)';
    val = sum(st .* pw);
    @show(st,val)
    return val
end

# fusiona los 4 vectores
function  FushNum(d,n,k,nv1,nv2,nv3,nv4);
    st1 = digits(nv1,base=2,pad=4)';
    st2 = digits(nv2,base=2,pad=4)';
    st3 = digits(nv3,base=2,pad=4)';
    st4 = digits(nv4,base=2,pad=4)';
    st = [st4 st3 st2 st1];
    pw = 2 .^collect(0:16-1)';
    val = sum(st .* pw);
    return val
end

# revisa si es una mubb n=2 d=2 k=4
function check(val1);
    t = digits(val1,base=2,pad=16)';
    ts = reshape(t,(4,4))';
    pw = 2 .^collect(0:4-1);
    nst = ts*pw;

    v = VecConFases.(d,k,ts*pw);
    
    o12 = abs(v[1]'*v[2]);
    o34 = abs(v[3]'*v[4]);

    M13 = abs(v[1]'*v[3]);
    M14 = abs(v[1]'*v[4]);
    M23 = abs(v[2]'*v[3]);
    M24 = abs(v[2]'*v[4]);    

    @show(o12,o34,M13,M14,M23,M24)

end


#
# Crea el vector las configuracion de las fases
# solo para un vectore la matriz
# d es la dimension
# n el numero de matrices
# k el nuemro de fases
# conf el numero de configuracion
function VecConFases(d,k,conf);
    t = digits(conf,base=k,pad=d)
    return exp.(im*2*pi.*t/k)./sqrt(d)
end


#
# Crea el vector de eigen estados de sigma_z  para 
# Nqb qbits de todo el sistema
# k es el numero de estado
# Nqb el numero de qbits
function VecConf(k,Nqb);
    t = digits(k,base=2,pad=Nqb).*2 .-1;
    return t
end
#
# Evalua  Mubsness, el producto entre dos 
# vectores cumple la relga de mubs
#  v1, v2 son vectores de norma 1
function IsMubs(v1,v2);
    d = size(v1,1);
    flag = 0;
    t = 1/sqrt(d) -abs(v1'*v2);
    if abs(t)< 1e-14;
	flag = 1
    else
	flag = 0;
    end
    return flag;
end
#
# Evalua ortogonalidad en vector
#
#  v1 ,v2 son vectores de norma 1
function  IsOrthogonality(v1,v2);
    t = abs(v1'*v2);
    flag = 0;
    if abs(t)< 1e-12;
	flag = 1;
    else
	flag = 0;
    end
    return flag;
end




