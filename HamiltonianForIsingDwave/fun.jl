using LinearAlgebra

#
# Evalua ortogonalidad en vectores
#
#  v1 ,v2 son vectores de norma 1
function  IsOrthogonality(v1,v2)
    t = abs(v1'*v2);
    flag = 0;
    if abs(t)< 1e-14;
	flag = 1;
    else
	flag = 0;
    end
    return flag;
end

#
# Evalua  Mubsness, el producto entre dos 
# vectores cumple la relga de mubs
#  v1, v2 son vectores de norma 1
function IsMubs(v1,v2)
    d,t = size(v1)
    flag = 0;
    t = 1/sqrt(d) -abs(v1'*v2);
    if abs(t)< 1e-14;
	flag = 1;
    else
	flag = 0;
    end
    return flag;

end

#
# Funcion para entronar los configuraciones de minima energia
#
function findM(d,n,k,a,b)
   # H = HOrth(d,n,k,a,b) 
    H = oneUpforVec(d,n,k,a,b);
    E = Energias(H);
    #@show H
    Em = minimum(E);
    minst = findall(isone,E.==Em);
    return minst

end


#
# Dado un numero de fases k busca construir las interacciones
# que castigan o premian el tener solo un estado up por todos
# los qubits que representan un vector
#  d es la dimension 
#  n el numero de matrices
#  k el numero de fases 
#  
function oneUpforVec(d,n,k,a,b)
    subd = k^d;
    Nqb = subd*n*d;
    H = zeros(Nqb,Nqb);
    Hsub = zeros(subd,subd);

    for i in 1:subd-1;
	Hsub[i,i+1:subd] .= b;
	Hsub[i,i] = a;
    end
    Hsub[subd,subd] = a;

    for i in 1:subd:Nqb
	H[i:subd+i-1,i:subd+i-1] = Hsub
    end
    return H
end

#
# Genera el las intereaccioens para  tener orgogonalidad 
#
function Ortogonalidad(d,k,a,b)
    kd = k^d;
    Hsub = zeros(kd,kd);
    for i in 0:kd-1;
	for j in i:kd-1;
	    Hsubt = ones(kd,kd).*-1;
	    v1 = VecConFases(d,k,i);  # genra par de vectoes 
	    v2 = VecConFases(d,k,j);
	    t = IsOrthogonality(v1,v2);
            if t==1   # ahora se se rellena si cumple la propidad
		Hsubt[i+1,:] = -1 .*Hsubt[i+1,:];
		Hsubt[:,j+1] = -1 .*Hsubt[:,j+1];
		Hsubt[Hsubt.==1] .= b;
		Hsubt[Hsubt.==-1] .= a;
		Hsub += Hsubt;
		# si no lo comple se puede  rellenar con el inverso
	    end
	end
    end
    return Hsub

end

function HOrth(d,n,k,a,b)
    Hs = Ortogonalidad(d,k,a,b);
    t = [0 1 ; 1 0];
    t2 = I(n);
    H = kron(t,Hs);
    H = kron(t2,H);
    return H

end

#
# Busqueda de las configuraciones de mínima energía
#
function Energias(H)
    Nqb, Nqb = size(H);
    E = zeros(2^Nqb)
    for i = 0:2^Nqb-1
	st = VecConf(i,Nqb);
	E[i+1] = st'*H*st;
    end
    return E
end

#
# Crea el vector de eigen estados de sigma_z  para 
# Nqb qbits de todo el sistema
# k es el numero de estado
# Nqb el numero de qbits
function VecConf(k,Nqb)
    t = digits(k,base=2,pad=Nqb)*2 .-1;
    return t
end

#
# Crea el vector las configuracion de las fases
# solo para un vectore la matriz
# d es la dimension
# n el numero de matrices
# k el nuemro de fases
# conf el numero de configuracion
function VecConFases(d,k,conf)
    t = digits(conf,base=k,pad=d)
    return exp.(im*2*pi.*t/k)./sqrt(d)
end

