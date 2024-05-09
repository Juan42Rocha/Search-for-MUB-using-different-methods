using LinearAlgebra
using Printf
using Statistics
using Random
using DelimitedFiles

#
# Funcion para entronar los configuraciones de minima energia
#
#function findM(d,n,k,a,b)
function findM(H,Jsi );
 # H1 = [-1  0.5 0    0    0 ;
 #      0  -1   0.5  0    0 ;
 #      0   0  -1    0.5  0 ;
 #      0   0   0    -1   0.5;
 #       0.5 0   0     0  -1];
 # H2 = [0  0   0    0    ;
 #      0   0  0.5  0.5   ;
 #      0   0   0    0     ;
 #      0   0   0    0   ]; 
 #    H = oneUpforVec(d,n,k,a,b);
  # H = HOrth(d,n,k,a,b);
   
    E = Energias(H,Jsi);
    #@show H
    Em = minimum(E);
    minst = findall(isone,abs.(E.-Em).<1e-12);
    return Em,minst;

end
function findMS(H,Jsi,Mts,temp);
    t = 0;
    Nqb, Nqb = size(H);
    E =[] #zeros(Nqb*Mts+1,1);
   # t = rand(collect(0:Nqb^2),1);
   # @show(t)
    st0 = rand([-1 1],Nqb); #VecConf.(t, Nqb )[1];
    st = rand([-1 1],Nqb);#VecConf(88, Nqb );
   # @show(st)
    Eo = Energiast(st0,H,Jsi);
    E  = push!(E,Eo);
    for k in 2:Mts*Nqb+1;
       st = change1qb(st0);
       En = Energiast(st,H,Jsi);
       if  En < Eo 
           st0 = st;
           Eo = En;
        E  = push!(E,Eo);
        #@show(k/(Nqb*Mts),Eo)
       elseif abs(En-Eo)< 1e-12;
           if rand([0 1],1)==0
               st0 = st;
               Eo = En;
           end  
       else
          if rand()> 1/(1+exp(-(En-Eo)/temp))
           st0 = st;
           Eo = En;
          end
       end
       push!(E, Eo);
    end
   
    return st0,E;
end
function change1qb(st);
    t = copy(st);
    iqb = rand(collect(1:size(st,1)),1);
    t[iqb] = -t[iqb];
    return t;
end

function findM2(H,Jsi,d,n,k)
    Nqb,Nqb = size(H);
    Emin = 0;
    E  = zeros(k^(d^2*n));
    E2  = zeros(k^(d^2*n));
    Ismin  = ones(Bool,k^(d^2*n),1);

    st1= zeros(k^(d^2*n),Nqb) ;
    st2= zeros(k^(d^2*n),Nqb) ;
    for i in 0:(k^(d^2*n)-1) 
       nst = digits(i,base=k*d*n,pad=d*n);
       st = buildst(nst,Nqb,n,d,k);
       E[i+1] = Energiast(st,H,Jsi);
       st1[i+1,:] = st';
       @show(st')
    end
         
    for iq in 1:Nqb
        stv = st2;
        stv[:,iq] = -1 .*stv[:,iq];
        for i in 0:(k^(d^2*n)-1) 
          E2[i+1] = Energiast(stv[i+1,:],H,Jsi);
        end
        Ismin =  (E .< E2) .& Ismin;
       @show(E .< E2)
    end

    return E,stv,Ismin
end
function buildst(nst,Nqb,n,d,k)
    st = ones(Nqb,1)*.-1;
    nst = nst.+1;
    aux = collect(0:d*n-1) .*k^d;
    t = nst + aux;
    st[t] .= -1 .*st[t];
    return st;
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
    # comentar o descomentar para  todos los vectores
     for i in 1:subd:Nqb
         H[i:subd+i-1,i:subd+i-1] = Hsub
     end
    return H;
end


#
# Evalua ortogonalidad en vectores
#
#  v1 ,v2 son vectores de norma 1
function  IsOrthogonality(v1,v2)
    t = abs(v1'*v2);
    flag = 0;
    if abs(t)< 1e-12;
	flag = 1;
    else
	flag = 0;
    end
    return flag;
end
#
# Genera el las intereaccioens para  tener orgogonalidad 
#
function Ortogonalidad(d,k,a,b)
    kd = k^d;  # numero de qubits para representar un vector
    Hsub = zeros(kd,kd);
    for i in 0:kd-1;
	for j in i:kd-1;
	#    Hsubt = ones(kd,kd).*-1;
	    v1 = VecConFases(d,k,i);  # genra par de vectoes 
	    v2 = VecConFases(d,k,j);
	    t = IsOrthogonality(v1,v2);
            if t==1   # ahora se se rellena si cumple la propidad
		# Hsubt[i+1,:] = -1 .*Hsubt[i+1,:];
		# Hsubt[:,j+1] = -1 .*Hsubt[:,j+1];
		# Hsubt[Hsubt.==1] .= b;
		# Hsubt[Hsubt.==-1] .= a;
		# Hsub += Hsubt;
                Hsub[i+1,j+1] = a;
                Hsub[j+1,i+1] = a;
		# si no lo comple se puede  rellenar con el inverso
	    else 
                Hsub[i+1,j+1] = b;
                Hsub[j+1,i+1] = b;
            end
	end
    end
    return Hsub

end
function HOrth(d,n,k,a,b)
    Hs = Ortogonalidad(d,k,a,b);
    t = triu(ones(d,d),1);
    t2 = I(n);
    H = kron(t,Hs);
    H = kron(t2,H);
    return H
end

#
# Evalua  Mubsness, el producto entre dos 
# vectores cumple la relga de mubs
#  v1, v2 son vectores de norma 1
function IsMubs(v1,v2)
    #d,t = size(v1)
    flag = 0;
    t = 1/sqrt(d) -abs(v1'*v2);
    if abs(t)< 1e-14;
	flag = 1;
    else
	flag = 0;
    end
    return flag;

end
function Mubsness(d,k,a,b)
    kd = k^d;  # numero de qubits para representar un vector
    Hsub = zeros(kd,kd);
    for i in 0:kd-1;
	for j in i:kd-1;
	#    Hsubt = ones(kd,kd).*-1;
	    v1 = VecConFases(d,k,i);  # genra par de vectoes 
	    v2 = VecConFases(d,k,j);
	    t = IsMubs(v1,v2);
            if t==1   # ahora se se rellena si cumple la propidad
		# Hsubt[i+1,:] = -1 .*Hsubt[i+1,:];
		# Hsubt[:,j+1] = -1 .*Hsubt[:,j+1];
		# Hsubt[Hsubt.==1] .= b;
		# Hsubt[Hsubt.==-1] .= a;
		# Hsub += Hsubt;
                Hsub[i+1,j+1] = a;
                Hsub[j+1,i+1] = a;
		# si no lo comple se puede  rellenar con el inverso
	    else 
                Hsub[i+1,j+1] = b;
                Hsub[j+1,i+1] = b;
            end
	end
    end
    return Hsub
end
function HMubs(d,n,k,a,b)
    Hs = Mubsness(d,k,a,b);
    t1 = ones(d,d);
    t2 = triu(ones(n,n),1);
    H = kron(t1,Hs);
    H = kron(t2,H);
    return H
end


#
# Crea el vector de eigen estados de sigma_z  para 
# Nqb qbits de todo el sistema
# k es el numero de estado
# Nqb el numero de qbits
function VecConf(k,Nqb)
    t = digits(k,base=2,pad=Nqb).*2 .-1;
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

function Ia(n)
    II = zeros(n,n);
    for i in  1:n
        II[i,n+1-i] = 1;
    end
    return II
end


#
# Busqueda de las configuraciones de mínima energía
# H es para   Js' H Js  = sum_{i,j}   H_{i,j}*s_i *s_j
# jsi es para   sum(Js)*jsi = sum_j  jsi*Js_j
function Energias(H,Jsi)
    Nqb, Nqb = size(H);
    E = zeros(2^Nqb);
    for i = 0:2^Nqb-1
        st = VecConf(i,Nqb);
        E[i+1] = Energiast(st,H,Jsi);
    end
    return E
end
function Energiast(st,H,Jsi);
    return sum(st'*H*st) + sum(st)*Jsi;
end

#
# Crea interacciones para  ising
#
function Hising(d,n,k,b);
    H = zeros(d^k,d^k);  
    for i in 1:d^k-1
        H[i,i+1] = b;
    end
    H[d^k,1] = b;
    H = kron(I(d),H);
    H = kron(I(n),H);
    return H;
end


#
# Crea interacciones para Ortogonalidad old
#
function Ortogonalidadold(d,k,a,b)
    kd = k^d;  # numero de qubits para representar un vector
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
                 Hsub[i+1,j+1] = a;
		# si no lo comple se puede  rellenar con el inverso
	    #else 
            #    Hsub[i+1,j+1] = 1;
            end
	end
    end
    return Hsub

end



function CheckOneSol(st,d,n,k);
    M = ConverStToVec(st,d,n,k);
    
    for ij in 1:n;
        

    end   

end
function ConverStToVec(st,d,n,k);
    v = reshape(st,(k^d,n*d))
    v = (v.+1)./2;
    aux = sum(sum(v,dims=1) .==1);
    M = [];    

    if aux ==n*d;
       PosOnes = findall(isone,v);
       for ij in 1:n;
           m = zeros(Complex,d,d);
           for ik in 1:d;
               idxv = PosOnes[ik][1];
#@show(idxv,d,n,k)
               m[:,ik] = VecConFases(d,k,idxv-1);
           end
           push!(M,m);
       end
    else
      a = 0;
    end

    return  M
end
function CheckOrthM(M,d,n,k);
    flag = 1;
    for ixn in n;
        M1 = M[ixn];
        t = sum(abs.(M1'*M1) - I(d));
        if abs(t)>1e-14;
            flag = 0*flag;
        end 
    end
    return flag;
end
function CheckMubsM(M,d,n,k);
   flag = 1; 

    for inx in 1:n-1;
        M1 = M[inx];
        for jnx in inx+1:n;
            M2 = M[jnx];
            t = sum( abs.(M1'*M2) - ones(d,d)./sqrt(d) );
            if abs(t)> 1e-14;
                flag = flag*0;
            end 
        end 
    end
    return flag;
end
