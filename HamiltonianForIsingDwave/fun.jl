using LinearAlgebra

#
# Evalua ortogonalidad en vectores
#
#
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
#
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




