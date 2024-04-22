file = "Horth2n0k4.dat";
h = load(file);

Nqb = size(h,1);
H = zeros(Nqb+1,Nqb+1).*nan;
H(1:Nqb,1:Nqb) = h;

[X,Y] = meshgrid(1:Nqb+1,Nqb+1:-1:1);

pcolor(X,Y,H)



