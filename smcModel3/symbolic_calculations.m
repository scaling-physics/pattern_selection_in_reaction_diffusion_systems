% Symbolic calculations
clear all
syms fu fv gu gv a b Gamma u0 v0 d

u0 = (b+1)/(a+b+1);
v0 = a/(a+b+1);
c = u0+v0;

fu = -a-2*a*u0-b;
fv = 1-2*a*u0;
gu = a+2*a*u0;
gv = 2*a*u0-1-b;

km2 = (Gamma/(d-1))*((d+1)*(-fv*gu/d)^(0.5) - fu + gv );
chr = latex(simplify(km2))