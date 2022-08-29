L = 4;
gamma = 3.6;
d = 25;
Dv = 0.012;
Gamma = gamma*L^2/Dv;
a = 3.75;
figure(1)
for b = 0.001:0.01:0.2
    u0 = (b+1)/(a+b+1);
    fu = -a - 2*a*u0-b;
    gv = 2*a*u0-1-b;
    m1 = fu + d*gv;
    m2 = 4*d*b*(a+b+1);
    k2minus = Gamma/(2*d)*(m1-sqrt(m1^2-m2));
    k2plus = Gamma/(2*d)*(m1+sqrt(m1^2-m2));
    n2minus = k2minus/ pi^2;
    n2plus = k2plus/ pi^2;
    nminus = sqrt(n2minus);
    nplus = sqrt(n2plus);
    hold on
    scatter(b,nplus,20,'filled','b')
    scatter(b,nminus,20,'filled','k')
end
ylabel('n')
xlabel('b')
legend('n_{max}','n_{min}')