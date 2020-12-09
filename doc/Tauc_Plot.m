%% Steven E. Bopp Materials Science & Engineering Nov 26, 2020

%% Tauc Plot for ZrN

ZrN_R = table2array(ZrN2R); ZrN_T = table2array(ZrN2T); ZrN_A = table2array(ZrN2A);

%Alpha_h_nu=ZrN_A(:,2).*(1240./((2*3.14159).*ZrN_A(:,1)));
Alpha_h_nu=ZrN_A(:,2).*(1240./ZrN_A(:,1)); % Absorption Coef. * 1240 ev*nm / wavelength 
h_nu=1240./ZrN_A(:,1); % Calculate energy from wavelength 1240 ev*nm / wavelength
r1 = Alpha_h_nu.^(1/(1/2)); % Tauc plot (αhν)^1/(1/2) Direct Allowed
r2 = Alpha_h_nu.^(1/(3/2)); % Tauc plot (αhν)^1/(3/2) Direct Forbidden
r3    = Alpha_h_nu.^(1/2);   % Tauc plot (αhν)^1/(2)   Indirect Allowed
r4    = Alpha_h_nu.^(1/3);   % Tauc plot (αhν)^1/(3)   Indirect Forbidden

x1 = linspace(0,3); % Make a vector for rootfinding as the x-value

r1_cut = r1(17:55); h_nu_cut_for_r1 = h_nu(17:55); % Choose linear bounds y=0.8471*x-1.7452, root at x=2.06021
p1=polyfit(h_nu_cut_for_r1,r1_cut,1); y1 = polyval(p1,x1); % First-order fit

r2_cut = r2(22:70); h_nu_cut_for_r2 = h_nu(22:70); % Choose linear bounds y=0.5072*x-0.5527, root at x=1.08971
p2=polyfit(h_nu_cut_for_r2,r2_cut,1); y2 = polyval(p2,x1); % First-order fit

r3_cut = r3(30:80); h_nu_cut_for_r3 = h_nu(30:80); % Choose linear bounds y=0.4138*x-0.2564, root at x=0.619623
p3=polyfit(h_nu_cut_for_r3,r3_cut,1); y3 = polyval(p3,x1); % First-order fit

r4_cut = r4(32:60); h_nu_cut_for_r4 = h_nu(32:60); % Choose linear bounds y=0.3201*x+0.0458, root at x=-0.14308
p4=polyfit(h_nu_cut_for_r4,r4_cut,1); y4 = polyval(p4,x1); % First-order fit

plot(h_nu,r1,h_nu,r2,h_nu,r3,h_nu,r4,x1,y1,'k--',x1,y2,'k--',x1,y3,'k--',x1,y4,'k--')
xlabel('hν'); ylabel('(αhν)^1/r (Arb. Units)'); title('Tauc Plot for Hf:ZrxNy'); 
legend('r=1/2 Direct Allowed        2.06 eV','r=3/2 Direct Forbidden     1.09 eV','r=2    Indirect Allowed      0.62 eV','r=3    Indirect Forbidden -0.14 eV','Location','northwest')
legend('boxoff')
xlim([0 3.25]); ylim([0,1]);