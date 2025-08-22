function alva = HIT(alva)

%-----------------------------homogenization------------------------------------------------------------------------------------------
% Define parameters %mm
D50 = alva.D_fifty;
alpha = 10;
nu1 = 0.3; 
E1 = 1300; 
hom_interlayer = alpha * D50;
n = alva.rfLayer;
E2 = alva.E(n);
nu2 = alva.nu(n);

RibThickness = 2.5; 
h1 = RibThickness; %mm
SoilCompositeThickness = alpha * D50; 
h2 =SoilCompositeThickness/2;%mm


[Ez, nuzx] = homogenize(E1, E2, nu1, nu2, h1, h2);


h1 = SoilCompositeThickness/2 + 2.5;
h2 = SoilCompositeThickness/2;
nu1 = nuzx; nu2 = 0.35;
E1 = Ez; E2 = 319;
[alva.E_hit, alva.nu_hit] = homogenize(E1, E2, nu1, nu2, h1, h2);
alva.T_hit = alpha*D50;


end