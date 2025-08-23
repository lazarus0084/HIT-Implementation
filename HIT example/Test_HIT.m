clc;
clear;
close all;
alva.N = 300; 
alva.n = 30; 
alva.analysis = 'Full';
alva.bond = 'Bonded';
%-----------------------------Load-----------------------------------------------------------------------------------------


alva.q = [0.56 
          0.56] ;

alva.a = [106.6 
          106.6] ;

alva.Xl = [-155 0
           155 0]; 
%-----------------------------Original System-----------------------------------------------------------------------------------------
alva.E =  [3000 2000 150 65 45];E =  [3000 2000 150 65 45];
E_mif = [3000 2000 1.2*150 65 45];
alva.nu = [0.35 0.35 0.35 0.35 0.45];nu = [0.35 0.35 0.35 0.35 0.45];

layers = [40 95 250 200];

alva.zi = cumsum(layers);
all_sections = [30 50 250 200
                30 60 250 200
                40 90 250 200
                40 95 250 200
                40 105 250 200
                40 115 250 200];

%-----------------------------Homogenization------------------------------------------------------------------------------------------

alva.rfLayer = 4; nn = alva.rfLayer; %reinforced layer

alva.D_fifty = 6; %mean particle size of reinforced layer

alva = HIT(alva); % E, nu and Thickness of homogenised layer

%------------------------------HIT ANALYSIS-----------------------------

d =  alva.T_hit; %inter_layer_thickness 
T = layers(alva.rfLayer); %Base Layer Thickness

resolution = 3; % in percentage
step = (T)*resolution/100;
y = d/2:step: (T - d/2); % placement depths
% y = T*[0.33 0.5 0.67];
 
for j = 1: size(all_sections, 1) 
for i = 1: size(y,2)
layers = all_sections(j,:); 
upper_layer = y(i) - d/2 ;

lower_layer = T - y(i) - d/2;

new_layers = [layers(1:nn-1) upper_layer d lower_layer layers(nn+1:end)];     
alva.zi = cumsum(new_layers);
alva.E = [E(1:nn) alva.E_hit E((nn):end)]; % Layer Young's moduli [MPa]
alva.nu = [nu(1:nn) alva.nu_hit nu((nn):end)]; % Layer Poisson's ratio [-]
alva.Xd = [0 0 (alva.zi(2) + 0.01);
           155 0 (alva.zi(2) + 0.01);
           0 0 (alva.zi(end) - 0.01);
           155 0 (alva.zi(end) - 0.01)];

alva = init_LET(alva);
rutting_strains = max(abs([alva.epsz(3) alva.epsz(4)]));

crack_strains = max(abs([alva.epsx(1) alva.epsx(2)]));
epsilon_t = crack_strains ; % Assume a reasonable strain value (in microstrain)
M_rm = 3000;  % Assume modulus ratio

C = 0.5;

% Compute reps
Nr_rf(i,j) = (1.4 * 10^-8 * (1/rutting_strains)^4.5337)/10^6;
Nf_rf(i,j) = (1.6064 * C * 1e-4 * (1 / epsilon_t)^3.89 * (1 / M_rm)^0.854)/10^6;

%---------------------MIF section--------------------------
layers = all_sections(j,:); 
alva.zi = cumsum(layers);
alva.E = E_mif; % Layer Young's moduli [MPa]
alva.nu = nu; % Layer Poisson's ratio [-]
alva.Xd = [0 0 (alva.zi(2) + 0.01);
           155 0 (alva.zi(2) + 0.01);
           0 0 (alva.zi(end) - 0.01);
           155 0 (alva.zi(end) - 0.01)];

alva = init_LET(alva);
rutting_strains = max(abs([alva.epsz(3) alva.epsz(4)]));

crack_strains = max(abs([alva.epsx(1) alva.epsx(2)]));
epsilon_t = crack_strains ; % Assume a reasonable strain value (in microstrain)
M_rm = 3000;  % Assume modulus ratio

C = 0.5;

% Compute reps
Nr_MIF(i,j) = (1.4 * 10^-8 * (1/rutting_strains)^4.5337)/10^6;
Nf_MIF(i,j) = (1.6064 * C * 1e-4 * (1 / epsilon_t)^3.89 * (1 / M_rm)^0.854)/10^6;
%---------------------unreinforced section--------------------------
layers = all_sections(j,:); 
alva.zi = cumsum(layers);
alva.E = E; % Layer Young's moduli [MPa]
alva.nu = nu; % Layer Poisson's ratio [-]
alva.Xd = [0 0 (alva.zi(2) + 0.01);
           155 0 (alva.zi(2) + 0.01);
           0 0 (alva.zi(end) - 0.01);
           155 0 (alva.zi(end) - 0.01)];

alva = init_LET(alva);
rutting_strains = max(abs([alva.epsz(3) alva.epsz(4)]));

crack_strains = max(abs([alva.epsx(1) alva.epsx(2)]));
epsilon_t = crack_strains ; % Assume a reasonable strain value (in microstrain)
M_rm = 3000;  % Assume modulus ratio

C = 0.5;

% Compute reps
Nr_unrf(i,j) = (1.4 * 10^-8 * (1/rutting_strains)^4.5337)/10^6;
Nf_unrf(i,j) = (1.6064 * C * 1e-4 * (1 / epsilon_t)^3.89 * (1 / M_rm)^0.854)/10^6;
end
end

x_axis = 100*y'./T;
output_folder = 'C:\Users\arman\Downloads\results';
% ----------------------- Plot Fatigue ---------------------------
for section = 1:6
    fig = figure;
    plot(x_axis, Nf_rf(:,section), 'b-', 'LineWidth', 2); hold on;
    plot(x_axis, Nf_unrf(:,section), 'r--', 'LineWidth', 2);
    plot(x_axis, Nf_MIF(:,section), 'g-.', 'LineWidth', 2);
    xlabel('Normalized Depth (% of base thickness)');
    ylabel('Fatigue Repetitions to Failure (Millions)');
    title(sprintf('Section %d: Fatigue - Reinforced vs Unreinforced vs MIF', section));
    legend('Reinforced', 'Unreinforced', 'MIF', 'Location', 'best');
    grid on;

    % Save fatigue figure
    filename = fullfile(output_folder, sprintf('Section_%d_Fatigue.png', section));
    saveas(fig, filename);
end

% ----------------------- Plot Rutting ---------------------------

for section = 1:6
    fig = figure;
    plot(x_axis, Nr_rf(:,section), 'b-', 'LineWidth', 2); hold on;
    plot(x_axis, Nr_unrf(:,section), 'r--', 'LineWidth', 2);
    plot(x_axis, Nr_MIF(:,section), 'g-.', 'LineWidth', 2);
    xlabel('Normalized Depth (% of base thickness)');
    ylabel('Rutting Repetitions to Failure (Millions)');
    title(sprintf('Section %d: Rutting - Reinforced vs Unreinforced vs MIF', section));
    legend('Reinforced', 'Unreinforced', 'MIF', 'Location', 'best');
    grid on;

    % Save rutting figure
    filename = fullfile(output_folder, sprintf('Section_%d_Rutting.png', section));
    saveas(fig, filename);

end
