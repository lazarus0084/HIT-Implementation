clc;
clear;
close all;
alva.N = 300;
alva.n = 30;
alva.analysis = 'Full';
alva.bond = 'Bonded';


d = 60;
M_rm = 3000;  % Assume modulus ratio
C = 0.5;
%-----------------------------Load-----------------------------------------------------------------------------------------
alva.q = [0.56
        0.56];
alva.a = [106.6
        106.6];
alva.Xl = [-155 0
         155 0];

correction_factor = 1.59;

% --- Material properties ---
alva.E =  [3000 2500 75 65 60];
alva.nu = [0.35 0.35 0.35 0.35 0.45];

% All layer sets
orginal_all_sections = [30 50 250 180
                        30 60 260 200
                        40 90 300 200
                        40 95 300 250
                        40 105 350 250];
aaa = cumsum(orginal_all_sections(1,:));
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
section = 5; 
base =  3;
reduction = 79;
d = 60; nn=3;
placement_ratio = 0.3;
all_sections =          [30 50 250 180
                        30 60 260 200
                        40 90 300 200
                        40 95 300 250
                        40 105 350 250];
all_sections(section, base) = all_sections(section, base) - reduction;

n_sections = size(orginal_all_sections, 1);

Nr_og = zeros(n_sections, 1); % Store rutting for each section
Nf_og = zeros(n_sections, 1); % Store rutting for each section
    for i = 1:size(orginal_all_sections, 1)
       layer = orginal_all_sections(i,:);  % Fixed variable name (hybrid_all_sections1)
        alva.zi = cumsum(layer);  % Cumulative sum for layer thickness
        
        % Define geometry points
        alva.Xd = [0 0 (alva.zi(2) + 0.01);
                   155 0 (alva.zi(2) + 0.01);
                   0 0 (alva.zi(end) - 0.01);
                   155 0 (alva.zi(end) - 0.01)];
        
        % Initialize LET (make sure init_LET is defined correctly)
        alva = init_LET(alva);  % Ensure the init_LET function is available
        
        % Calculate maximum vertical strain
        rutting_strains = max(abs([alva.epsz(3), alva.epsz(4)]));
        
        % Calculate rutting life (million repetitions)
       Nr_og(i) = (1.4e-8 * (1 / rutting_strains)^4.5337) / 1e6; % original
    epsilon_t = max(abs([alva.epsx(1) alva.epsx(2)]));

M_rm = 3000;  % Assume modulus ratio
C = 0.5;
% Compute reps
Nf_og(i) = (1.6064 * C * 1e-4 * (1 / epsilon_t)^3.89 * (1 / M_rm)^0.854)/10^6;
    end


% Optionally, display rutting value2
disp('Rutting Reps:');
disp(Nr_og);
disp('Fatigue Reps:');
disp(Nf_og);


for j = 1: size(all_sections, 1)
T = all_sections(j,nn); %Base Layer Thickness
y = placement_ratio*T;
  
layers = all_sections(j,:);
upper_layer = y - d/2 ;
lower_layer = T - y - d/2;
new_layers(j,:) = [layers(1:nn-1) upper_layer d lower_layer layers(nn+1:end)];  
end



   
    layer = new_layers(section,:);
    
    alva.zi = cumsum(layer);
    alva.E = [3000	2500 75	152.03	75	65	60];
    alva.nu = [0.35	0.35 0.35 0.35 0.20 0.35 0.45];
    % Define geometry points
    alva.Xd = [0 0 (alva.zi(2) + 0.01);
               155 0 (alva.zi(2) + 0.01);
               0 0 (alva.zi(end) - 0.01);
               155 0 (alva.zi(end) - 0.01)];
    
    % Initialize LET (your function)
    alva = init_LET(alva);
    
    % Calculate maximum vertical strain
    rutting_strains = max(abs([alva.epsz(3), alva.epsz(4)]));
    
    % Calculate rutting life (million repetitions)
   rutting_values = ((1.4e-8 * (1 / rutting_strains)^4.5337) / 1e6)*correction_factor % reduced layer
    
    epsilon_t = max(abs([alva.epsx(1) alva.epsx(2)]));

  %tbr = rutting_values/(Nr_og(section)*correction_factor);
 



% Optionally, display rutting value2
bcr = reduction*100/orginal_all_sections(section,3);
disp('BCR with correction factor (%):');
disp(bcr);
%disp('TBR with correction factor (%):');

%disp(tbr);
