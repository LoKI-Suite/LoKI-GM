% Script para gerar secção eficaz para teste Maxwellian (nu = const)
% Para obter nu = const, precisamos de sigma(v) ∝ 1/v
% nu = N * sigma(v) * v = N * sigma_0 * (v_0/v) * v = N * sigma_0 * v_0 = const
% Esta configuração resulta numa distribuição Maxwelliana, não Druyvesteyn

clear; clc;

fprintf('=== Generating Maxwellian (nu=const) Cross Section ===\n');

% Constantes
eV_to_J = 1.60218e-19;  % conversão eV para Joules
me = 9.10938e-31;       % massa do eletrão (kg)

% Parâmetros da secção eficaz
% Escolher valores REALISTAS para garantir nu = const
% Usar valores similares ao H2 real (~1e-19 m²)
v_0 = 1e6;              % velocidade de referência (m/s) - REDUZIDO
sigma_0 = 1e-19;        % secção eficaz de referência (m^2) - REDUZIDO para ser mais realista

% Criar array de energia
energy_eV = logspace(-3, 2, 200);  % de 1 meV a 10 eV
energy_J = energy_eV * eV_to_J;

% Calcular velocidade correspondente: E = (1/2) * m * v^2
velocity = sqrt(2 * energy_J / me);

% Secção eficaz: sigma(v) = sigma_0 * (v_0 / v)
% Mas na verdade, na base de dados é sigma(E), então:
% sigma(E) = sigma_0 * (v_0 / v(E)) = sigma_0 * v_0 * sqrt(m / (2*E))
cross_section = sigma_0 * v_0 * sqrt(me ./ (2 * energy_J));

% Converter para unidades típicas do LXCat (m^2)
% Já está em m^2

% Criar directório se não existir
if ~exist('../Input/Maxwellian_const_v', 'dir')
    mkdir('../Input/Maxwellian_const_v');
end

% Escrever ficheiro no formato LXCat
filename = '../Input/Maxwellian_const_v/constant_nu_elastic.txt';
fid = fopen(filename, 'w');

% Cabeçalho no estilo LXCat (seguindo o formato do H2_LXCat.txt)
fprintf(fid, 'LXCat, www.lxcat.net\n');
fprintf(fid, 'Generated on %s. All rights reserved.\n\n', datestr(now, 'dd mmm yyyy'));
fprintf(fid, 'RECOMMENDED REFERENCE FORMAT\n');
fprintf(fid, '- Maxwellian (nu=const) test database, custom cross section for constant collision frequency.\n\n');
fprintf(fid, 'CROSS SECTION DATA FORMAT\n');
fprintf(fid, 'Custom elastic cross section designed to ensure constant collision frequency nu = N*sigma*v = const.\n');
fprintf(fid, 'This is achieved by making sigma(v) proportional to 1/v.\n');
fprintf(fid, 'Note: This configuration results in a Maxwellian distribution, not a Druyvesteyn distribution.\n\n');
fprintf(fid, '********************************************************** H2 **********************************************************\n\n');
fprintf(fid, 'ELASTIC\n');
fprintf(fid, 'H2\n');
fprintf(fid, ' 2.743480e-4\n');  % m/M ratio (electron mass / H2 mass)
fprintf(fid, 'SPECIES: e / H2\n');
fprintf(fid, 'PROCESS: E + H2 -> E + H2, Elastic (Maxwellian nu=const test)\n');
fprintf(fid, 'PARAM.:  m/M = 0.000274348, sigma(v) proportional to 1/v\n');
fprintf(fid, 'COMMENT: [E + H2(X) -> E + H2(X), Elastic] Custom cross section for Maxwellian (nu=const) test.\n');
fprintf(fid, 'COMMENT: Designed to ensure constant collision frequency nu = N * sigma * v = const.\n');
fprintf(fid, 'COMMENT: Results in a Maxwellian distribution, not Druyvesteyn.\n');
fprintf(fid, 'UPDATED: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'COLUMNS: Energy (eV) | Cross section (m2)\n');
fprintf(fid, '-----------------------------\n');

% Escrever dados (formato exato do LXCat)
for i = 1:length(energy_eV)
    fprintf(fid, ' %.6e\t%.6e\n', energy_eV(i), cross_section(i));
end
fprintf(fid, '-----------------------------\n');

fclose(fid);

fprintf('Cross section file created: %s\n', filename);

% Plotar para verificação
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
loglog(energy_eV, cross_section * 1e20, 'b-', 'LineWidth', 2);
xlabel('Energy (eV)');
ylabel('Cross Section (10^{-20} m^2)');
title('Cross Section: \sigma(E)');
grid on;

subplot(1,3,2);
loglog(velocity, cross_section * 1e20, 'r-', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('Cross Section (10^{-20} m^2)');
title('Cross Section: \sigma(v)');
grid on;

subplot(1,3,3);
% Verificar que nu = N * sigma * v = constante
N = 133.32 / 1.380649e-23 / 300;  % densidade do gás a 133.32 Pa e 300 K (m^-3)
nu = N * cross_section .* velocity;
semilogx(energy_eV, nu, 'g-', 'LineWidth', 2);
xlabel('Energy (eV)');
ylabel('Collision Frequency \nu (s^{-1})');
title('Collision Frequency: \nu = N \sigma v');
grid on;

saveas(gcf, '../Input/Maxwellian_const_v/cross_section_verification.png');
fprintf('Verification plot saved: ../Input/Maxwellian_const_v/cross_section_verification.png\n');

% Verificar que nu é aproximadamente constante
fprintf('\nVerification:\n');
fprintf('  Min collision frequency: %.3e s^-1\n', min(nu));
fprintf('  Max collision frequency: %.3e s^-1\n', max(nu));
fprintf('  Mean collision frequency: %.3e s^-1\n', mean(nu));
fprintf('  Relative variation: %.2f%%\n', (max(nu)-min(nu))/mean(nu)*100);

fprintf('\n=== Cross section generation completed ===\n');

