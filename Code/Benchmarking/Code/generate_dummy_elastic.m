% Script para gerar secção eficaz "dummy" extremamente pequena
% Para teste Maxwellian e-e (queremos APENAS colisões e-e, sem colisões com gás)

clear; clc;

fprintf('=== Generating Dummy Elastic Cross Section ===\n');

% Criar directório se não existir
if ~exist('../Input/Dummy', 'dir')
    mkdir('../Input/Dummy');
end

% Criar array de energia
energy_eV = logspace(-3, 1, 200);  % de 1 meV a 10 eV

% Secção eficaz EXTREMAMENTE pequena (praticamente zero, mas não zero)
% Isto permite que o LoKI reconheça o H2, mas as colisões são desprezáveis
cross_section = ones(size(energy_eV)) * 1e-30;  % 10^-30 m² (10^10 vezes menor que H2 real)

% Escrever ficheiro no formato LXCat
filename = '../Input/Dummy/H2_dummy_elastic.txt';
fid = fopen(filename, 'w');

% Cabeçalho no estilo LXCat (seguindo o formato do H2_LXCat.txt)
fprintf(fid, 'LXCat, www.lxcat.net\n');
fprintf(fid, 'Generated on %s. All rights reserved.\n\n', datestr(now, 'dd mmm yyyy'));
fprintf(fid, 'RECOMMENDED REFERENCE FORMAT\n');
fprintf(fid, '- Dummy database for e-e collision test.\n\n');
fprintf(fid, 'CROSS SECTION DATA FORMAT\n');
fprintf(fid, 'Dummy elastic cross section (effectively zero) to allow H2 recognition without gas collisions.\n');
fprintf(fid, 'This allows testing pure e-e collision effects.\n\n');
fprintf(fid, '********************************************************** H2 **********************************************************\n\n');
fprintf(fid, 'ELASTIC\n');
fprintf(fid, 'H2\n');
fprintf(fid, ' 2.743480e-4\n');  % m/M ratio (electron mass / H2 mass)
fprintf(fid, 'SPECIES: e / H2\n');
fprintf(fid, 'PROCESS: E + H2 -> E + H2, Elastic (dummy - effectively zero)\n');
fprintf(fid, 'PARAM.:  m/M = 0.000274348, dummy cross section\n');
fprintf(fid, 'COMMENT: [E + H2(X) -> E + H2(X), Elastic] Dummy cross section for e-e collision test.\n');
fprintf(fid, 'COMMENT: Cross section is effectively zero (1e-30 m²) to eliminate gas collisions.\n');
fprintf(fid, 'UPDATED: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'COLUMNS: Energy (eV) | Cross section (m2)\n');
fprintf(fid, '-----------------------------\n');

% Escrever dados (formato exato do LXCat)
for i = 1:length(energy_eV)
    fprintf(fid, ' %.6e\t%.6e\n', energy_eV(i), cross_section(i));
end
fprintf(fid, '-----------------------------\n');

fclose(fid);

fprintf('Dummy cross section file created: %s\n', filename);
fprintf('  Cross section value: %.3e m² (effectively zero)\n', cross_section(1));
fprintf('  Ratio to H2 real (~1e-19): %.3e (10^%.0f smaller)\n', ...
    cross_section(1)/1e-19, log10(cross_section(1)/1e-19));

fprintf('\n=== Dummy cross section generation completed ===\n');

