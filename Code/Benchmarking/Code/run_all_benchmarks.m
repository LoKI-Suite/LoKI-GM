% ===============================================================
% run_all_benchmarks.m
% Script principal para executar todos os testes de benchmarking
% do LoKI-B: comparação malha variável vs uniforme e testes
% analíticos (Maxwelliana e Druyvesteyn)
% ===============================================================

clear; clc; close all;

% Save current directory and change to Code/ directory
% This ensures all relative paths work correctly for lokibcl
original_dir = pwd;
cd('..');  % Move to Code/ directory
addpath(pwd);  % Add Code/ to path for lokibcl function

fprintf('=================================================================\n');
fprintf('           LoKI-B COMPREHENSIVE BENCHMARKING SUITE             \n');
fprintf('=================================================================\n\n');
fprintf('Working directory: %s\n\n', pwd);

% ===============================================================
% CONFIGURAÇÃO DOS TESTES
% ===============================================================

% Testes a executar (pode desativar alguns para testes rápidos)
run_grid_comparison = true;      % Comparação malha variável vs uniforme
run_maxwellian_elastic = true;   % Teste Maxwelliano com colisões elásticas
run_maxwellian_ee = true;        % Teste Maxwelliano com colisões e-e  
run_maxwellian_const_v = true;   % Teste Maxwellian (nu=const)

% Parâmetros para teste de convergência (grid comparison)
N_values = [50, 100, 200, 400];  % Valores de cellNumber a testar
delta_u_values = [5e-4, 1e-3, 5e-3]; % Valores de firstEnergyStep a testar

% Criar diretório principal de resultados
output_base = 'Benchmarking/Output/comprehensive_benchmark';
if ~exist(output_base, 'dir')
    mkdir(output_base);
end

% Guardar configuração
config = struct();
config.date = datestr(now);
config.N_values = N_values;
config.delta_u_values = delta_u_values;
config.tests_run = struct('grid_comparison', run_grid_comparison, ...
                          'maxwellian_elastic', run_maxwellian_elastic, ...
                          'maxwellian_ee', run_maxwellian_ee, ...
                          'maxwellian_const_v', run_maxwellian_const_v);
save(fullfile(output_base, 'benchmark_config.mat'), 'config');

% ===============================================================
% TESTE 1: COMPARAÇÃO MALHA VARIÁVEL VS UNIFORME
% ===============================================================

if run_grid_comparison
    fprintf('\n=== TESTE 1: Comparação Malha Variável vs Uniforme ===\n\n');
    
    % 1A: Delta u fixo, variar N
    fprintf('1A: Delta u fixo (1e-3 eV), variando N...\n');
    for i = 1:length(N_values)
        N = N_values(i);
        
        % Malha variável
        fprintf('  Executando malha VARIÁVEL com N=%d...\n', N);
        try
            run_single_benchmark('benchmark_fixed_delta_u', 'variable', N, 1e-3, output_base);
            fprintf('  ✓ OK\n');
        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            fprintf('    Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
        
        % Malha uniforme
        fprintf('  Executando malha UNIFORME com N=%d...\n', N);
        try
            run_single_benchmark('benchmark_fixed_delta_u', 'uniform', N, 1e-3, output_base);
            fprintf('  ✓ OK\n');
        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            fprintf('    Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
    end
    
    % 1B: N fixo, variar delta u
    fprintf('\n1B: N fixo (200), variando delta u...\n');
    for i = 1:length(delta_u_values)
        delta_u = delta_u_values(i);
        
        % Malha variável
        fprintf('  Executando malha VARIÁVEL com delta_u=%.1e...\n', delta_u);
        try
            run_single_benchmark('benchmark_fixed_N', 'variable', 200, delta_u, output_base);
            fprintf('  ✓ OK\n');
        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            fprintf('    Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
        
        % Malha uniforme
        fprintf('  Executando malha UNIFORME com spacing=%.1e...\n', delta_u);
        try
            run_single_benchmark('benchmark_fixed_N', 'uniform', 200, delta_u, output_base);
            fprintf('  ✓ OK\n');
        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            fprintf('    Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
    end
end

% ===============================================================
% TESTE 2: MAXWELLIANA COM COLISÕES ELÁSTICAS (E/N = 0)
% ===============================================================

if run_maxwellian_elastic
    fprintf('\n=== TESTE 2: Maxwelliana com Colisões Elásticas ===\n\n');
    
    % Malha variável
    fprintf('  Executando malha VARIÁVEL...\n');
    try
        run_analytical_test('maxwellian_elastic', 'variable', output_base);
        fprintf(' OK\n');
    catch ME
        fprintf(' ERRO: %s\n', ME.message);
    end
    
    % Malha uniforme
    fprintf('  Executando malha UNIFORME...\n');
    try
        run_analytical_test('maxwellian_elastic', 'uniform', output_base);
        fprintf(' OK\n');
    catch ME
        fprintf(' ERRO: %s\n', ME.message);
    end
end

% ===============================================================
% TESTE 3: MAXWELLIANA COM COLISÕES E-E
% ===============================================================

if run_maxwellian_ee
    fprintf('\n=== TESTE 3: Maxwelliana com Colisões e-e ===\n\n');
    
    % Verificar se a cross section dummy existe
    fprintf('  Verificando cross section dummy...\n');
    if ~exist('Input/Dummy/H2_dummy_elastic.txt', 'file')
        fprintf('  AVISO: Cross section dummy não encontrada!\n');
        fprintf('  Por favor, execute manualmente antes de continuar:\n');
        fprintf('    >> generate_dummy_elastic\n\n');
        fprintf('  Pulando teste Maxwellian e-e...\n');
        run_maxwellian_ee = false;  % Pular este teste
    else
        fprintf('  ✓ Cross section dummy encontrada\n');
    end
end

if run_maxwellian_ee
    fprintf('  Executando malha VARIÁVEL...\n');
    try
        run_analytical_test('maxwellian_ee', 'variable', output_base);
        fprintf(' OK\n');
    catch ME
        fprintf(' ERRO: %s\n', ME.message);
    end
    
    fprintf('  Executando malha UNIFORME...\n');
    try
        run_analytical_test('maxwellian_ee', 'uniform', output_base);
        fprintf(' OK\n');
    catch ME
        fprintf(' ERRO: %s\n', ME.message);
    end
end

% ===============================================================
% TESTE 4: MAXWELLIAN (nu=const)
% ===============================================================

if run_maxwellian_const_v
    fprintf('\n=== TESTE 4: Maxwellian (nu=const) ===\n\n');
    
    % Verificar se a secção eficaz especial existe
    fprintf('  Verificando secção eficaz para nu=const...\n');
    if ~exist('Input/Maxwellian_const_v/constant_nu_elastic.txt', 'file')
        fprintf('  AVISO: Secção eficaz Maxwellian_const_v não encontrada!\n');
        fprintf('  Por favor, execute manualmente antes de continuar:\n');
        fprintf('    >> generate_maxwellian_const_v_cross_section\n\n');
        fprintf('  Pulando teste Maxwellian_const_v...\n');
    else
        fprintf('  ✓ Secção eficaz encontrada\n');
        
        % Malha variável
        fprintf('  Executando malha VARIÁVEL...\n');
        try
            run_analytical_test('maxwellian_const_v', 'variable', output_base);
            fprintf('  ✓ OK\n');
        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            fprintf('    Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
        
        % Malha uniforme
        fprintf('  Executando malha UNIFORME...\n');
        try
            run_analytical_test('maxwellian_const_v', 'uniform', output_base);
            fprintf('  ✓ OK\n');
        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            fprintf('    Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
    end
end

% ===============================================================
% FINALIZAÇÃO
% ===============================================================

fprintf('\n=================================================================\n');
fprintf('           BENCHMARKING COMPLETO                                \n');
fprintf('=================================================================\n');
fprintf('Resultados salvos em: %s\n', output_base);
fprintf('\nExecute analyze_all_benchmarks.m para analisar os resultados.\n\n');

% ===============================================================
% FUNÇÕES AUXILIARES
% ===============================================================

function run_single_benchmark(test_type, grid_type, N, delta_u, output_base)
    % Executa um teste individual de benchmarking
    
    % Determinar template de input
    if strcmp(test_type, 'benchmark_fixed_delta_u')
        if strcmp(grid_type, 'variable')
            template = 'Input/benchmark/benchmark_fixed_delta_u.in';
        else
            template = 'Input/benchmark/benchmark_uniform_fixed_delta_u.in';
        end
    elseif strcmp(test_type, 'benchmark_fixed_N')
        if strcmp(grid_type, 'variable')
            template = 'Input/benchmark/benchmark_fixed_N.in';
        else
            template = 'Input/benchmark/benchmark_uniform_fixed_N.in';
        end
    else
        error('Tipo de teste desconhecido: %s', test_type);
    end
    
    % Ler template
    fid = fopen(template, 'r');
    if fid == -1
        error('Não foi possível abrir template: %s', template);
    end
    content = fread(fid, '*char')';
    fclose(fid);
    
    % Modificar parâmetros
    content = regexprep(content, 'cellNumber: \d+', sprintf('cellNumber: %d', N));
    if strcmp(test_type, 'benchmark_fixed_N')
        if strcmp(grid_type, 'variable')
            content = regexprep(content, 'firstEnergyStep: [\d\.e-]+', sprintf('firstEnergyStep: %.1e', delta_u));
        else
            % Para uniforme, ajustar maxEnergy para manter resolução similar
            maxEnergy = delta_u * N;
            content = regexprep(content, 'maxEnergy: [\d\.]+', sprintf('maxEnergy: %.1f', maxEnergy));
        end
    end
    
    % Criar nome de arquivo temporário
    temp_file = sprintf('Input/benchmark/temp_%s_%s_N%d_du%.0e.in', test_type, grid_type, N, delta_u);
    
    % Modificar pasta de output no conteúdo
    output_folder = sprintf('%s_%s_N%d_du%.0e', test_type, grid_type, N, delta_u);
    content = regexprep(content, 'folder: \w+', sprintf('folder: %s/%s', output_base, output_folder));
    
    % Escrever arquivo temporário
    fid = fopen(temp_file, 'w');
    fprintf(fid, '%s', content);
    fclose(fid);
    
    % Executar simulação
    [~, filename, ext] = fileparts(temp_file);
    fprintf('    [DEBUG] Executando: %s%s\n', filename, ext);
    fprintf('    [DEBUG] Output folder: %s\n', output_folder);
    % lokibcl expects path relative to Input/ folder from Code/ directory
    lokibcl(fullfile('benchmark', [filename ext]));
    
    % Remover arquivo temporário
    delete(temp_file);
end

function run_analytical_test(test_name, grid_type, output_base)
    % Executa um teste analítico (Maxwelliana ou Maxwellian_const_v)
    
    % Determinar arquivo de input
    input_file = sprintf('Input/benchmark/benchmark_%s_%s.in', test_name, grid_type);
    
    if ~exist(input_file, 'file')
        error('Arquivo de input não encontrado: %s', input_file);
    end
    
    % Ler arquivo
    fid = fopen(input_file, 'r');
    content = fread(fid, '*char')';
    fclose(fid);
    
    % Modificar pasta de output
    output_folder = sprintf('benchmark_%s_%s', test_name, grid_type);
    content = regexprep(content, 'folder: benchmark_\w+', sprintf('folder: %s/%s', output_base, output_folder));
    
    % Criar arquivo temporário
    temp_file = sprintf('Input/benchmark/temp_%s_%s.in', test_name, grid_type);
    fid = fopen(temp_file, 'w');
    fprintf(fid, '%s', content);
    fclose(fid);
    
    % Executar simulação
    [~, filename, ext] = fileparts(temp_file);
    fprintf('    [DEBUG] Input original: %s\n', input_file);
    fprintf('    [DEBUG] Executando: %s%s\n', filename, ext);
    fprintf('    [DEBUG] Output folder: %s\n', output_folder);
    % lokibcl expects path relative to Input/ folder from Code/ directory
    lokibcl(fullfile('benchmark', [filename ext]));
    
    % Remover arquivo temporário
    delete(temp_file);
end

