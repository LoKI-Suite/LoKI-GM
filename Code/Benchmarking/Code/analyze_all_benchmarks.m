% ===============================================================
% analyze_all_benchmarks.m
% Script para analisar todos os resultados de benchmarking
% Compara com soluções analíticas e gera relatórios detalhados
% ===============================================================

clear; clc; close all;

fprintf('=================================================================\n');
fprintf('          ANÁLISE DE RESULTADOS DE BENCHMARKING                \n');
fprintf('=================================================================\n\n');

% Diretório base dos resultados
% IMPORTANTE: Os resultados estão em Output/Output/comprehensive_benchmark/ 
% (há um nesting duplo devido ao regexprep no run_all_benchmarks.m)
output_base = 'Output/Output/comprehensive_benchmark';
config_base = 'Output/comprehensive_benchmark';

if ~exist(output_base, 'dir')
    fprintf('Aviso: Pasta %s não encontrada.\n', output_base);
    fprintf('Procurando em Output/benchmark_* ...\n');
    output_base = 'Output';  % Fallback para resultados diretos
end

% Carregar configuração
config_file = fullfile(config_base, 'benchmark_config.mat');
if exist(config_file, 'file')
    load(config_file, 'config');
    fprintf('Configuração carregada de: %s\n', config.date);
else
    fprintf('Aviso: Arquivo de configuração não encontrado. Usando padrões.\n');
    config.N_values = [50, 100, 200, 400];
    config.delta_u_values = [5e-4, 1e-3, 5e-3];
end

% Criar diretório para figuras
figures_dir = fullfile(config_base, 'figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% ===============================================================
% ANÁLISE 1: COMPARAÇÃO MALHA VARIÁVEL VS UNIFORME
% ===============================================================

fprintf('\n=== ANÁLISE 1: Comparação Malha Variável vs Uniforme ===\n\n');

results_grid = analyze_grid_comparison(output_base, config);

if ~isempty(results_grid)
    % Plotar resultados de convergência
    plot_grid_convergence(results_grid, figures_dir);
end

% ===============================================================
% ANÁLISE 2: TESTES ANALÍTICOS
% ===============================================================

fprintf('\n=== ANÁLISE 2: Comparação com Soluções Analíticas ===\n\n');

% 2A: Maxwelliana com colisões elásticas
fprintf('2A: Maxwelliana com colisões elásticas...\n');
results_max_elastic = analyze_maxwellian_elastic(output_base);
if ~isempty(results_max_elastic)
    plot_maxwellian_comparison(results_max_elastic, figures_dir, 'elastic');
end

% 2B: Maxwelliana com colisões e-e
fprintf('2B: Maxwelliana com colisões e-e...\n');
results_max_ee = analyze_maxwellian_ee(output_base);
if ~isempty(results_max_ee)
    plot_maxwellian_comparison(results_max_ee, figures_dir, 'ee');
end

% 2C: Maxwellian_const_v
fprintf('2C: Maxwellian_const_v...\n');
results_maxwellian_const_v = analyze_maxwellian_const_v(output_base);
if ~isempty(results_maxwellian_const_v)
    plot_maxwellian_const_v_comparison(results_maxwellian_const_v, figures_dir);
end

% ===============================================================
% RELATÓRIO FINAL
% ===============================================================

fprintf('\n=== Gerando Relatório Final ===\n');
generate_report(config_base, results_grid, results_max_elastic, results_max_ee, results_druy);

fprintf('\n=================================================================\n');
fprintf('          ANÁLISE COMPLETA                                      \n');
fprintf('=================================================================\n');
fprintf('Figuras salvas em: %s\n', figures_dir);
fprintf('Relatório salvo em: %s\n\n', fullfile(config_base, 'benchmark_report.txt'));

% ===============================================================
% FUNÇÕES DE ANÁLISE
% ===============================================================

function results = analyze_grid_comparison(output_base, config)
    % Analisa resultados de comparação de malhas
    
    results = struct();
    results.fixed_delta_u = struct('variable', [], 'uniform', []);
    results.fixed_N = struct('variable', [], 'uniform', []);
    
    % Análise: delta u fixo, N variável
    for i = 1:length(config.N_values)
        N = config.N_values(i);
        
        % Malha variável
        folder_var = sprintf('%s/benchmark_fixed_delta_u_variable_N%d_du1e-03', output_base, N);
        data_var = read_swarm_parameters(folder_var);
        if ~isempty(data_var)
            results.fixed_delta_u.variable = [results.fixed_delta_u.variable; struct('N', N, 'data', data_var)];
        end
        
        % Malha uniforme
        folder_uni = sprintf('%s/benchmark_fixed_delta_u_uniform_N%d_du1e-03', output_base, N);
        data_uni = read_swarm_parameters(folder_uni);
        if ~isempty(data_uni)
            results.fixed_delta_u.uniform = [results.fixed_delta_u.uniform; struct('N', N, 'data', data_uni)];
        end
    end
    
    % Análise: N fixo, delta u variável
    for i = 1:length(config.delta_u_values)
        delta_u = config.delta_u_values(i);
        
        % Malha variável
        folder_var = sprintf('%s/benchmark_fixed_N_variable_N200_du%.0e', output_base, delta_u);
        data_var = read_swarm_parameters(folder_var);
        if ~isempty(data_var)
            results.fixed_N.variable = [results.fixed_N.variable; struct('delta_u', delta_u, 'data', data_var)];
        end
        
        % Malha uniforme
        folder_uni = sprintf('%s/benchmark_fixed_N_uniform_N200_du%.0e', output_base, delta_u);
        data_uni = read_swarm_parameters(folder_uni);
        if ~isempty(data_uni)
            results.fixed_N.uniform = [results.fixed_N.uniform; struct('delta_u', delta_u, 'data', data_uni)];
        end
    end
    
    fprintf('  Resultados carregados: %d testes com delta_u fixo, %d com N fixo\n', ...
        length(results.fixed_delta_u.variable) + length(results.fixed_delta_u.uniform), ...
        length(results.fixed_N.variable) + length(results.fixed_N.uniform));
end

function results = analyze_maxwellian_elastic(output_base)
    % Analisa teste Maxwelliano com colisões elásticas
    
    results = struct();
    T_gas_K = 300;  % Temperatura do gás
    T_gas_eV = kelvin_to_eV(T_gas_K);
    
    % Ler resultados de malha variável
    folder_var = sprintf('%s/benchmark_maxwellian_elastic_variable', output_base);
    [eedf_var, energy_var] = read_eedf(folder_var);
    data_var = read_swarm_parameters(folder_var);
    
    % Ler resultados de malha uniforme
    folder_uni = sprintf('%s/benchmark_maxwellian_elastic_uniform', output_base);
    [eedf_uni, energy_uni] = read_eedf(folder_uni);
    data_uni = read_swarm_parameters(folder_uni);
    
    if isempty(eedf_var) && isempty(eedf_uni)
        fprintf('  Aviso: Nenhum resultado encontrado para Maxwelliana elástica\n');
        results = [];
        return;
    end
    
    % Calcular solução analítica
    energy_range = [0, max([max(energy_var), max(energy_uni)])];
    [f0_analytical, energy_analytical] = analytical_maxwellian(T_gas_eV, energy_range, 1000);
    
    % Calcular erros relativos
    if ~isempty(eedf_var)
        f0_var_interp = interp1(energy_analytical, f0_analytical, energy_var, 'linear');
        error_var = calculate_relative_error(eedf_var, f0_var_interp);
    else
        error_var = NaN;
    end
    
    if ~isempty(eedf_uni)
        f0_uni_interp = interp1(energy_analytical, f0_analytical, energy_uni, 'linear');
        error_uni = calculate_relative_error(eedf_uni, f0_uni_interp);
    else
        error_uni = NaN;
    end
    
    results.T_gas_eV = T_gas_eV;
    results.analytical.energy = energy_analytical;
    results.analytical.f0 = f0_analytical;
    results.variable.energy = energy_var;
    results.variable.f0 = eedf_var;
    results.variable.data = data_var;
    results.variable.error = error_var;
    results.uniform.energy = energy_uni;
    results.uniform.f0 = eedf_uni;
    results.uniform.data = data_uni;
    results.uniform.error = error_uni;
    
    fprintf('  Erro relativo - Variável: %.2f%%, Uniforme: %.2f%%\n', error_var*100, error_uni*100);
end

function results = analyze_maxwellian_ee(output_base)
    % Analisa teste Maxwelliano com colisões e-e
    
    results = struct();
    E_N_Td = 10;  % Campo reduzido usado no teste
    
    % Ler resultados
    folder_var = sprintf('%s/benchmark_maxwellian_ee_variable', output_base);
    [eedf_var, energy_var] = read_eedf(folder_var);
    data_var = read_swarm_parameters(folder_var);
    
    folder_uni = sprintf('%s/benchmark_maxwellian_ee_uniform', output_base);
    [eedf_uni, energy_uni] = read_eedf(folder_uni);
    data_uni = read_swarm_parameters(folder_uni);
    
    if isempty(eedf_var) && isempty(eedf_uni)
        fprintf('  Aviso: Nenhum resultado encontrado para Maxwelliana e-e\n');
        results = [];
        return;
    end
    
    % Estimar Te dos resultados (energia média)
    % Para Maxwellian: <E> = (3/2) * Te → Te = (2/3) * <E>
    if ~isempty(data_uni) && isfield(data_uni, 'mean_energy')
        T_e_eV = data_uni.mean_energy * (2/3);
    elseif ~isempty(data_var) && isfield(data_var, 'mean_energy')
        T_e_eV = data_var.mean_energy * (2/3);
    
    else
        T_e_eV = 1.0;  % Valor padrão
        fprintf('  Aviso: Usando Te = 1 eV como padrão\n');
    end
    
    % Calcular solução analítica
    energy_range = [0, max([max(energy_var), max(energy_uni)])];
    [f0_analytical, energy_analytical] = analytical_maxwellian(T_e_eV, energy_range, 1000);
    
    


    % Calcular erros
    if ~isempty(eedf_var)
        f0_var_interp = interp1(energy_analytical, f0_analytical, energy_var, 'linear');
        error_var = calculate_relative_error(eedf_var, f0_var_interp);
    else
        error_var = NaN;
    end
    
    if ~isempty(eedf_uni)
        f0_uni_interp = interp1(energy_analytical, f0_analytical, energy_uni, 'linear');
        error_uni = calculate_relative_error(eedf_uni, f0_uni_interp);
    else
        error_uni = NaN;
    end
    
    results.T_e_eV = T_e_eV;
    results.E_N_Td = E_N_Td;
    results.analytical.energy = energy_analytical;
    results.analytical.f0 = f0_analytical;
    results.variable.energy = energy_var;
    results.variable.f0 = eedf_var;
    results.variable.data = data_var;
    results.variable.error = error_var;
    results.uniform.energy = energy_uni;
    results.uniform.f0 = eedf_uni;
    results.uniform.data = data_uni;
    results.uniform.error = error_uni;
    
    fprintf('  Te estimado: %.2f eV\n', T_e_eV);
    fprintf('  Erro relativo - Variável: %.2f%%, Uniforme: %.2f%%\n', error_var*100, error_uni*100);
end

function results = analyze_maxwellian_const_v(output_base)
    % Analisa teste Maxwellian_const_v (nu=const)
    
    results = struct();
    E_N_Td = 100;  % Campo reduzido usado no teste
    
    % Ler resultados
    folder_var = sprintf('%s/benchmark_maxwellian_const_v_variable', output_base);
    [eedf_var, energy_var] = read_eedf(folder_var);
    data_var = read_swarm_parameters(folder_var);
    
    folder_uni = sprintf('%s/benchmark_maxwellian_const_v_uniform', output_base);
    [eedf_uni, energy_uni] = read_eedf(folder_uni);
    data_uni = read_swarm_parameters(folder_uni);
    
    if isempty(eedf_var) && isempty(eedf_uni)
        fprintf('  Aviso: Nenhum resultado encontrado para Maxwellian_const_v\n');
        results = [];
        return;
    end
    
    % % Estimar Te dos resultados
    % if ~isempty(data_var) && isfield(data_var, 'mean_energy')
    %     T_eff_eV = data_var.mean_energy * (2/3);
    % elseif ~isempty(data_uni) && isfield(data_uni, 'mean_energy')
    %     T_eff_eV = data_uni.mean_energy * (2/3);
    % else
    %     T_eff_eV = 2.0;  % Valor padrão
    %     fprintf('  Aviso: Usando T_eff = 2 eV como padrão\n');
    % end
    
    % Calcular solução analítica
    energy_range = [0, max([max(energy_var), max(energy_uni)])];
    [f0_analytical, energy_analytical] = analytical_maxwellian_const_v(energy_range, 2000);
    
    % Calcular erros
    if ~isempty(eedf_var)
        f0_var_interp = interp1(energy_analytical, f0_analytical, energy_var, 'linear');
        error_var = calculate_relative_error(eedf_var, f0_var_interp);
    else
        error_var = NaN;
    end
    
    if ~isempty(eedf_uni)
        f0_uni_interp = interp1(energy_analytical, f0_analytical, energy_uni, 'linear');
        error_uni = calculate_relative_error(eedf_uni, f0_uni_interp);
    else
        error_uni = NaN;
    end
    
    %results.T_eff_eV = T_eff_eV;
    results.E_N_Td = E_N_Td;
    results.analytical.energy = energy_analytical;
    results.analytical.f0 = f0_analytical;
    results.variable.energy = energy_var;
    results.variable.f0 = eedf_var;
    results.variable.data = data_var;
    results.variable.error = error_var;
    results.uniform.energy = energy_uni;
    results.uniform.f0 = eedf_uni;
    results.uniform.data = data_uni;
    results.uniform.error = error_uni;
    
    %fprintf('  T_eff estimado: %.2f eV\n', T_eff_eV);
    fprintf('  Erro relativo - Variável: %.2f%%, Uniforme: %.2f%%\n', error_var*100, error_uni*100);
end

% ===============================================================
% FUNÇÕES AUXILIARES
% ===============================================================

function data = read_swarm_parameters(folder)
    % Lê parâmetros de swarm do arquivo lookUpTableSwarm.txt
    
    swarm_file = fullfile(folder, 'lookUpTableSwarm.txt');
    if ~exist(swarm_file, 'file')
        data = [];
        return;
    end
    
    try
        tbl = readtable(swarm_file, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
        data = struct();
        data.mean_energy = tbl{1, 9};      % MeanE(eV)
        data.drift_velocity = tbl{1, 4};   % DriftVelocity(ms^-1)
        data.diffusion_coeff = tbl{1, 2};  % RedDiff((ms)^-1)
        data.mobility = tbl{1, 3};         % Mobility((ms)^-1Td^-1)
    catch
        data = [];
    end
end

function [eedf, energy] = read_eedf(folder)
    % Lê EEDF do arquivo eedf.txt
    
    eedf_file = fullfile(folder, 'eedf.txt');
    if ~exist(eedf_file, 'file')
        eedf = [];
        energy = [];
        return;
    end
    
    try
        tbl = readtable(eedf_file, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
        energy = tbl{:, 1};  % Energy(eV)
        eedf = tbl{:, 2};    % EEDF(eV^-(3/2))
    catch
        eedf = [];
        energy = [];
    end
end

function rel_error = calculate_relative_error(f_numerical, f_analytical)
    % Calcula erro relativo médio entre soluções numérica e analítica
    
    % % Remover zeros e pontos com valor muito baixo
    % idx = (f_analytical > max(f_analytical)*1e-6) & (f_numerical > 0);
    % 
    % if sum(idx) < 10
    %     rel_error = NaN;
    %     return;
    % end
    
    % Erro relativo ponto a ponto
    pointwise_error = abs(f_numerical - f_analytical) ./ f_analytical;
    
    % Erro médio (com peso proporcional ao valor)
    %weights = f_analytical / sum(f_analytical);
    tamanho = size(f_analytical,1);
    rel_error = sum(pointwise_error) / tamanho;
end

function T_eV = kelvin_to_eV(T_K)
    % Converte temperatura de Kelvin para eV
    k_B = 1.380649e-23;  % Constante de Boltzmann (J/K)
    e = 1.602176634e-19; % Carga elementar (C)
    T_eV = k_B * T_K / e;
end

% ===============================================================
% FUNÇÕES DE PLOTAGEM
% ===============================================================

function plot_grid_convergence(results, figures_dir)
    % Plota resultados de convergência das malhas
    
    figure('Position', [100, 100, 1400, 900]);
    
    % Extrair dados - delta u fixo
    if ~isempty(results.fixed_delta_u.variable)
        N_var = [results.fixed_delta_u.variable.N];
        mean_E_var = arrayfun(@(x) x.data.mean_energy, results.fixed_delta_u.variable);
        drift_v_var = arrayfun(@(x) x.data.drift_velocity, results.fixed_delta_u.variable);
    else
        N_var = []; mean_E_var = []; drift_v_var = [];
    end
    
    if ~isempty(results.fixed_delta_u.uniform)
        N_uni = [results.fixed_delta_u.uniform.N];
        mean_E_uni = arrayfun(@(x) x.data.mean_energy, results.fixed_delta_u.uniform);
        drift_v_uni = arrayfun(@(x) x.data.drift_velocity, results.fixed_delta_u.uniform);
    else
        N_uni = []; mean_E_uni = []; drift_v_uni = [];
    end
    
    % Plot 1: Mean Energy vs N
    subplot(2,3,1);
    if ~isempty(N_var)
        semilogx(N_var, mean_E_var, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Variable');
        hold on;
    end
    if ~isempty(N_uni)
        semilogx(N_uni, mean_E_uni, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Uniform');
    end
    xlabel('N (cellNumber)');
    ylabel('Mean Energy (eV)');
    title('Convergence: Mean Energy');
    legend('Location', 'best');
    grid on;
    
    % Plot 2: Drift Velocity vs N
    subplot(2,3,2);
    if ~isempty(N_var)
        semilogx(N_var, drift_v_var, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Variable');
        hold on;
    end
    if ~isempty(N_uni)
        semilogx(N_uni, drift_v_uni, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Uniform');
    end
    xlabel('N (cellNumber)');
    ylabel('Drift Velocity (m/s)');
    title('Convergence: Drift Velocity');
    legend('Location', 'best');
    grid on;
    
    % Plot 3: Relative Error vs N
    subplot(2,3,3);
    if length(mean_E_var) > 1
        ref_E_var = mean_E_var(end);
        rel_err_var = abs(mean_E_var - ref_E_var) / ref_E_var * 100;
        semilogy(N_var, rel_err_var, 'go-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Variable');
        hold on;
    end
    if length(mean_E_uni) > 1
        ref_E_uni = mean_E_uni(end);
        rel_err_uni = abs(mean_E_uni - ref_E_uni) / ref_E_uni * 100;
        semilogy(N_uni, rel_err_uni, 'ms-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Uniform');
    end
    xlabel('N (cellNumber)');
    ylabel('Relative Error (%)');
    title('Relative Error vs N_{max}');
    legend('Location', 'best');
    grid on;
    
    % Extrair dados - N fixo
    if ~isempty(results.fixed_N.variable)
        du_var = [results.fixed_N.variable.delta_u];
        mean_E_var2 = arrayfun(@(x) x.data.mean_energy, results.fixed_N.variable);
        drift_v_var2 = arrayfun(@(x) x.data.drift_velocity, results.fixed_N.variable);
    else
        du_var = []; mean_E_var2 = []; drift_v_var2 = [];
    end
    
    if ~isempty(results.fixed_N.uniform)
        du_uni = [results.fixed_N.uniform.delta_u];
        mean_E_uni2 = arrayfun(@(x) x.data.mean_energy, results.fixed_N.uniform);
        drift_v_uni2 = arrayfun(@(x) x.data.drift_velocity, results.fixed_N.uniform);
    else
        du_uni = []; mean_E_uni2 = []; drift_v_uni2 = [];
    end
    
    % Plot 4: Mean Energy vs delta u
    subplot(2,3,4);
    if ~isempty(du_var)
        semilogx(du_var, mean_E_var2, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Variable');
        hold on;
    end
    if ~isempty(du_uni)
        semilogx(du_uni, mean_E_uni2, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Uniform');
    end
    xlabel('\Delta u (eV)');
    ylabel('Mean Energy (eV)');
    title('Convergence: Mean Energy vs \Delta u');
    legend('Location', 'best');
    grid on;
    
    % Plot 5: Drift Velocity vs delta u
    subplot(2,3,5);
    if ~isempty(du_var)
        semilogx(du_var, drift_v_var2, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Variable');
        hold on;
    end
    if ~isempty(du_uni)
        semilogx(du_uni, drift_v_uni2, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Uniform');
    end
    xlabel('\Delta u (eV)');
    ylabel('Drift Velocity (m/s)');
    title('Convergence: Drift Velocity vs \Delta u');
    legend('Location', 'best');
    grid on;
    
    % Salvar figura
    saveas(gcf, fullfile(figures_dir, 'grid_convergence_analysis.png'));
    savefig(gcf, fullfile(figures_dir, 'grid_convergence_analysis.fig'));
    fprintf('  Figura salva: grid_convergence_analysis.png\n');
end

function plot_maxwellian_comparison(results, figures_dir, test_type)
    % Plota comparação com Maxwelliana analítica
    
    figure('Position', [100, 100, 1400, 500]);
    
    % Plot 1: EEDF comparison
    subplot(1,3,1);
    
    % Debug info
    fprintf('    [DEBUG] Analytical EEDF range: %.3e - %.3e\n', min(results.analytical.f0), max(results.analytical.f0));
    if ~isempty(results.variable.f0)
        fprintf('    [DEBUG] Variable EEDF range: %.3e - %.3e\n', min(results.variable.f0), max(results.variable.f0));
    end
    if ~isempty(results.uniform.f0)
        fprintf('    [DEBUG] Uniform EEDF range: %.3e - %.3e\n', min(results.uniform.f0), max(results.uniform.f0));
    end
    
    semilogy(results.analytical.energy, results.analytical.f0, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Analytical');
    hold on;
    if ~isempty(results.variable.f0)
        semilogy(results.variable.energy, results.variable.f0, 'b--', 'LineWidth', 2, 'DisplayName', 'Variable');
    end
    if ~isempty(results.uniform.f0)
        semilogy(results.uniform.energy, results.uniform.f0, 'r:', 'LineWidth', 2, 'DisplayName', 'Uniform');
    end
    xlabel('Energy (eV)');
    ylabel('EEDF (eV^{-3/2})');
    if strcmp(test_type, 'elastic')
        title(sprintf('Maxwellian (T = %.3f eV) - Elastic Collisions', results.T_gas_eV));
    else
        title(sprintf('Maxwellian (T_e = %.2f eV) - e-e Collisions', results.T_e_eV));
    end
    legend('Location', 'best');
    grid on;
    
    % Plot 2: Relative difference
    subplot(1,3,2);
    if ~isempty(results.variable.f0)
        f0_var_interp = interp1(results.analytical.energy, results.analytical.f0, results.variable.energy, 'linear');
        rel_diff_var = (f0_var_interp - results.variable.f0) ./ f0_var_interp * 100;
        %idx_var = results.analytical.f0 > max(results.analytical.f0)*1e-6;
        plot(results.variable.energy, rel_diff_var, 'b-', 'LineWidth', 2, 'DisplayName', 'Variable');
        hold on;
    end
    if ~isempty(results.uniform.f0)
        f0_uni_interp = interp1(results.analytical.energy, results.analytical.f0, results.uniform.energy, 'linear');
        rel_diff_uni = (f0_uni_interp - results.uniform.f0) ./ f0_uni_interp * 100;
        %idx_uni = results.analytical.f0 > max(results.analytical.f0)*1e-6;
        plot(results.uniform.energy, rel_diff_uni, 'r--', 'LineWidth', 2, 'DisplayName', 'Uniform');
    end
    xlabel('Energy (eV)');
    ylabel('Relative Error (%)');
    title('Error vs Analytical Solution');
    legend('Location', 'best');
    grid on;
    yline(0, 'k:', 'LineWidth', 1);
    
    % Plot 3: Error summary
    subplot(1,3,3);
    errors = [];
    labels = {};
    if ~isempty(results.variable.error) && ~isnan(results.variable.error)
        errors = [errors, results.variable.error*100];
        labels = [labels, {'Variable'}];
    end
    if ~isempty(results.uniform.error) && ~isnan(results.uniform.error)
        errors = [errors, results.uniform.error*100];
        labels = [labels, {'Uniform'}];
    end
    if ~isempty(errors)
        bar(errors);
        set(gca, 'XTickLabel', labels);
        ylabel('Mean Relative Error (%)');
        title('Error Summary');
        grid on;
    end
    
    % Salvar figura
    filename = sprintf('maxwellian_%s_comparison.png', test_type);
    saveas(gcf, fullfile(figures_dir, filename));
    savefig(gcf, fullfile(figures_dir, strrep(filename, '.png', '.fig')));
    fprintf('  Figura salva: %s\n', filename);
end

function plot_maxwellian_const_v_comparison(results, figures_dir)
    % Plota comparação com Maxwellian_const_v analítica
    
    figure('Position', [100, 100, 1400, 500]);
    
    % Plot 1: EEDF comparison
    subplot(1,3,1);
    
    % Debug info
    fprintf('    [DEBUG] Analytical EEDF range: %.3e - %.3e\n', min(results.analytical.f0), max(results.analytical.f0));
    if ~isempty(results.variable.f0)
        fprintf('    [DEBUG] Variable EEDF range: %.3e - %.3e (variation: %.2f%%)\n', ...
            min(results.variable.f0), max(results.variable.f0), ...
            (max(results.variable.f0) - min(results.variable.f0))/mean(results.variable.f0)*100);
    end
    if ~isempty(results.uniform.f0)
        fprintf('    [DEBUG] Uniform EEDF range: %.3e - %.3e (variation: %.2f%%)\n', ...
            min(results.uniform.f0), max(results.uniform.f0), ...
            (max(results.uniform.f0) - min(results.uniform.f0))/mean(results.uniform.f0)*100);
    end
    
    semilogy(results.analytical.energy, results.analytical.f0, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Analytical');
    hold on;
    if ~isempty(results.variable.f0)
        semilogy(results.variable.energy, results.variable.f0, 'b--', 'LineWidth', 2, 'DisplayName', 'Variable');
    end
    if ~isempty(results.uniform.f0)
        semilogy(results.uniform.energy, results.uniform.f0, 'r:', 'LineWidth', 2, 'DisplayName', 'Uniform');
    end
    xlabel('Energy (eV)');
    ylabel('EEDF (eV^{-3/2})');
    %title(sprintf('Druyvesteyn (T_{eff} = %.2f eV)', results.T_eff_eV));
    legend('Location', 'best');
    grid on;
    
    % Plot 2: Relative difference
    subplot(1,3,2);
    if ~isempty(results.variable.f0)
        f0_var_interp = interp1(results.analytical.energy, results.analytical.f0, results.variable.energy, 'linear', 0);
        rel_diff_var = (f0_var_interp - results.variable.f0) ./ f0_var_interp * 100;
        idx_var = results.analytical.f0 > max(results.analytical.f0)*1e-6;
        plot(results.analytical.energy(idx_var), rel_diff_var(idx_var), 'b-', 'LineWidth', 2, 'DisplayName', 'Variable');
        hold on;
    end
    if ~isempty(results.uniform.f0)
        f0_uni_interp = interp1(results.analytical.energy, results.analytical.f0, results.uniform.energy, 'linear', 0);
        rel_diff_uni = (f0_uni_interp - results.uniform.f0) ./ f0_uni_interp * 100;
        idx_uni = results.analytical.f0 > max(results.analytical.f0)*1e-6;
        plot(results.analytical.energy(idx_uni), rel_diff_uni(idx_uni), 'r--', 'LineWidth', 2, 'DisplayName', 'Uniform');
    end
    xlabel('Energy (eV)');
    ylabel('Relative Error (%)');
    title('Error vs Analytical Solution');
    legend('Location', 'best');
    grid on;
    yline(0, 'k:', 'LineWidth', 1);
    
    % Plot 3: Error summary
    subplot(1,3,3);
    errors = [];
    labels = {};
    if ~isempty(results.variable.error) && ~isnan(results.variable.error)
        errors = [errors, results.variable.error*100];
        labels = [labels, {'Variable'}];
    end
    if ~isempty(results.uniform.error) && ~isnan(results.uniform.error)
        errors = [errors, results.uniform.error*100];
        labels = [labels, {'Uniform'}];
    end
    if ~isempty(errors)
        bar(errors);
        set(gca, 'XTickLabel', labels);
        ylabel('Mean Relative Error (%)');
        title('Error Summary');
        grid on;
    end
    
    % Salvar figura
    saveas(gcf, fullfile(figures_dir, 'maxwellian_const_v_comparison.png'));
    savefig(gcf, fullfile(figures_dir, 'maxwellian_const_v_comparison.fig'));
    fprintf('  Figura salva: maxwellian_const_v_comparison.png\n');
end

function generate_report(output_base, results_grid, results_max_elastic, results_max_ee, results_druy)
    % Gera relatório em texto
    
    report_file = fullfile(output_base, 'benchmark_report.txt');
    fid = fopen(report_file, 'w');
    
    fprintf(fid, '=================================================================\n');
    fprintf(fid, '          RELATÓRIO DE BENCHMARKING - LoKI-B                    \n');
    fprintf(fid, '=================================================================\n');
    fprintf(fid, 'Data: %s\n\n', datestr(now));
    
    % Resultados de comparação de malhas
    fprintf(fid, '--- COMPARAÇÃO MALHA VARIÁVEL VS UNIFORME ---\n\n');
    if ~isempty(results_grid)
        fprintf(fid, 'Testes com delta_u fixo:\n');
        fprintf(fid, '  Variável: %d testes\n', length(results_grid.fixed_delta_u.variable));
        fprintf(fid, '  Uniforme: %d testes\n\n', length(results_grid.fixed_delta_u.uniform));
        
        fprintf(fid, 'Testes com N fixo:\n');
        fprintf(fid, '  Variável: %d testes\n', length(results_grid.fixed_N.variable));
        fprintf(fid, '  Uniforme: %d testes\n\n', length(results_grid.fixed_N.uniform));
    else
        fprintf(fid, '  Sem resultados disponíveis.\n\n');
    end
    
    % Resultados analíticos
    fprintf(fid, '--- COMPARAÇÃO COM SOLUÇÕES ANALÍTICAS ---\n\n');
    
    if ~isempty(results_max_elastic)
        fprintf(fid, 'Maxwelliana com colisões elásticas:\n');
        fprintf(fid, '  Temperatura do gás: %.4f eV (%.1f K)\n', results_max_elastic.T_gas_eV, 300);
        fprintf(fid, '  Erro relativo (Variável): %.2f%%\n', results_max_elastic.variable.error*100);
        fprintf(fid, '  Erro relativo (Uniforme): %.2f%%\n\n', results_max_elastic.uniform.error*100);
    end
    
    if ~isempty(results_max_ee)
        fprintf(fid, 'Maxwelliana com colisões e-e:\n');
        fprintf(fid, '  Temperatura eletrônica: %.2f eV\n', results_max_ee.T_e_eV);
        fprintf(fid, '  Erro relativo (Variável): %.2f%%\n', results_max_ee.variable.error*100);
        fprintf(fid, '  Erro relativo (Uniforme): %.2f%%\n\n', results_max_ee.uniform.error*100);
    end
    
    if ~isempty(results_maxwellian_const_v)
        fprintf(fid, 'Maxwellian_const_v:\n');
      %  fprintf(fid, '  Temperatura efetiva: %.2f eV\n', results_maxwellian_const_v.T_eff_eV);
        fprintf(fid, '  Erro relativo (Variável): %.2f%%\n', results_maxwellian_const_v.variable.error*100);
        fprintf(fid, '  Erro relativo (Uniforme): %.2f%%\n\n', results_maxwellian_const_v.uniform.error*100);
    end
    
    fprintf(fid, '=================================================================\n');
    fprintf(fid, 'Fim do relatório\n');
    
    fclose(fid);
end

