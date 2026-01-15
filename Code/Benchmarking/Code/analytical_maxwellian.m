function [f0, energy_eV] = analytical_maxwellian(T_eV, energy_range, n_points)
% ANALYTICAL_MAXWELLIAN Calcula a EEDF Maxwelliana analítica
%
% Inputs:
%   T_eV        - Temperatura em eV
%   energy_range- [min_energy, max_energy] em eV (opcional, default [0, 5*T_eV])
%   n_points    - Número de pontos (opcional, default 500)
%
% Outputs:
%   f0          - EEDF Maxwelliana normalizada (eV^(-3/2))
%   energy_eV   - Array de energia (eV)
%
% A distribuição Maxwelliana é dada por:
%   f0(u) = (2/sqrt(pi)) * (1/T_e)^(3/2) * sqrt(u) * exp(-u/T_e)
%
% onde u é a energia em eV e T_e é a temperatura em eV

    if nargin < 2 || isempty(energy_range)
        energy_range = [0, 5*T_eV];
    end
    if nargin < 3 || isempty(n_points)
        n_points = 500;
    end
    
    % Validação de entrada
    if T_eV <= 0
        error('Temperatura deve ser positiva');
    end
    if energy_range(1) < 0
        energy_range(1) = 0;
    end
    if energy_range(2) <= energy_range(1)
        error('energy_range inválido: max deve ser > min');
    end
    
    % Criar array de energia
    energy_eV = linspace(energy_range(1), energy_range(2), n_points);
    
    % Calcular EEDF Maxwelliana
    % f0(u) = (2/sqrt(pi)) * (1/T_e)^(3/2) * sqrt(u) * exp(-u/T_e)
    
    % Evitar divisão por zero em u=0
    %f0 = zeros(size(energy_eV));
    %idx = energy_eV > 0;
    
    prefactor = (2/sqrt(pi)) * (1/T_eV)^(3/2);
    f0 = prefactor * exp(-energy_eV/T_eV);

    C = dot(f0, sqrt(energy_eV))*energy_range(2)/n_points;

    f0 = f0 / C;

    C = dot(f0, sqrt(energy_eV))*energy_range(2)/n_points;

    
    % Para u=0, f0(0) = 0 (já definido)
    
    % Normalização (verificação)
    % Integral de f0 * 4*pi*sqrt(2*m*u)/h^3 sobre u deve ser n
    % Para eV^(-3/2), a normalização já está incluída no prefactor
    
end

% Função auxiliar para temperatura em Kelvin
function T_eV = kelvin_to_eV(T_K)
    % Converte temperatura de Kelvin para eV
    % T_eV = k_B * T_K / e
    k_B = 1.380649e-23;  % Constante de Boltzmann (J/K)
    e = 1.602176634e-19; % Carga elementar (C)
    T_eV = k_B * T_K / e;
end

