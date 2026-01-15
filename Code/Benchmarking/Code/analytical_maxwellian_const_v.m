function [f0, energy_eV] = analytical_maxwellian_const_v(energy_range, n_points)
% ANALYTICAL_MAXWELLIAN_CONST_V Calcula a EEDF Maxwelliana analítica
%   (obtida quando frequência de colisão é constante, nu = const)
%
% Inputs:
%   energy_range- [min_energy, max_energy] em eV
%   n_points    - Número de pontos
%
% Outputs:
%   f0          - EEDF Maxwelliana normalizada (eV^(-3/2))
%   energy_eV   - Array de energia (eV)
%
% A distribuição Maxwelliana é obtida quando:
%   - Há campo elétrico E/N não nulo
%   - Colisões elásticas com frequência constante (nu = const, sigma ∝ 1/v)
%   - Sem colisões inelásticas ou e-e
%
% A EEDF Maxwelliana é dada por:
%   f0(u) = C * exp(-u / T_eff)
%
% onde:
%   T_eff é a temperatura efetiva relacionada com E/N e propriedades do gás

    % if nargin < 2 || isempty(energy_range)
    %     energy_range = [0, 3*T_eff_eV];  % Druyvesteyn decai mais rápido
    % end
    % if nargin < 3 || isempty(n_points)
    %     n_points = 500;
    % end
    
    % % Validação de entrada
    % if T_eff_eV <= 0
    %     error('Temperatura efetiva deve ser positiva');
    % end
    % if energy_range(1) < 0
    %     energy_range(1) = 0;
    % end
    % if energy_range(2) <= energy_range(1)
    %     error('energy_range inválido: max deve ser > min');
    % end
    
    % Criar array de energia
    energy_eV = linspace(energy_range(1), energy_range(2), n_points);
    
    % Calcular EEDF Druyvesteyn
    % f0(u) = C * sqrt(u) * exp(-B * u^2)
    
    % Parâmetro B (ajustado para normalização adequada)

    N = 133.32 / 1.380649e-23 / 300;
    T_eff_eV = 1/3 * 1.602176634e-19 * (2e-3/6e23) * (100e-21/(1e-13*9.10938e-31))^2;
    
    % Calcular distribuição (sem normalização ainda)
    % f0_unnorm = zeros(size(energy_eV));
    % idx = energy_eV > 0;
    f0_unnorm = exp(-energy_eV/T_eff_eV);
    
    % Constante de normalização
    % Integral: int_0^inf sqrt(u) * exp(-B*u^2) du = (1/4) * sqrt(pi/B^3)
    % Para normalizar como EEDF: C = (2/sqrt(pi)) * B^(3/4)

    C = dot(f0_unnorm, sqrt(energy_eV))*energy_range(2)/n_points;
    %C = (2/sqrt(pi)) * (B^(3/4));
    
    % EEDF normalizada
    f0 = f0_unnorm/C;
    
    % Para u=0, f0(0) = 0 (já definido)
    
end

% % Função auxiliar para calcular T_eff a partir de E/N
% function T_eff = calculate_Teff_from_EN(E_N_Td, mass_gas_amu, nu, T_gas_K)
%     % Calcula temperatura efetiva Druyvesteyn a partir de E/N
%     %
%     % Inputs:
%     %   E_N_Td      - Campo reduzido em Townsend (1 Td = 1e-21 V·m^2)
%     %   mass_gas_amu- Massa do gás em u.m.a.
%     %   nu          - Frequência de colisão (s^-1)
%     %   T_gas_K     - Temperatura do gás (K)
%     %
%     % Output:
%     %   T_eff       - Temperatura efetiva em eV
% 
%     % Constantes
%     Td_to_SI = 1e-21;           % 1 Td = 1e-21 V·m^2
%     e = 1.602176634e-19;        % Carga elementar (C)
%     me = 9.10938e-31;           % Massa do eletrão (kg)
%     k_B = 1.380649e-23;         % Constante de Boltzmann (J/K)
%     amu_to_kg = 1.66053906660e-27; % u.m.a. para kg
% 
%     % Converter E/N para SI
%     E_N_SI = E_N_Td * Td_to_SI;  % V·m^2
% 
%     % Para Druyvesteyn com nu = const:
%     % T_eff ≈ (e^2 * E^2) / (3 * m * nu^2)
%     % Mas E depende de N, então usar E/N:
%     % T_eff ≈ (e * E/N)^2 / (3 * m/M * k_B * T_gas)
% 
%     % Aproximação simplificada (mais rigorosa requer solução da equação de Boltzmann)
%     M_gas = mass_gas_amu * amu_to_kg;
% 
%     % Temperatura efetiva (eV)
%     T_eff = (e * E_N_SI)^2 / (3 * (me/M_gas) * k_B * T_gas_K * e);
% 
% end
% 
