function [rateCoeff, dependent] = oxygenMolecularVT(~, ~, ~, reaction, ~, chemistry)

% oxygenMolecularVT evaluates the rate coefficient for V-T transitions between oxygen molecules
% described in:
% Adriana Annusova et al 2018 Plasma Sources Sci. Technol. 27 045006
% https://doi.org/10.1088/1361-6595/aab47d
% including also dissociation due to VT transitions (coefficients extrapolated from VTs, by JSilva (2019))
% The rate coefficients are to be used at gas temperatures in the interval 300-1000 K
% The use of this function is RESTRICTED to a set of O2(X,v=0-41) levels
%  - the tabulated coefficients are for vibrational levels in this range
%  - the function assumes level O2(X,v=42) as the dissociation continuum
% Last changes and checks by LLAlves (June 2020)

    % local copies and definitions of different parameters used in the function
    Tg = chemistry.workCond.gasTemperature; %Tg in K
    KbTg = Constant.boltzmannInEV*Tg; %KbTg in eV
    omega = 2*pi*299792458*1580.161*100; %harmonicOscillatorFrequency in rad/s
  
    persistent coeffM; %coefficients' matrix for v >= 33
    persistent coeffA; %coefficients' array for dissociation
   
    % identify current V-T reaction by evaluating the vibrational levels (O2(X,v) + O2(X) <-> O2(X,v-1) + O2(X)) 
    % check first if this is a dissociation reaction (O2(X,41) + O2(X) -> 2O(3P) + O2(X))

 
    if sum(reaction.productStoiCoeff) > 1

        w = str2double(reaction.reactantArray(1).vibLevel);
        if isempty(coeffA) %extrapolated coefficients for v=42 (reactions are written as deexcitations (w = v-1 = 41))
            % These extrapolated coefficients WILL NOT WORK AT VERY-HIGH GAS TEMPERATURES
            coeffA = [-4.232558e-21 1.7724117e-17 -6.912383e-15 3.9297334e-11];
        end
        a = coeffA(1);
        b = coeffA(2);
        c = coeffA(3);
        d = coeffA(4);

        directRateCoeff = 0.25*(a*Tg^3 + b*Tg^2 + c*Tg + d)*1e-6;
% generalization of the last level leading to dissociation (LLA, June 2020)
%         exponent = Constant.planckReducedInEV*omega*(w-42);
        exponent = Constant.planckReducedInEV*omega*(-1);
        rateCoeff = directRateCoeff*exp(exponent/KbTg);

    else % Molecular VT reactions not leading to dissociation
        
        v = str2double(reaction.reactantArray(1).vibLevel);
        if v <= 6 && v >= 1
            %1st part of Table A2
            rateCoeff = 0.25*exp((-2.16870e-3*v^3 + 1.02950e-2*v^2-9.03306e-2*v + 6.02557)*log(Tg)+...
                (2.7813e-2*v^3 - 2.76608e-1*v^2 + 2.04169*v - 76.0376))*1e-6;
        elseif v >= 7 && v <= 32
            %1st part of Table A2 (cont)
            rateCoeff = 0.25*(2.44106e-33*Tg^6.25015)*exp((1.76911e-7*Tg^2 - 4.82968e-4*Tg + ...
                5.40946e-1)*v)*1e-6;
        elseif v<= 41 && v >= 33
            if isempty(coeffM)
                %2nd part of Table A2
                coeffM = [8.11526e-20 -1.62934e-16 1.34668e-13 -2.19694e-11
                6.59091e-20 -1.40487e-16 1.27056e-13 -1.76905e-11
                7.60101e-20 -1.66732e-16 1.43746e-13 -1.6131e-11
                7.57187e-20 -1.61678e-16 1.32137e-13 -7.15596e-12
                3.60748e-20 -7.34407e-17 6.7491e-14  1.04153e-11
                3.40294e-20 -6.43874e-17 5.61750e-14 1.57281e-11
                1.76281e-20 -3.04842e-17 3.28204e-14 2.21318e-11
                1.11133e-20 -1.56706e-17 2.12930e-14 2.56952e-11
                1.24648e-20 -1.79119e-17 2.19332e-14 2.61423e-11];
            end

            index = v - 32; %vib level starts at 33 here
            a = coeffM(index, 1);
            b = coeffM(index, 2);
            c = coeffM(index, 3);
            d = coeffM(index, 4);

            rateCoeff = 0.25*(a*Tg^3 + b*Tg^2 + c*Tg + d)*1e-6;

        else
            error('Invalid vibrational level')
        end

    end

    % set function dependencies
    dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end
