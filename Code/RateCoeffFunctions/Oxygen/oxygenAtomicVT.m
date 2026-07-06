function [rateCoeff, dependent] = oxygenAtomicVT(~, ~, ~, reaction, ~, chemistry)
% oxygenAtomicVT evaluates the rate coefficient for V-T transitions between oxygen molecules and atoms described in:
% Esposito et al 2008 Chemical Physics Volume 351 91-98
% https://doi.org/10.1016/j.chemphys.2008.04.004
% (note that one should refer to this paper to understand how the coefficients c_ijk are extracted from the matrices in the code)

% The rate coefficients are to be used at gas temperatures in the interval 50-10000 K
% The use of this function is RESTRICTED to a set of O2(X,v=0-41) levels
%  - the tabulated coefficients are for vibrational levels in this range
%  - the function assumes level O2(X,v=42) as the dissociation continuum
% Last changes and checks by LLAlves (June 2020)

% see also
% Adriana Annusova et al 2018 Plasma Sources Sci. Technol. 27 045006
% https://doi.org/10.1088/1361-6595/aab47d


    % local copies and definitions of different parameters used in the function
    Tg = chemistry.workCond.gasTemperature;
    KbTg = Constant.boltzmannInEV*Tg; %KbTg in eV
    omega = 2*pi*299792458*1580.161*100; %harmonicOscillatorFrequency in rad/s
    
    persistent cCoeffM1;
    persistent cCoeffM2;
    persistent cCoeffM3;
    persistent bCoeffM;
    persistent aCoeff;

    if isempty(cCoeffM1)
        %Table 2 (Esposito et al 2008): Coefficient c1jk for the calculation of a1 (v, Dv)
        %j=1:5, k=1:4; each line corresponds to a pair (j,k): (1,1), (1,2), ... (2,1),(2,2),...
        cCoeffM1 = [-26.18993227 -19.45675583 -3.993663 2.821759
        0 -3.38007590 2.684515 4.083138
        0 89.85158436 -1.009927e5 -8.809991e1
        0 5.853646e-2 -0.283676 2.369644e-2
        -1.69828917 -1.48712025 -3.030965 3.841080
        0 -0.54018334 -4.594443 -0.998637
        0 -53.33457116 3.590903e4 -3.438479e1
        0 -9.543789e-2 7.104764e-2 1.222381e-2
        3.349076e19 1.505136e21 5.492061e21 4.330479e21
        0 1.621622e21 1.212196e21 -1.677646e22
        0 -1.066183e21 9.092488e21 -5.573334e21
        0 2.169160e20 1.540038e20 -3.089812e19
        -3.946126e20 -1.532916e20 1.308503e20 -1.194045e20
        0 -4.105380e19 1.831856e20 -5.121704e19
        0 1.185042e21 1.079540e22 4.013656e21
        0 -1.748646e18 -5.608629e18 -1.052730e18
        1.391056e19 -4.838473e18 2.160753e19 -9.939380e18
        0 9.529399e18 -1.465914e19 4.180038e18
        0 5.290114e19 5.483520e21 -2.265448e18
        0 6.021495e16 1.142128e17 1.050411e16];
    end
    
    if isempty(cCoeffM2)
        %Table 3: Coefficient c2jk for the calculation of a2 (v, Dv)
        %j=1:5, k=1:4; each line corresponds to a pair (j,k): (1,1), (1,2), ... (2,1),(2,2),...
        cCoeffM2 = [7.83331061 14.47821980 -7.575157e1 -1.066105e2
        0 -43.32225364 -9.234026 6.618737e1
        0 -348.16803354 2.807916e5 2.630406e2
        0 -0.16418596 0.333397 2.791153e-2
        3.71221451 32.66821082 -7.713850 -2.658559e1
        0 -0.18462189 2.545634e1 -9.167211
        0 119.03129957 -9.592245e4 1.678357e2
        0 0.22258852 -0.179262 -0.106474
        3.573261e20 1.522507e21 4.002520e21 1.312884e22
        0 2.654567e21 -8.192010e21 -5.437653e21
        0 -3.528669e21 1.462011e23 5.735816e21
        0 -2.861293e20 -1.821224e20 1.568233e20
        6.433503e20 -1.533872e21 -2.912948e21 4.530951e21
        0 -1.522587e20 6.399791e19 2.662341e20
        0 -3.124046e21 -3.531505e22 -2.932068e22
        0 2.322209e19 1.964744e19 6.788371e18
        -2.901352e19 -3.762650e19 -7.070723e19 -3.473472e19
        0 1.955942e19 4.805948e19 -5.623449e18
        0 8.847719e18 5.201014e21 2.765213e20
        0 -8.252347e17 -4.170396e17 -6.030509e16];
    end
    
    if isempty(cCoeffM3)
        %Table 4: Coefficient c3jk for the calculation of a3 (v, Dv)
        %j=1:5, k=1:4; each line corresponds to a pair (j,k): (1,1), (1,2), ... (2,1),(2,2),...
        cCoeffM3 = [0.37163948 0.88657647 -1.271181 -3.476825
        0 -1.14258799 0.317834 -0.415675
        0 -4.53080420 5.830886e3 2.341590e1
        0 4.732032e-3 7.186328e-3 -1.760866e-3
        0.10587091 -0.18566919 -0.709028 -0.562411
        0 -0.19256423 4.753706e-2 0.166319
        0 3.95819117 -1.757570e3 2.356659
        0 1.353007e-2 1.465161e-3 -7.409818e-4
        -5.312491e19 -2.027215e21 -8.686388e20 1.908092e21
        0 2.381051e21 9.428176e20 1.107010e21
        0 -1.248596e20 1.738719e22 -1.769244e22
        0 -4.014395e19 -1.251868e19 1.272578e18
        3.754092e19 1.819921e19 1.877581e19 -2.139241e19
        0 3.708631e19 1.097908e19 2.171483e19
        0 2.031805e19 -1.006633e22 -2.478535e19
        0 -6.021031e17 2.734279e16 7.612186e16
        -1.189832e18 1.530423e18 2.238756e18 1.542813e18
        0 -1.859900e18 -7.688157e17 -7.531694e17
        0 -8.207403e18 -2.061752e21 -4.924709e18
        0 9.708459e15 -3.427940e15 -5.546602e14]; 
    end
    
    % identify current V-T reaction by evaluating the vibrational levels (O2(X,v) + O <-> O2(X,w) + O ; w < v)
    % check first if this is a dissociation reaction (O2(X,41) + O -> 2O(3P) + O)
    
    if sum(reaction.productStoiCoeff) > 1
        v = 42; % reactions are written as deexcitations (w < v = 42)
        w = str2double(reaction.reactantArray(1).vibLevel);
    else
        v = str2double(reaction.reactantArray(1).vibLevel);  
        w = str2double(reaction.productArray(1).vibLevel);
    end
    
    deltav = v-w;
    
    if deltav == 1
        columnIndex = 1;
    elseif deltav >= 2 && deltav <= 10
        columnIndex = 2;
    elseif deltav >= 11 && deltav <= 20
        columnIndex = 3;
    elseif deltav >= 21 && deltav <= 30
        columnIndex = 4;
    else
        error('Quanta jump not allowed');
    end
    
    if isempty(bCoeffM)
        bCoeffM = zeros(3,5);
    end
    
    for j=1:5
        bCoeffM(1,j) = cCoeffM1(4*(j-1)+1, columnIndex)+...
            cCoeffM1(4*(j-1)+2, columnIndex)*log(deltav)+...
            cCoeffM1(4*(j-1)+3, columnIndex)*deltav*exp(-deltav)+...
            cCoeffM1(4*(j-1)+4, columnIndex)*deltav^2;
                
        bCoeffM(2,j) = cCoeffM2(4*(j-1)+1, columnIndex)+...
            cCoeffM2(4*(j-1)+2, columnIndex)*log(deltav)+...
            cCoeffM2(4*(j-1)+3, columnIndex)*deltav*exp(-deltav)+...
            cCoeffM2(4*(j-1)+4, columnIndex)*deltav^2;
                
        bCoeffM(3,j) = cCoeffM3(4*(j-1)+1, columnIndex)+...
            cCoeffM3(4*(j-1)+2, columnIndex)*log(deltav)+...
            cCoeffM3(4*(j-1)+3, columnIndex)*deltav*exp(-deltav)+...
            cCoeffM3(4*(j-1)+4, columnIndex)*deltav^2;
    end
       
    if isempty(aCoeff)
        aCoeff = zeros(1,3);
    end
    
    for i=1:3
        aCoeff(i) = bCoeffM(i,1) + bCoeffM(i,2)*log(v)+...
            (bCoeffM(i,3)+bCoeffM(i,4)*v+bCoeffM(i,5)*v^2)/...
            (10^21+exp(v));  
    end
    
    degF =1/((3*(5+3*exp(-227/Tg)+exp(-325.9/Tg)))); % as done in LoKI1.2.0
    directRateCoeff = degF*exp(aCoeff(1) + aCoeff(2)/log(Tg) + aCoeff(3)*log(Tg))*1e-6;
    
    if v ~= 42
        rateCoeff = directRateCoeff;
    else
        exponent = Constant.planckReducedInEV*omega*(w-42);
        rateCoeff = directRateCoeff*exp(exponent/KbTg);
    end
    
    % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);
end