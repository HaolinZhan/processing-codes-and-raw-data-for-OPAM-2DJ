function signals = addWindow(signal, sw, param, option)
% Add window function for FID
if nargin ~= 4
    error('Wrong input...');
end
signal = reshape(signal, length(signal), 1);
M = length(signal);
k = (1:M).';
dt = 1 / sw;
option = upper(option);

switch option
    case {'EM', 'LB'}% exponential multiplication
        if length(param) == 1
            W = param;
        else
            error('wrong input')
        end
        win_EM = exp(-pi * W * k * dt);
        signals = signal .* win_EM;
    case 'GM'%Gaussian multiplication
        if length(param) == 1
            W = param;
        else
            error('wrong input')
        end
        win_GM = exp(-W * (k * dt).^2);
        signals = signal .* win_GM;
    case 'SGM'%Shifted gaussian multiplication
        if length(param) == 2
            W = param(1);
            P = param(2);
        else
            error('wrong input')
        end
        win_SGM = exp(-W * (k*dt-P*M*dt).^2);
        signals = signal .* win_SGM;
    case 'DM'%Double exponential multiplication(lorentz-to-gauss transforamtion)
        if length(param) == 2
            WE = param(1);
            WG = param(2);
        else
            error('wrong input')
        end
        win_DM = exp(pi * WE * k * dt) .* exp(-WG * (k*dt).^2);
        signals = signal .* win_DM;
    case 'SB'%Sine bell
        if length(param) == 1
            P = param;
        else
            error('wrong input')
        end
        win_SB = sin((pi-P)/M * k + P);
        signals = signal .* win_SB;
end