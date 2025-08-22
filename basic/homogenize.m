function [Ez, nuzx] = homogenize(E1, E2, nu1, nu2, h1, h2)
    % Compute volume fractions
    f1 = h1 / (h1 + h2);
    f2 = 1 - f1;

    % Compute E_z (iso-stress modulus)
    Ez_num = E1 * E2 * (f1 * (1 - nu2) * E1 + f2 * (1 - nu1) * E2);
    Ez_den = E1 * E2 * (f1^2 * (1 - nu2) + f2^2 * (1 - nu1)) + ...
             f1 * f2 * ((1 + nu2) * (1 - 2 * nu2) * E1^2 + ...
                        4 * nu1 * nu2 * E1 * E2 + ...
                        (1 + nu1) * (1 - 2 * nu1) * E2^2);
    Ez = Ez_num / Ez_den;

    % Compute nu_zx (iso-stress Poisson's ratio)
    nuzx_num = E1 * E2 * (f1 * (1 - nu2) * nu1 + f2 * (1 - nu1) * nu2);
    nuzx_den = Ez_den; % Same denominator as Ez
    nuzx = nuzx_num / nuzx_den;

end
