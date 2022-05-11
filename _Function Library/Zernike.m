% Ralf Mouthaan
% University of Cambridge
% December 2020
%
% Function to return zernike

function [Z] = Zernike(m, n, rho, phi)

    if size(rho) ~= size(phi)
        error('\phi and \rho must be of the same size');
    end
    
    if n < abs(m)
        error('n must be >= m')
    end
    
    if mod(n-m,2) == 1 && n ~= 0 
        error('This is equivalent to Z_{0,0}')
    end
    
    %%
    
    Z.m = abs(m);
    Z.n = n;
    Z.rho = rho;
    Z.phi = phi;
    
    R = zeros(size(rho));
    
    if mod(Z.n-Z.m,2) == 0
        for k = 0:(Z.n-Z.m)/2 % What if this is not an integer?

            term1 = (-1)^k;
            term2 = factorial(Z.n-k);
            term3 = factorial(k);
            term4 = factorial((Z.n+Z.m)/2-k);
            term5 = factorial((Z.n-Z.m)/2-k);
            term6 = rho.^(Z.n-2*k);

            R = R + term1*term2/term3/term4/term5*term6;

        end
    end
    
    if m >= 0
        Z.F = R.*cos(Z.m*phi);
    else
        Z.F = R.*sin(Z.m*phi);
    end

end