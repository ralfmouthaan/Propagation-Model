function F = propLens(F, x, lambda, f)
    % PROPLENS Propagation through a lens
    % Author : RPM
    % Based on Eq. 5-10 of Goodman
    % Inputs:
    %   - F : fields in source plane
    %   - x : coordinates in source plane
    %   - lambda : wavelength
    %   - f : lens focal length
    % Outputs:
    %   - F = fields in observation plane

    k = 2*pi/lambda;
    [x_mesh, y_mesh] = meshgrid(x, x);
    F = F .* exp(-1i*k/2/f*(x_mesh.^2 + y_mesh.^2));
    
end