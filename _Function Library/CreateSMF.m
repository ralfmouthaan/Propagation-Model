function SMF = CreateSMF(x)
    % CREATESMF Create SMF structure
    % Author : RPM
    % Inputs:
    %   x : coordinates
    %   w0 : width of fiber
    % Outputs:
    %   SMF.x : coordinates
    %   SMF.F : fields at SMF facet
    %   SMF.w0 : fiber width

    % Coordinates
    SMF.x = x;
    SMF.w0 = 8.2e-6;

    % Gaussian beam from SMF
    [x_mesh, y_mesh] = meshgrid(x, x);
    SMF.F = exp(2*(-x_mesh.^2 - y_mesh.^2)/SMF.w0^2);
    SMF.F = SMF.F/max(max(SMF.F));   
    
    SMF.NA = 0.14;

end