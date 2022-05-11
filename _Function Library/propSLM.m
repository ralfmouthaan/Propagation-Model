function [F, SLMPattern] = propSLM(F, x, lambda, SLM)
    % PROPSLM Propagation through SLM
    % Author: RPM
    % Inputs:
    %   - F : field in source plane
    %   - x : coordinates in source plane
    %   - lambda : wavelength
    %   - SLM : SLM structure
    % Outputs :
    %   - F : fields in observation plane
    % SLM struct : 
    %   - SLM.Holo : Displayed hologram without adjustments like tilt/focus
    %   - SLM.Holox : Coordinates in source plane
    %   - SLM.Tiltx : x tilt
    %   - SLM.Tilty : y tilt
    %   - SLM.FocalLength : Focal length
    
    if (x(2) - x(1)) > (SLM.x(2) - SLM.x(1))
        warning('Field is undersampling SLM')
    end
    
    k = 2*pi/lambda;
    [x_mesh, y_mesh] = meshgrid(x + SLM.Offsetx, x + SLM.Offsety);
    SLMPattern = SLM.Holo;
    SLMPattern = interp2(SLM.x, SLM.x', SLMPattern, x, x', 'nearest');
    SLMPattern(isnan(SLMPattern)) = 1;
    SLMPattern = SLMPattern.*exp(-1i*k/2*(x_mesh*SLM.Tiltx + y_mesh*SLM.Tilty));
    SLMPattern = SLMPattern.*exp(-1i*k/2*SLM.Focus*(x_mesh.^2 + y_mesh.^2));
    F = F.*SLMPattern;

end