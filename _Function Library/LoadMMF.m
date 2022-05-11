function Fibre = LoadMMF (f, lambda, TargetFile)
    % LOADMMF Load Fibre data
    % Author: RPM
    % Inputs:
    %   f : Focal length of lens that focuses light onto facet
    %   lambda : wavelength
    %   TargetFile : mat file with hologram target data as created from
    %           genTarget.m
    %   ModeNo : Index of target to be loaded from file
    % Outputs:
    %   Fibre.F : Field distribution in fibre
    %   Fibre.L : L index for given mode
    %   Fibre.M : M index for given mode
    %   Fibre.O : O index for given mode
    %   Fibre.beta : mode propagation constant
    %   Fibre.x : x coordinates along one axis of target
    
    % Note: JC's code does not output dimensions, so that is done here
    % based on lambda and f. Ideally, the dimensions would just be given in
    % the mat file.
    
    % Determine pixel pitch based on whether it's a Hamamatsu or HoloEye.
    if contains(TargetFile, 'Ham')
        PixelPitch = 20e-6;
    elseif contains(TargetFile, 'Holo')
        PixelPitch = 8e-6;
    elseif contains(lower(TargetFile), 'tx')
        PixelPitch = 20e-6; % NB, transmit may not always be hamamatsu
    elseif contains(lower(TargetFile), 'rx')
        PixelPitch = 8e-6; % NB, receive may not always be holoeye
    end

    % Load from file
    Fibre = load(TargetFile, 'TARGET', 'L', 'M', 'O', 'beta');
    Fibre.F = Fibre.TARGET;
    Fibre = rmfield(Fibre, 'TARGET');
    
    % Calculate target dimensions based on JC's heuristic (see genTarget.m)
    Fibre.x = f*lambda / PixelPitch / size(Fibre.F, 2) ...
        * (-(size(Fibre.F, 2)-1)/2:(size(Fibre.F,2)-1)/2);
    
    Fibre.NA = 0.2;

end