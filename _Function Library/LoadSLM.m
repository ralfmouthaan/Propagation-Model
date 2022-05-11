function SLM = LoadSLM(HologramFile, ModeNo)
    % LOADSLM Load SLM data
    %   Hologram : mat file with hologram data as created by C++
    %           HologramGenerator code
    %   ModeNo : Index of hologram to be loaded from file
    % Outputs:
    %   SLM struct:
    %       SLM.Holo : matrix containing hologram displayed on SLM
    %       SLM.L : L subscript for displayed hologram
    %       SLM.M : M subscript for displayed hologram
    %       SLM.O : O subscript for displayed hologram
    %       SLM.PixelPitch : Width of pixels on SLM
    %       SLM.x : x coordinates along one axis of hologram
    
    % Load from file, tidy up a bit.
    SLM = load(HologramFile, 'ILLUMINATION', 'L', 'M', 'O');
    SLM.Holo = squeeze(SLM.ILLUMINATION(ModeNo,:,:));
    SLM.L = SLM.L(ModeNo);
    SLM.M = SLM.M(ModeNo);
    SLM.O = SLM.O(ModeNo);
    SLM = rmfield(SLM, 'ILLUMINATION');

    % Determine pixel pitch based on whether it's a Hamamatsu or HoloEye.
    if contains(HologramFile, 'Ham')
        SLM.PixelPitch = 20e-6;
    elseif contains(HologramFile, 'Hol')
        SLM.PixelPitch = 8e-6;
    end
    
    % Calculate dimensions from pixel pitch and size.
    SLM.x = SLM.PixelPitch * (-(size(SLM.Holo, 1)-1)/2 : ...
        (size(SLM.Holo, 1)-1)/2);
    
end