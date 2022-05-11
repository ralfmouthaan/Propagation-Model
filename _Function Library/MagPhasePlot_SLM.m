function MagPhasePlot_SLM(F, x)

    x = x*1e3;

    figure('Position', [500 500 250 250]);
    imagesc(x, x, abs(F));
    axis square
    colormap jet
    colorbar
    title('Magnitude');
    xlabel('mm'); ylabel('mm')
    xlim([-1 1]); ylim([-1 1])

    figure('Position', [500 500 250 250]);
    imagesc(x, x, mod(angle(F), 2*pi)-pi);
    axis square
    colormap gray;
    colorbar
    caxis([-pi pi])
    title('Phase')
    xlabel('mm'); ylabel('mm');
    xlim([-1 1]); ylim([-1 1]);

end