function MagPhasePlot_Image(F, x)

        x = x*1e6;

        figure('Position', [500 500 250 250]);
        imagesc(x, x, abs(F));
        xlim([-50 50]); ylim([-50 50]);
        axis square
        colormap jet;
        colorbar
        xlabel('\mum'); ylabel('\mum');
        title('Magnitude')
        
        figure('Position', [500 500 250 250]);
        imagesc(x, x, angle(F));
        xlim([-50 50]); ylim([-50 50]);
        axis square
        colormap jet;
        colorbar
        title('Phase')
        xlabel('\mum'); ylabel('\mum');

end