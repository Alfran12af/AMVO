function plot_field( u )


    figure
    imagesc(u')
    set(gca, 'YDir', 'normal');

    colorbar;
    colormap('jet');
    axis equal tight;
 
end