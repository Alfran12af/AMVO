function animate_NS(x, y, U, V, P, t, rho)
% ANIMATE_NS - Post-Processing and Animation for 2D Navier-Stokes
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function performs the post-processing of the stored 2D velocity 
%   and pressure fields. It computes secondary physical variables such as 
%   vorticity or velocity magnitude and generates a dashboard animation.
%
% Inputs:
%   x, y    - Coordinate matrices for the NxN interior domain
%   U, V, P - 3D matrices (NxNxFrames) containing the temporal history
%   t       - Vector containing the physical time for each stored frame
%   rho     - Density of the fluid

 
    n_frames = size(U, 3);
    N = size(U, 1);

    % Auxiliar vectors
    x_vec = x(:,1); 
    y_vec = y(1,:);
    h = x_vec(2) - x_vec(1);

    % Pre-locating variables
    vorticity = zeros(N, N, n_frames);
    velocity  = zeros(N, N, n_frames);
    K_energy  = zeros(N, N, n_frames);
    total_K   = zeros(1, n_frames);
    
    for k = 1:n_frames
        % Velocity module
        velocity(:,:,k) = sqrt(U(:,:,k).^2 + V(:,:,k).^2);
        
        % Vorticity (2D Rotational)
        [dvdx, ~] = gradient(V(:,:,k), h);
        [~, dudy] = gradient(U(:,:,k), h);
        vorticity(:,:,k) = dvdx - dudy;

        % Kinetic energy
        K_energy(:,:,k) = 0.5 * rho*(U(:,:,k).^2 + V(:,:,k).^2);
        total_K(k) = sum(sum(K_energy(:,:,k))) * h^2;
    end
    
    % Limit defition
    f_sat = 0.7; % Saturation factor
    lim_vort = [min(vorticity(:)) max(vorticity(:))] * f_sat;
    lim_u    = [min(U(:)) max(U(:))] * f_sat;
    lim_v    = [min(V(:)) max(V(:))] * f_sat;
    lim_p    = [min(P(:)) max(P(:))] * f_sat;
    lim_k    = [0, max(K_energy(:)) * 0.6];

    % Lagrangian particles grid set up
    [xp_grid, yp_grid] = meshgrid(linspace(0.1, 0.9, 15), linspace(0.1, 0.9, 15));
    xp = xp_grid(:); 
    yp = yp_grid(:);

    % Colormap
    c1 = [0, 0, 1];
    c2 = [1, 1, 1];
    c3 = [1, 0, 0];
    custom_bwr = [linspace(c1(1),c2(1),128)', linspace(c1(2),c2(2),128)', linspace(c1(3),c2(3),128)';
                  linspace(c2(1),c3(1),128)', linspace(c2(2),c3(2),128)', linspace(c2(3),c3(3),128)'];

    % Video set up
    v = VideoWriter('NavierStokes_Result.mp4', 'MPEG-4');
    v.FrameRate = 8;
    v.Quality = 100;
    open(v);

    % Visualization
    h_fig = figure('Color', 'w', 'Name', 'NS 2D: Visualization');
    set(h_fig, 'Position', [50, 50, 1400, 900]);

    tic

    for k = 1:n_frames

        if ~ishandle(h_fig), break; end
        clf(h_fig);

        if k > 1 
            dt_frame = t(k) - t(k-1); 
        else
            dt_frame = 0; 
        end



        % ---------------------------------------------------------
        % PANEL 1: Vorticity + streamlines
        % ---------------------------------------------------------
        subplot(2,3,1);
        contourf(x_vec, y_vec, vorticity(:,:,k), 50, 'LineColor', 'none');
        hold on;
        h_s = streamslice(x_vec, y_vec, U(:,:,k), V(:,:,k), 2);
        set(h_s, 'Color', [0.1 0.1 0.1, 0.5], 'LineWidth', 0.5);
        colormap(gca, custom_bwr); % jet(256)
        clim(lim_vort); 
        colorbar; 
        axis equal tight;
        title('Vorticity and streamlines', 'FontSize', 11);
        ylabel('y [m]'); xlabel('x [m]');



        % ---------------------------------------------------------
        % PANEL 2: Pressure + velocity vectors (Quiver)
        % ---------------------------------------------------------
        subplot(2,3,2);
        p = pcolor(x_vec, y_vec, P(:,:,k));
        set(p, 'EdgeColor', 'none', 'FaceColor', 'interp');
        hold on;

        % Vectors: quiver
        skip = max(1, floor(N/12)); 
        quiver(x(1:skip:end, 1:skip:end), y(1:skip:end, 1:skip:end), ...
               U(1:skip:end, 1:skip:end, k), V(1:skip:end, 1:skip:end, k), ...
               1.5, 'w', 'LineWidth', 1);
        colormap(gca, parula(256));
        clim(lim_p); 
        colorbar; 
        axis equal tight;
        title('Preassure and velocity vectors', 'FontSize', 11);
        axis equal;
        xlim([0 1]);
        ylim([0 1]);
        ylabel('y [m]'); xlabel('x [m]');



        % ---------------------------------------------------------
        % PANEL 3: Velocity U
        % ---------------------------------------------------------
        subplot(2,3,3);
        contourf(x_vec, y_vec, U(:,:,k), 50, 'LineColor', 'none');
        colormap(gca, parula(256));
        clim(lim_u); 
        colorbar; 
        axis equal tight;
        title('Horizontal velocity (U)');
        ylabel('y [m]'); xlabel('x [m]');



        % ---------------------------------------------------------
        % PANEL 4: Velocity V
        % ---------------------------------------------------------
        subplot(2,3,4);
        contourf(x_vec, y_vec, V(:,:,k), 50, 'LineColor', 'none');
        colormap(gca, parula(256));
        clim(lim_v); 
        colorbar; 
        axis equal tight;
        title('Vertical velocity (V)');
        ylabel('y [m]'); xlabel('x [m]');



        % ---------------------------------------------------------
        % PANEL 5: Kinetic Energy Map
        % ---------------------------------------------------------
        subplot(2,3,5);
        contourf(x_vec, y_vec, K_energy(:,:,k), 50, 'LineColor', 'none');
        set(gca, 'YDir', 'normal');
        colormap(gca, 'hot');
        colorbar;
        clim(lim_k);
        axis equal tight;
        title('Kinetic energy density');
        ylabel('y [m]'); xlabel('x [m]');



        % ---------------------------------------------------------
        % PANEL 6: Lagrangian particle tracking
        % ---------------------------------------------------------
        subplot(2,3,6);
        if dt_frame > 0
            up = interp2(x_vec, y_vec, U(:,:,k), xp, yp, 'linear', 0);
            vp = interp2(x_vec, y_vec, V(:,:,k), xp, yp, 'linear', 0);
            xp = mod(xp + up * dt_frame, 1); 
            yp = mod(yp + vp * dt_frame, 1);
        end
        scatter(xp, yp, 15, 'filled', 'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', 'k');
        hold on;
        h_s = streamslice(x_vec, y_vec, U(:,:,k), V(:,:,k), 2);
        set(h_s, 'Color', [0.1 0.1 0.1, 0.1], 'LineWidth', 0.5);
        axis equal tight; 
        xlim([0 1]); 
        ylim([0 1]);
        title('Lagrangian particle tracking');        
        ylabel('y [m]'); xlabel('x [m]');

        drawnow;
        sgtitle(['Navier-Stokes simulation dashboard | t = ', num2str(t(k), '%.3f'), 's']);


        % Video
        frame = getframe(h_fig);
        writeVideo(v, frame);


    end

    toc
    close(v)

    % ---------------------------------------------------------
    % Figure: Energy decay
    % ---------------------------------------------------------
    figure('Color', 'w'); 
    plot(t, total_K, 'r-', 'LineWidth', 2.5);
    grid on;
    xlabel('Time [s]', 'FontSize', 12);
    ylabel('Total kinetic energy', 'FontSize', 12);
    title('Dissipation of total kinetic energy', 'FontSize', 14);
    xlim([0 t(end)]);
    ylim([0 max(total_K)*1.1]); 

end