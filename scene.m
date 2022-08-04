function scene
    
    %   Keyboard
    %   --------
    %
    %    Move    Turn
    %      e       i
    %    s d f   j k l
    
    %   Unit-Vec Schematic              Values
    %   ------------------              ------
    %            z         
    %            |    .                 u from 'o' to '.'
    %            |   /|        
    %            |  u |                 ux = cos(el)*cos(az)
    %            | /  |                 uy = cos(el)*sin(az)
    %            o/|el| ___ y           uz = sin(el)
    %           /_\|  |         
    %          /az \  |                 u = [ux,uy,uz]
    %         /     \ |        
    %        /       \|       
    %       x 
    
    % ------------------------- Camera -------------------------
    
    cam_pos = [15,15,5];        
    cam_el  = -.2;              
    cam_az  = -2;              
    
    cam_el_speed   = 0;           
    cam_az_speed   = 0;
    cam_fwd_speed  = 0;
    cam_side_speed = 0;
    
    unit_vec = @(el,az)[            
        cos(el)*cos(az), ...
        cos(el)*sin(az), ...
        sin(el)];
    
    % -------------------------- Gui -------------------------
    
    fig = figure( ...
        'menubar', 'none', ...
        'numbertitle','off', 'name','3d scene', ...
        'color', [.2, .2, .2], ...
        'keypressfcn', @key_down, ...
        'keyreleasefcn', @key_up, ...
        'closerequestfcn', @close_fcn);
    
    ax = axes( ...
        'projection', 'perspective', ...
        'cameraviewangle', 60, ...
        'dataaspectratio', [1,1,1], ...
        'visible', 'off');
     
    annotation('textbox', ...
        'position', [0,0,1,1], ...
        'edgecolor', 'none', ...
        'color', [.7, .8, .9], ...
        'fontsize', 16, ...
        'interpreter', 'latex', ...
        'string', ...
            {'. e . . . . i .';
             's d f . . j k l'});
    
    % ------------------------- Graphics -----------------------
    
    % patches
    for i = 1:100   
        patch(...
            'xdata', 10*rand*[1,1,1] + rand(1,3), ...
            'ydata', 10*rand*[1,1,1] + rand(1,3), ...
            'zdata', 10*rand*[1,1,1] + rand(1,3), ...
            'facecolor', rand(1,3));
    end
    
    % surface
    [ x_matrix, y_matrix ] = meshgrid(...
            linspace(0,20), ...
            linspace(0,20));
         
    plot_surf = surface(...
        'xdata', linspace(0,20), ...
        'ydata', linspace(0,20), ...
        'zdata', cos(x_matrix), ...
        'cdata', rand(100,100, 3), ...
        'facelighting', 'gouraud');
    
    % lights
    light;
    light('style', 'local', 'position', [5,5,5]);
    
    % animated line
    anim_line = animatedline(...
        'color', [.5, .6, .7], ...
        'linestyle', 'none', ...
        'marker', 'o');
   
    % lines 
    line([0, 0], [0, 0], [0,10], 'linewidth', 3);
    line([0, 0], [0,10], [0, 0], 'linewidth', 3);
    line([0,10], [0, 0], [0, 0], 'linewidth', 3);
    
    % text (within axes)
    text(11,0,0, 'x', 'fontsize', 30, 'color', [.5, .6, .7]);
    text(0,11,0, 'y', 'fontsize', 30, 'color', [.5, .6, .7]);
    text(0,0,11, 'z', 'fontsize', 30, 'color', [.5, .6, .7]);
    
    % ---------------------------- Loop --------------------------------
    
    loop_on = true;
    tic;
        
    while loop_on
        
        % Update graphics objects 
        anim_line.addpoints(...
            toc*cos(toc), sin(toc)*toc, toc);
        
        plot_surf.ZData = ...
            cos(x_matrix + toc).*sin(y_matrix + 2*toc);
        
        % Update camera
        cam_el = cam_el + cam_el_speed;
        cam_az = cam_az + cam_az_speed;
        
        cam_pos = cam_pos + ...
            cam_fwd_speed*unit_vec(cam_el, cam_az) + ...
            cam_side_speed*unit_vec( 0, cam_az + pi/2);
        
        ax.CameraPosition = cam_pos;
        ax.CameraTarget = cam_pos + unit_vec(cam_el, cam_az);
        
        drawnow; 
    end
    
    % ------------------------ Functions --------------------------
    
    function key_down(~,e)
        switch e.Key
            case 'e', cam_fwd_speed  =  .1;
            case 'd', cam_fwd_speed  = -.1;
            case 's', cam_side_speed =  .1;
            case 'f', cam_side_speed = -.1;
                
            case 'i', cam_el_speed =  .02; 
            case 'k', cam_el_speed = -.02; 
            case 'j', cam_az_speed =  .02;
            case 'l', cam_az_speed = -.02;
        end
    end
    function key_up(~,~)
        cam_el_speed   = 0;
        cam_az_speed   = 0;
        cam_fwd_speed  = 0;
        cam_side_speed = 0;
    end
    
    function close_fcn(~,~)
        loop_on = false;
        delete(fig);
    end
end