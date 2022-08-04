function scene_demo
    
    cam_pos = [0,0,0];
    cam_el = .5;
    cam_az = .5;
    cam_fwd_speed = 0.01;

    unit_vec = @(el,az)[
        cos(el)*cos(az),...
        cos(el)*sin(az),...
        sin(el)];

    fig = figure(...
        'color', [0.2,0.2,0.2],...
        'KeyPressFcn', @key_down, ...
        'KeyReleaseFcn', @key_up, ...
        'CloseRequestFcn',@close_fcn);
    ax = axes(...
        'projection','perspective',...
        'CameraViewAngle',60,...
        'DataAspectRatio',[1,1,1],...
        'visible', 'off');
    rendererinfo(ax)
    for i = 1:100
    patch(...
        'xdata', 10*rand*[1,1,1] + rand(1,3), ...
        'YData', 10*rand*[1,1,1] + rand(1,3), ...
        'ZData', 10*rand*[1,1,1] + rand(1,3), ...)
        'facecolor', rand(1,3))
%     pause(.02);
    end

    loop_on = true;
    while loop_on
        cam_pos = cam_pos + ...
            cam_fwd_speed * unit_vec(cam_el,cam_az);
        ax.CameraPosition = cam_pos;
        ax.CameraTarget = cam_pos + unit_vec(cam_el,cam_az);
        drawnow;
    end

    function close_fcn(~,~)
        loop_on = false;
        delete(fig)
    end

end