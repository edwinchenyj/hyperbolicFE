%####################################################################
% This function tests calc_meshdata by integrating a function over
% the unit square in 2D. The domain is divided into triangles, and
% the integration quadrature is defined on the unit triangle, so the
% integration over the domain is performed by chage of variables.
% This function also draws the integration poins for demonstration.
%####################################################################

function test_2d_drawpoints

    cd 2d
    [p,e,t] = initmesh('geom','Hmax',0.2);
    cd ..
    
    % meshdata
    md = calc_meshdata(2,p,e,t);

    % integration quadrature in 2d
    [ip,w] = inttri(2);
    nip = size(ip,2);

    % plot the quad points in the unit triangle
    subplot(1,2,1)
    plot([0 1 0 0],[0 0 1 0],'--'); hold on;
    plot(ip(1,:),ip(2,:),'*'); axis equal tight;

    % plot the mesh
    subplot(1,2,2)
    pdemesh(md.p,md.e,md.t); hold on; axis equal tight;

    % calculate the integral of the function 'fun'
    integral = 0;
    for i=1:size(md.t,2)
        % get element info
        B_K      = md.B_K(:,:,i);
        b_K      = md.b_K(:,i);
        B_K_detA = abs(md.B_K_det(i));

        % calculate the affine mapping F_K( ip )
        F_K_ip = B_K * ip + b_K(:,ones(1,nip));
        plot(F_K_ip(1,:),F_K_ip(2,:),'*');
        
        % the values of the function to be integrated
        funval = fun(F_K_ip);
        
        % calculate the integral contribution
        integral = integral + sum( w' .* B_K_detA .* funval );
    end

    disp(['The integral value:  ',num2str(integral)])
    disp(['Domain area:         ',num2str(sum(md.areas))])

end


%####################################################################
% this is the function we are integrating
%####################################################################

function val = fun(p)
    %val = ones(1,size(p,2));    %fun = 1
    val = p(1,:);               %fun = x
    %val = p(2,:).^2;            %fun = y^2
end