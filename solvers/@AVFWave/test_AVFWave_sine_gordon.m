function test_AVFWave_sine_gordon

for test_alphas = 1:5
    for test_ks = 4
        close all
        
        
        d = 2*pi; % domain from -d to +d
        Alpha = 10^(test_alphas);
        n = 50; % number of grid points
        h = 2 * d / (n); % size of the x interval
        I = speye(n);
        % u = zeros( 2 * n, 1);
        A = sparse( 2 * n, 2 * n );
        A( 1 : n, (n + 1): 2 * n) = speye( n );
        DD = (-2 * diag( ones(n, 1)) + diag( ones(n-1,1),1) + diag(ones(n-1,1),-1));
        DD(1,end) = 1; DD(end,1) = 1; % discretized second derivative
        sDD = sparse(DD);
        A( (n + 1):2*n, 1:n) = sDD/h/h;
        B = sparse(2*n, 2*n);
        B(1:n, (n+1):2*n) = I;
        B( (n+1):2*n,1:n) = -I;
        C = sparse(2*n,2*n); C(1:n,1:n) = I;
        
        
        t0 = 0; % initial time
        k = 10^(-test_ks + 1); % time step
        T = t0;
        t = t0;
        tend = 5;
        N = tend/k; % number of steps
        
        x0 = zeros(1,n);
        v0 = 8./ cosh(2*((-d+h) : h : d));
        % x0 = sin((-d+h) : h : d);
        % v0 = zeros(1,n);
        u0 = [x0, v0]; % initial condition
        
        %% ODE function
        % u' = f(t, u)
        %
        % Here f is autonomous, so f(t, u) := f(u)
        f = @(t, u) (A * u + B * Alpha * sin(C*u)); % u, the ode state, is a row vector
        Jv = @(t, u) A - Alpha*spdiags(cos(u(1:n)),-n,sparse(2*n,2*n));
        
        
        

        %% AVF with quadrature approximation
        avf_quad = AVFWave(f, u0, 1, 0);
        % NonlinearSovle indicates which method to use for the nonlinear solve
        % 0: use matlab build-in fsolve
        % 1: use fixed-point iteration
        % Integration indicates the integration method to use for the
        % discrete gradient
        % 0: use 3 point Gauss Legendre quadrature
        % 1: use matlab build-in function integral
        % 2: use user provided analytical integrated function

        %  function obj = AVFWave(Hf, IC, NonlinearSolve, Integration, varargin)
        % If the integral of the RHS is not provided explicitly, use
        % Integration = 0 or 1 for numerical
        % integration. If the analytical integrated function is
        % provided, Integration = 2
        % varargin list if using numerical integral:
        % {}
        % varargin list if not using numerical integral:
        % {H_intf}
        
        % initial energy
        H = energy(u0,h,Alpha);
        avf_quad_H = H;
        
        %% Simulation and Draw loop
        hf1 = figure;
        hf1.Position = [100 150 570 510];
        hf2 = figure;
        hf2.Position = [700 150 570 510];
        
        Time = 0:k:N*k;
        for i = 1:N
            
            
            avf_quad.step(k);
            avf_quad_state = avf_quad.GetCurrentState;
            avf_quad_position = avf_quad_state(1:n);
            figure(hf1)
            
            plot((-d+h) : h : d, avf_quad_position)
            title(['avf_quad' ' SG' ' alpha ' num2str(Alpha)])
            figure(hf2)
            
            avf_quad_H(i+1) = energy(avf_quad_state,h,Alpha);
            scatter(Time(1:i+1),avf_quad_H(1:i+1),'.');
            title(['avf_quad' ' SG' ' alpha ' num2str(Alpha)])
            xlim([t0,N*k])
            
            %     drawnow
        end
        
        avf_quad.SaveStates;
    end
end

end


function H = energy(u,h,Alpha)
% A sub-routine to calculate an energy of the system. Notice there may be a
% few type of energies, but here we only use one that is the first integral
% of the Hamiltonian system
n = length(u)/2;
L = (-2 * diag( ones(n, 1)) + diag( ones(n-1,1),1) + diag(ones(n-1,1),-1));
L(1,end) = 1; L(end,1) = 1;
% h is the spacing in x
x = u(1:end/2);
v = u(end/2 + 1:end);

H = (1/2 * (v * v') + 1/2/h/h * x * (-L) * x' + Alpha * (sum(1 -cos(x))))*h;
end