function test_ERAVF_nonlinear_wave
for test_alphas = 2
    for test_ks = 3
        close all
        d = 2*pi; % domain from -d to +d
        Alpha = 10^(test_alphas);
        n = 100; % number of grid points
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
        v0 = Alpha * 8./ cosh(2*((-d+h) : h : d));
        % x0 = sin((-d+h) : h : d);
        % v0 = zeros(1,n);
        u0 = [x0, v0]; % initial condition
        
        %% ODE function
        % u' = f(t, u)
        %
        % Here f is autonomous, so f(t, u) := f(u)
f = @(t, u) (A * u + B * Alpha * C*(u.^3)); % u, the ode state, is a row vector
Jv = @(t, u) A - Alpha*spdiags(3*(u(1:n).^2),-n,sparse(2*n,2*n));

        %% ERAVF
        ERavf = ERAVF(f, u0, Jv, 1, 1, false);
        
        % function obj = ERAVF(Hf, IC, HJ, NonlinearSolve, Integration, usingKSM)
        % NonlinearSolve = 0; % indicating which method to use for the nonlinear solve
        % % 0: use matlab build-in fsolve
        % % 1: use fixed-point iteration
        % Integration = 0; % indicate which integration method to use for the discrete gradient
        % % 0: use 3 point Gauss Legendre quadrature
        % % 1: use matlab build-in function integral
        

        % initial energy
        H = energy(u0,h,Alpha);
        ERavf_H = H;
        
        %% Simulation and Draw loop
        hf1 = figure;
        hf1.Position = [100 150 570 510];
        hf2 = figure;
        hf2.Position = [700 150 570 510];
        
        Time = 0:k:N*k;
        for i = 1:N
            
            
            ERavf.step(k);
            ERavf_state = ERavf.GetCurrentState;
            ERavf_position = ERavf_state(1:n);
            figure(hf1)
            
            plot((-d+h) : h : d, ERavf_position)
            title(['ERAVF' ' NLW' ' alpha ' num2str(Alpha)])
            figure(hf2)
            
            ERavf_H(i+1) = energy(ERavf_state,h,Alpha);
            scatter(Time(1:i+1),ERavf_H(1:i+1),'.');
            title(['ERAVF' ' NLW' ' alpha ' num2str(Alpha)])
            xlim([t0,N*k])
            
                drawnow
        end
        
        ERavf.SaveStates;
        
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

H = (1/2 * (v * v') + 1/2/h/h * x * (-L) * x' + Alpha * (sum(1/4*x.^4)))*h;
end