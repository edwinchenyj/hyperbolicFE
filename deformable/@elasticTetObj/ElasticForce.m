function f = ElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

f = zeros(3*obj.N,1);
if all(obj.material_type == 2)
    if isempty(obj.K0)
        index = 1;
        for t = 1:obj.NT
            
            if length(obj.material_type) == 1
                mu = obj.mu;
                lambda = obj.lambda;
                
            else
                mu = obj.mu(t);
                lambda = obj.lambda(t);
            end
            tT = obj.T(9*(t-1)+1:9*t,:);
            W = obj.W(t);
            %                 tF = obj.F(3*(t-1)+1:3*t,:);
            % for linear elasticity
            C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
            
            Kt = W * tT'*C*tT;
            Kt = 1/2 * (Kt + Kt');
            sA(index:index+143) = Kt(:);
            index = index + 144;
        end
        % global stiffness matrix
        obj.K0 = sparse(obj.ii, obj.jj, sA, obj.Dim*obj.N, obj.Dim*obj.N);
        
    end
    f = -obj.K0 * (obj.x - obj.X);
    
elseif length(obj.material_type) == 1
    
    switch obj.material_type
        case 1
            for t = 1:obj.NT
                
                f_new = zeros(3*obj.N,1);
                
                %             type = obj.material_type;
                tF = obj.F(3*(t-1)+1:3*t,:);
                
                mu = obj.mu;
                lambda = obj.lambda;
                
                tFINV = obj.FINV(3*(t-1)+1:3*t,:);
                
                 J = det(tF);
                 P = mu *(tF - tFINV') + lambda * log(J) * tFINV';
%                 P = neohookean_P_mex(mu,lambda,tF,tFINV);
                
                H = -obj.W(t) * P * (obj.DmINV(3*(t-1)+1:3*t,:)');
                i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3); l = obj.elem(t, 4);
                
                f_new(3*(i-1)+1:3*i) = f_new(3*(i-1)+1:3*i)+H(:,1);
                f_new(3*(j-1)+1:3*j) = f_new(3*(j-1)+1:3*j)+H(:,2);
                f_new(3*(k-1)+1:3*k) = f_new(3*(k-1)+1:3*k)+H(:,3);
                f_new(3*(l-1)+1:3*l) = f_new(3*(l-1)+1:3*l) - H(:,1) - H(:,2) - H(:,3);
                f = f + f_new;
            end
        case 3
            for t = 1:obj.NT
                
                f_new = zeros(3*obj.N,1);
                
                %             type = obj.material_type;
                tF = obj.F(3*(t-1)+1:3*t,:);
                
                mu = obj.mu;
                lambda = obj.lambda;
                
%                 P = (tF)*(mu * ((tF') * tF - obj.Iv)) ...
%                     + (tF) * (lambda * trace((1/2) * ((tF') * tF - obj.Iv)) * obj.Iv);
                P = stvk_P_mex(mu,lambda,tF,obj.Iv);
                
                
                H = -obj.W(t) * P * (obj.DmINV(3*(t-1)+1:3*t,:)');
                i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3); l = obj.elem(t, 4);
                
                f_new(3*(i-1)+1:3*i) = f_new(3*(i-1)+1:3*i)+H(:,1);
                f_new(3*(j-1)+1:3*j) = f_new(3*(j-1)+1:3*j)+H(:,2);
                f_new(3*(k-1)+1:3*k) = f_new(3*(k-1)+1:3*k)+H(:,3);
                f_new(3*(l-1)+1:3*l) = f_new(3*(l-1)+1:3*l) - H(:,1) - H(:,2) - H(:,3);
                f = f + f_new;
            end
    end
else
    for t = 1:obj.NT
        
        switch obj.material_type(t)
            case 1
                
                f_new = zeros(3*obj.N,1);
                
                %             type = obj.material_type;
                tF = obj.F(3*(t-1)+1:3*t,:);
                
                mu = obj.mu(t);
                lambda = obj.lambda(t);
                
                tFINV = obj.FINV(3*(t-1)+1:3*t,:);
                
%                 J = det(tF);
%                 P = mu *(tF - tFINV') + lambda * log(J) * tFINV';
                P = neohookean_P_mex(mu,lambda,tF,tFINV);
                
                H = -obj.W(t) * P * (obj.DmINV(3*(t-1)+1:3*t,:)');
                i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3); l = obj.elem(t, 4);
                
                f_new(3*(i-1)+1:3*i) = f_new(3*(i-1)+1:3*i)+H(:,1);
                f_new(3*(j-1)+1:3*j) = f_new(3*(j-1)+1:3*j)+H(:,2);
                f_new(3*(k-1)+1:3*k) = f_new(3*(k-1)+1:3*k)+H(:,3);
                f_new(3*(l-1)+1:3*l) = f_new(3*(l-1)+1:3*l) - H(:,1) - H(:,2) - H(:,3);
                f = f + f_new;
                
                
                
            case 3
                
                f_new = zeros(3*obj.N,1);
                
                %             type = obj.material_type;
                tF = (obj.F(3*(t-1)+1:3*t,:));
                
                mu = obj.mu(t);
                lambda = obj.lambda(t);
                % TODO: add a switch for mex
                P = (tF)*(mu * ((tF') * tF - obj.Iv)) ...
                    + (tF) * (lambda * trace((1/2) * ((tF') * tF - obj.Iv)) * obj.Iv);
%                 isreal(tF)
%                 tF
%                 P = stvk_P_mex(mu,lambda,tF,obj.Iv);
                
                H = -obj.W(t) * P * (obj.DmINV(3*(t-1)+1:3*t,:)');
                i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3); l = obj.elem(t, 4);
                
                f_new(3*(i-1)+1:3*i) = f_new(3*(i-1)+1:3*i)+H(:,1);
                f_new(3*(j-1)+1:3*j) = f_new(3*(j-1)+1:3*j)+H(:,2);
                f_new(3*(k-1)+1:3*k) = f_new(3*(k-1)+1:3*k)+H(:,3);
                f_new(3*(l-1)+1:3*l) = f_new(3*(l-1)+1:3*l) - H(:,1) - H(:,2) - H(:,3);
                f = f + f_new;
                
        end
    end
end
end