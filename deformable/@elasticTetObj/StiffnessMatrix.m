function K = StiffnessMatrix(obj)
% construct and return the stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input
index = 1;
sA = zeros(size(obj.ii));
if length(obj.material_type) == 1
    switch obj.material_type % TODO: change the material type initialization
        case 1
            for t = 1:obj.NT
                
                mu = obj.mu;
                lambda = obj.lambda;
                tT = obj.T(9*(t-1)+1:9*t,:);
                W = obj.W(t);
                tF = obj.F(3*(t-1)+1:3*t,:);
                
                % the fourth order tensor
                tFINV = obj.FINV(3*(t-1)+1:3*t,:);
%                 Kk = obj.Kmm* kron(tFINV',tFINV);
                % for Neo-hookean
%                 C = mu * obj.Im + mu * Kk...
%                     - lambda * (log(det(tF))*Kk)...
%                     + lambda*(obj.Kmm*(tFINV(:)*reshape(transpose(tFINV),1,9)));
                C = neohookean_c_mex(mu,lambda,tF,tFINV,obj.Im,obj.Kmm);
                % element stiffness matrix
                
                Kt = W * tT'*C*tT;
                Kt = 1/2 * (Kt + Kt');
                %         if det(Kt) < 0
                [U,D] = eig(Kt);
                for i_s = 1:size(D,1)
                    if D(i_s,i_s) < 1e-6
                        %                     disp(D(i_s,i_s))
                        %                     disp('element eigenvalue < 0')
                        D(i_s,i_s) = 1e-3;
                        
                    end
                end
                Kt = U * D *(U');
                %             eig(Kt)
                %             assert(~any(eig(Kt) < 0));
                %         end
                sA(index:index+143) = Kt(:);
                index = index + 144;

                
            end
            % global stiffness matrix
            K = sparse(obj.ii, obj.jj, sA, obj.Dim*obj.N, obj.Dim*obj.N);
        case 2
            if isempty(obj.K0)
                for t = 1:obj.NT
                    
                    mu = obj.mu;
                    lambda = obj.lambda;
                    tT = obj.T(9*(t-1)+1:9*t,:);
                    W = obj.W(t);
                    %                 tF = obj.F(3*(t-1)+1:3*t,:);
                    % for linear elasticity
                    C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
                    
                    Kt = W * tT'*C*tT;
                    Kt = 1/2 * (Kt + Kt');
                    sA(index:index+143) = Kt(:);
                    index = index + 144;
                    
%                     for ti = 1:4
%                         for tj = 1:4
%                             % same as, but faster:
%                             %             sA(index:index+8) = reshape(Kt(3*(ti-1)+1:3*ti,3*(tj-1)+1:3*tj),9,1);
%                             sA(index:index+8) = Kt(obj.IndK(4*(ti-1) + tj,:));
%                             
%                             index = index + 9;
%                         end
%                     end
                end
                % global stiffness matrix
                obj.K0 = sparse(obj.ii, obj.jj, sA, obj.Dim*obj.N, obj.Dim*obj.N);
                
            end
            K = obj.K0;
        case 3
            for t = 1:obj.NT
                
                mu = obj.mu;
                lambda = obj.lambda;
                tT = obj.T(9*(t-1)+1:9*t,:);
                W = obj.W(t);
                tF = obj.F(3*(t-1)+1:3*t,:);
%                 FTF = (tF')*tF;
                % StVK
                C = stvk_c_mex(mu,lambda,tF,obj.Iv,obj.Im,obj.Kmm);
                
                % element stiffness matrix
                
                Kt = W * tT'*C*tT;
                Kt = 1/2 * (Kt + Kt');
                %         if det(Kt) < 0
                [U,D] = eig(Kt);
%                 [U,D] = kt_eig_mex(Kt);
                
                for i_s = 1:size(D,1)
                    if D(i_s,i_s) < 1e-6
                        D(i_s,i_s) = 1e-3;
                        
                    end
                end
                Kt = U * D * (U');
                
                sA(index:index+143) = Kt(:);
                index = index + 144;

         
            end
            K = sparse(obj.ii, obj.jj, sA, obj.Dim*obj.N, obj.Dim*obj.N);
    end
else
    for t = 1:obj.NT
        switch obj.material_type(t)
            case 1
                
                
                mu = obj.mu(t);
                lambda = obj.lambda(t);
                tT = obj.T(9*(t-1)+1:9*t,:);
                W = obj.W(t);
                tF = obj.F(3*(t-1)+1:3*t,:);
                
                % the fourth order tensor
                tFINV = obj.FINV(3*(t-1)+1:3*t,:);
%                 Kk = obj.Kmm* kron(tFINV',tFINV);
%                 % for Neo-hookean
%                 C = mu * obj.Im + mu * Kk...
%                     - lambda * (log(det(tF))*Kk)...
%                     + lambda*(obj.Kmm*(tFINV(:)*reshape(transpose(tFINV),1,9)));
                C = neohookean_c_mex(mu,lambda,tF,tFINV,obj.Im,obj.Kmm);
                % element stiffness matrix
                
                Kt = W * tT'*C*tT;
                Kt = 1/2 * (Kt + Kt');
                %         if det(Kt) < 0
                [U,D] = eig(Kt);
                for i_s = 1:size(D,1)
                    if D(i_s,i_s) < 1e-6
                        %                     disp(D(i_s,i_s))
                        %                     disp('element eigenvalue < 0')
                        D(i_s,i_s) = 1e-3;
                        
                    end
                end
                Kt = U * D *(U');
                
                sA(index:index+143) = Kt(:);
                index = index + 144;

                
            case 2
                
                mu = obj.mu(t);
                lambda = obj.lambda(t);
                tT = obj.T(9*(t-1)+1:9*t,:);
                W = obj.W(t);
                %                 tF = obj.F(3*(t-1)+1:3*t,:);
                % for linear elasticity
                C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
                
                Kt = W * tT'*C*tT;
                Kt = 1/2 * (Kt + Kt');
                
                sA(index:index+143) = Kt(:);
                index = index + 144;

                
            case 3
                
                
                mu = obj.mu(t);
                lambda = obj.lambda(t);
                tT = obj.T(9*(t-1)+1:9*t,:);
                W = obj.W(t);
                tF = obj.F(3*(t-1)+1:3*t,:);
                % TODO: add a switch for mex
                % StVK 
                C = mu * (kron((tF')*tF,obj.Iv) + kron(tF',tF) * obj.Kmm...
                    + (kron(obj.Iv,tF*(tF'))- obj.Im)) ...
                    + lambda * (trace(1/2 * ((tF')*tF-obj.Iv)) ...
                    + reshape(tF,[],1) * reshape(tF,[],1)');
%                 C = stvk_c_mex(mu,lambda,tF,obj.Iv,obj.Im,obj.Kmm);
                
                % element stiffness matrix
                
                Kt = W * tT'*C*tT;
                Kt = 1/2 * (Kt + Kt');
                %         if det(Kt) < 0
                [U,D] = eig(Kt);
%                 [U,D] = kt_eig_mex(Kt);
                
                for i_s = 1:size(D,1)
                    if D(i_s,i_s) < 1e-6
                        D(i_s,i_s) = 1e-3;
                        
                    end
                end
                Kt = U * D * (U');
%                 Kt = ek_projection_mex(Kt);
%                 Kt = ek_projection(Kt);
                
                sA(index:index+143) = Kt(:);
                index = index + 144;

        end
    end
    % global stiffness matrix
    K = sparse(obj.ii, obj.jj, sA, obj.Dim*obj.N, obj.Dim*obj.N);
end
end