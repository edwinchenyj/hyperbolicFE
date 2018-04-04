function K = StiffnessMatrix(obj)
% construct and return the stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input
assert(obj.finalized)

index = 1;
for s = 1:obj.NS
    
    ks = obj.ks(s);
    l0 = obj.W(s);
    F = obj.F(3*(s-1)+1:3*s,:);
    
    Kii = 1/2 * ks* (2 * F *(F') + ((F') * F - 1)*eye(3));
%     Ks = [Kii -Kii; -Kii Kii]; 
%     obj.C = C;
%     disp('dFdx = [')
%     disp(tT)
%     disp('];')
%     obj.dFdx = tT;
%     disp('C= [')
%     disp(C)
%     disp('];')
    % element stiffness matrix
    
    for si = 1:2
        for sj = 1:2
            
            if si == sj
                sA(index:index+8) = Kii(:);
            else
                sA(index:index+8) = -Kii(:);
            end
            % the same as: 
            % sA(index:index+8) = reshape(Ks(3*(si-1)+1:3*si,3*(sj-1)+1:3*sj),9,1);
            index = index + 9;
        end
    end
end

% global stiffness matrix
K = sparse(obj.ii, obj.jj, sA, 3*obj.N, 3*obj.N);
end