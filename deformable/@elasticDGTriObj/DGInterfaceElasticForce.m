function DGf_k = DGInterfaceElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

DGf_k = zeros(2*obj.DGN,1);


% traverse the interface
for e = 1:size(obj.HalfEdge,1)
    if (obj.HalfEdge(e,3) ~= 0) % if it is an interior edge
        elem_minus = obj.MapHalfEdgeToElement(e);
        elem_plus = obj.MapHalfEdgeToElement(obj.MapHalfEdge(e));
        
        F_minus = obj.F(2*(elem_minus-1)+1:2*elem_minus,:);
        F_plus = obj.F(2*(elem_plus-1)+1:2*elem_plus,:);
        
        b_minus = obj.DGb(2*(elem_minus-1)+1:2*elem_minus);
        b_plus = obj.DGb(2*(elem_plus-1)+1:2*elem_plus);
        
        D_minus = obj.DmINV(2*(elem_minus-1)+1:2*elem_minus,:);
        D_plus = obj.DmINV(2*(elem_plus-1)+1:2*elem_plus,:);
        
        Len = obj.HalfEdgeLength(e);
        P0 = obj.DGnodeM(obj.DGEdge(e,1),:)';
        P1 = obj.DGnodeM(obj.DGEdge(e,2),:)';
        
        verts_minus = obj.DGelem(elem_minus,:);
        verts_plus = obj.DGelem(elem_plus,:);
        
        eta_t = obj.eta * Len * (1/obj.W(elem_minus) + 1/obj.W(elem_plus));
        
        bary = [1 0 0]; % used to calculate derivatives of b
        
        for i = 1:3
            for j = 1:2
                % the element corresponding to the inside edge
                A = (obj.I2(:,j)*obj.G(i,:)*D_minus)' * F_minus + F_minus' * (obj.I2(:,j)*obj.G(i,:)*D_minus)...
                    - (obj.I2(:,j)*obj.G(i,:)*D_minus)' * F_plus - F_plus' * (obj.I2(:,j)*obj.G(i,:)*D_minus);
                B =  2*(b_minus - b_plus)' * (obj.I2(:,j)*obj.G(i,:)*D_minus) ...
                    + 2 *(obj.I2(:,j) * bary(i) - obj.I2(:,j)*obj.G(i,:)*D_minus * obj.DGnodeM(verts_minus(1),:)')' * (F_minus - F_plus);
                c = (obj.I2(:,j) * bary(i) - obj.I2(:,j)*obj.G(i,:)*D_minus * (obj.DGnodeM(verts_minus(1),:)'))' * (b_minus - b_plus)...
                    + ((b_minus - b_plus)') * (obj.I2(:,j) * bary(i) - obj.I2(:,j)*obj.G(i,:)*D_minus * (obj.DGnodeM(verts_minus(1),:)'));
                f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);    
                DGf_k(2*(verts_minus(i)-1)+j) = DGf_k(2*(verts_minus(i)-1)+j) - 1/2 *  eta_t * f;
                
                
                
                % the element corresponding to the outside edge
                A = (obj.I2(:,j)*obj.G(i,:)*D_plus)' * F_plus + F_plus' * (obj.I2(:,j)*obj.G(i,:)*D_plus)...
                    - (obj.I2(:,j)*obj.G(i,:)*D_plus)' * F_minus - F_minus' * (obj.I2(:,j)*obj.G(i,:)*D_plus);
                B = 2 * (b_minus - b_plus)' * (-obj.I2(:,j)*obj.G(i,:)*D_plus)...
                    - 2 * (obj.I2(:,j) * bary(i) - obj.I2(:,j)*obj.G(i,:)*D_plus * obj.DGnodeM(verts_plus(1),:)')' * (F_minus - F_plus);
                c = - (obj.I2(:,j) * bary(i) - obj.I2(:,j)*obj.G(i,:)*D_plus * (obj.DGnodeM(verts_plus(1),:)'))' * (b_minus - b_plus)...
                    - (b_minus - b_plus)' * (obj.I2(:,j) * bary(i) - obj.I2(:,j)*obj.G(i,:)*D_plus * (obj.DGnodeM(verts_plus(1),:)'));
                f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                DGf_k(2*(verts_plus(i)-1)+j) = DGf_k(2*(verts_plus(i)-1)+j) - 1/2 * eta_t * f;
%                 
            end
        end
    end
end

end