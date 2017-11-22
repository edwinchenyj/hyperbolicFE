function DGf_k = DGInterfaceElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

DGf_k = zeros(2*obj.N,1);
if isempty(obj.DGKi)
if obj.DGBZ
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
            P0 = obj.nodeM(obj.DGEdge(e,1),:)';
            P1 = obj.nodeM(obj.DGEdge(e,2),:)';
            
            
            outward_n = obj.HalfEdgeUnitOutNormal(e,:)'; % column vector for the normal
            
            
            verts_minus = obj.elem(elem_minus,:);
            verts_plus = obj.elem(elem_plus,:);
            
            eta_t = obj.eta * Len * (1/obj.W(elem_minus) + 1/obj.W(elem_plus));
            
            bary = [1 0 0]; % used to calculate derivatives of b
            
            for i = 1:3
                for j = 1:2
                    % the element corresponding to the inside edge
                    A = (obj.Iv(:,j)*obj.G(i,:)*D_minus)' * F_minus + F_minus' * (obj.Iv(:,j)*obj.G(i,:)*D_minus)...
                        - (obj.Iv(:,j)*obj.G(i,:)*D_minus)' * F_plus - F_plus' * (obj.Iv(:,j)*obj.G(i,:)*D_minus);
                    B =  2*(b_minus - b_plus)' * (obj.Iv(:,j)*obj.G(i,:)*D_minus) ...
                        + 2 *(obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (F_minus - F_plus);
                    c = (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * (obj.nodeM(verts_minus(1),:)'))' * (b_minus - b_plus)...
                        + ((b_minus - b_plus)') * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * (obj.nodeM(verts_minus(1),:)'));
                    f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                    DGf_k(2*(verts_minus(i)-1)+j) = DGf_k(2*(verts_minus(i)-1)+j) - 1/2 *  eta_t * f;
                    
                    % the element corresponding to the outside edge
                    A = (obj.Iv(:,j)*obj.G(i,:)*D_plus)' * F_plus + F_plus' * (obj.Iv(:,j)*obj.G(i,:)*D_plus)...
                        - (obj.Iv(:,j)*obj.G(i,:)*D_plus)' * F_minus - F_minus' * (obj.Iv(:,j)*obj.G(i,:)*D_plus);
                    B = 2 * (b_minus - b_plus)' * (-obj.Iv(:,j)*obj.G(i,:)*D_plus)...
                        - 2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (F_minus - F_plus);
                    c = - (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * (obj.nodeM(verts_plus(1),:)'))' * (b_minus - b_plus)...
                        - (b_minus - b_plus)' * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * (obj.nodeM(verts_plus(1),:)'));
                    f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                    DGf_k(2*(verts_plus(i)-1)+j) = DGf_k(2*(verts_plus(i)-1)+j) - 1/2 * eta_t * f;
                    
                    if obj.DGIP
                        
                        % the element corresponding to the inside edge
                        B = 1/2 * kron(outward_n, reshape(obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j),2,2))' * reshape(obj.Stress(elem_minus) + obj.Stress(elem_plus),[],1)...
                            + 1/2 * (kron(outward_n,F_minus) - kron(outward_n,F_plus))' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j);
                        c = 1/2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * (obj.nodeM(verts_minus(1),:)'))' * kron(outward_n,obj.Iv)' * reshape(obj.Stress(elem_minus) + obj.Stress(elem_plus),[],1)...
                            + 1/2 *((b_minus - b_plus)') * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j);
                        f = obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                        DGf_k(2*(verts_minus(i)-1)+j) = DGf_k(2*(verts_minus(i)-1)+j) + 1/2 * f;

                        % the element corresponding to the outside edge
                        B = -1/2 * kron(outward_n, reshape(obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j),2,2))' * reshape(obj.Stress(elem_minus) + obj.Stress(elem_plus),[],1)...
                            + 1/2 * (kron(outward_n,F_minus) - kron(outward_n,F_plus))' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j);
                        c = -1/2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * (obj.nodeM(verts_plus(1),:)'))' * kron(outward_n,obj.Iv)' * reshape(obj.Stress(elem_minus) + obj.Stress(elem_plus),[],1)...
                            + 1/2 *((b_minus - b_plus)') * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j);
                        f = obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                        DGf_k(2*(verts_plus(i)-1)+j) = DGf_k(2*(verts_plus(i)-1)+j) + 1/2 * f;
                        
                    end
                end
            end
        end
    end
end
else
    if ~obj.DGIP
        DGf_k = -obj.DGKi * (obj.x - obj.X);
    end
end
end