function DGK_interface = DGInterfaceStiffnessMatrix(obj)
% construct and return the interface DG stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input

% TODO: vectorize dF/dx using obj.T
DGK_interface = sparse(size(obj.x,1),size(obj.x,1));
% traverse the interface
for e = 1:size(obj.HalfEdge,1)
    if (obj.HalfEdge(e,3) ~= 0) % if it is an interior edge
        elem_minus = obj.MapHalfEdgeToElement(e);
        elem_plus = obj.MapHalfEdgeToElement(obj.MapHalfEdge(e));
        
        F_minus = obj.F(2*(elem_minus-1)+1:2*elem_minus,:);
        F_plus = obj.F(2*(elem_plus-1)+1:2*elem_plus,:);
        
        %         b_minus = obj.DGb(2*(elem_minus-1)+1:2*elem_minus);
        %         b_plus = obj.DGb(2*(elem_plus-1)+1:2*elem_plus);
        
        D_minus = obj.DmINV(2*(elem_minus-1)+1:2*elem_minus,:);
        D_plus = obj.DmINV(2*(elem_plus-1)+1:2*elem_plus,:);
        
        outward_n = obj.HalfEdgeUnitOutNormal(e,:)'; % column vector for the normal
            
        Len = obj.HalfEdgeLength(e);
        P0 = obj.nodeM(obj.DGEdge(e,1),:)';
        P1 = obj.nodeM(obj.DGEdge(e,2),:)';
        
        verts_minus = obj.elem(elem_minus,:);
        verts_plus = obj.elem(elem_plus,:);
        
        eta_t = obj.eta * (Len) * (1/obj.W(elem_minus) + 1/obj.W(elem_plus));
        
        bary = [1 0 0];
        for i = 1:3
            for j = 1:2
                for k = 1:3
                    for l = 1:2
                        A = (obj.Iv(:,j)*obj.G(i,:)*D_minus)' * (obj.Iv(:,l)*obj.G(k,:)*D_minus) + (obj.Iv(:,l)*obj.G(k,:)*D_minus)' * (obj.Iv(:,j)*obj.G(i,:)*D_minus);
                        B = 2 * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (obj.Iv(:,j)*obj.G(i,:)*D_minus)...
                            + 2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (obj.Iv(:,l)*obj.G(k,:)*D_minus);

                        c = ((obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * obj.nodeM(verts_minus(1),:)'))...
                            + (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)');
                        f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                        DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_minus(i)-1)+j) = DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_minus(i)-1)+j) + 1/2 * eta_t * f;
                        
                        A = -(obj.Iv(:,j)*obj.G(i,:)*D_minus)' * (obj.Iv(:,l)*obj.G(k,:)*D_plus) - (obj.Iv(:,l)*obj.G(k,:)*D_plus)' * (obj.Iv(:,j)*obj.G(i,:)*D_minus);
                        B = -2 * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,j)*obj.G(i,:)*D_minus)...
                            + 2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (-obj.Iv(:,l)*obj.G(k,:)*D_plus);
                        c = - (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * obj.nodeM(verts_plus(1),:)')...
                            - (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * obj.nodeM(verts_minus(1),:)');
                        f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                        DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_minus(i)-1)+j) = DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_minus(i)-1)+j) + 1/2 * eta_t * f;
                    
                        A = -(obj.Iv(:,j)*obj.G(i,:)*D_plus)' * (obj.Iv(:,l)*obj.G(k,:)*D_minus) - (obj.Iv(:,l)*obj.G(k,:)*D_minus)' * (obj.Iv(:,j)*obj.G(i,:)*D_plus);
                        B = 2 * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (-obj.Iv(:,j)*obj.G(i,:)*D_plus)...
                            - 2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,l)*obj.G(k,:)*D_minus);
                        c = - (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * obj.nodeM(verts_minus(1),:)')...
                            -(obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * obj.nodeM(verts_minus(1),:)')' * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)');
                            
                        f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                        DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_plus(i)-1)+j) = DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_plus(i)-1)+j) + 1/2 * eta_t * f;
                        
                        A = (obj.Iv(:,j)*obj.G(i,:)*D_plus)' * (obj.Iv(:,l)*obj.G(k,:)*D_plus) + (obj.Iv(:,l)*obj.G(k,:)*D_plus)' * (obj.Iv(:,j)*obj.G(i,:)*D_plus);
                        B = 2 * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,j)*obj.G(i,:)*D_plus)...
                            + 2 * (obj.Iv(:,j) *bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,l)*obj.G(k,:)*D_plus);
                        c = (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * obj.nodeM(verts_plus(1),:)')...
                            + (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * obj.nodeM(verts_plus(1),:)')' * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * obj.nodeM(verts_plus(1),:)');

                        f = obj.quadratic_line_int(P0,P1,A) + obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                        
                        DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_plus(i)-1)+j) = DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_plus(i)-1)+j) + 1/2 * eta_t * f;
                        
                        
                        if obj.DGIP
                            % the element corresponding to the inside edge
                            B = 1/2 * kron(outward_n, reshape(obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j),2,2))' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(k-1)+l)...
                                + 1/2 * kron(outward_n, reshape(obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(k-1)+l),2,2))' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j);
                            c = 1/2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * (obj.nodeM(verts_minus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(k-1)+l)...
                                + 1/2 *(obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * (obj.nodeM(verts_minus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j);
                            f = obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                            DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_minus(i)-1)+j) = DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_minus(i)-1)+j) - 1/2 * f;
                        
                            B = 1/2 * kron(outward_n, reshape(obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j),2,2))' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(k-1)+l)...
                                - 1/2 * kron(outward_n, reshape(obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(k-1)+l),2,2))' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j);
                            c = 1/2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_minus * (obj.nodeM(verts_minus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(k-1)+l)...
                                - 1/2 *(obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * (obj.nodeM(verts_plus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(i-1)+j);
                            f = obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                            DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_minus(i)-1)+j) = DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_minus(i)-1)+j) - 1/2 * f;
                    
                            % the element corresponding to the outside edge
                            B = -1/2 * kron(outward_n, reshape(obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j),2,2))' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(k-1)+l)...
                                + 1/2 * kron(outward_n, reshape(obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(k-1)+l),2,2))' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j);
                            c = -1/2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * (obj.nodeM(verts_plus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_minus) * obj.T(4*(elem_minus-1)+1:4*elem_minus,2*(k-1)+l)...
                                + 1/2 *(obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_minus * (obj.nodeM(verts_minus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j);
                            f = obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                            DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_plus(i)-1)+j) = DGK_interface(2*(verts_minus(k)-1)+l,2*(verts_plus(i)-1)+j) - 1/2 * f;
                        
                            
                            B = -1/2 * kron(outward_n, reshape(obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j),2,2))' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(k-1)+l)...
                                - 1/2 * kron(outward_n, reshape(obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(k-1)+l),2,2))' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j);
                            c = -1/2 * (obj.Iv(:,j) * bary(i) - obj.Iv(:,j)*obj.G(i,:)*D_plus * (obj.nodeM(verts_plus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(k-1)+l)...
                                - 1/2 * (obj.Iv(:,l) * bary(k) - obj.Iv(:,l)*obj.G(k,:)*D_plus * (obj.nodeM(verts_plus(1),:)'))' * kron(outward_n,obj.Iv)' * obj.FourthOrderTensor(elem_plus) * obj.T(4*(elem_plus-1)+1:4*elem_plus,2*(i-1)+j);
                            f = obj.linear_line_int(P0,P1,B) + obj.const_line_int(P0,P1,c);
                            DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_plus(i)-1)+j) = DGK_interface(2*(verts_plus(k)-1)+l,2*(verts_plus(i)-1)+j) - 1/2 * f;
                        
                        end
                    end
                end
            end
        end
    end
end

assert(max(max(DGK_interface - DGK_interface')) < 1e-6); % stiffness matrix should be symmetric

end