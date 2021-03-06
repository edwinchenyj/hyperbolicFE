draw = true;
rerun_flag = true;
save_state = true;
dt_list = [1/100];
T_list = [2];
mesh_shape = 'circle';
simulation_type_list = {'CG'};
solver_list = {'BE'};
maxA_list = [0.01 0.001];
Y_list = [100000, 1000000];
P_list = [0.48 0.4];
rho_list = [1000];
a_list = [0];
b_list = [0];
constraint_list = {'none'};
material_list = {'neo-hookean'};
gravity = 'off';

for i_maxA = 1:length(maxA_list)
    maxA = maxA_list(i_maxA);
    for i_Y = 1:length(Y_list)
        Y = Y_list(i_Y);
        for i_P = 1:length(P_list)
            P = P_list(i_P);
            for i_rho = 1:length(rho_list)
                rho = rho_list(i_rho);
                for i_a = 1:length(a_list)
                    a = a_list(i_a);
                    for i_b = 1:length(b_list)
                        b = b_list(i_b);
                        for i_dt = 1:length(dt_list)
                            dt = dt_list(i_dt);
                            for i_T = 1:length(T_list)
                                T = T_list(i_T);
                                for i_simulation_type = 1:length(simulation_type_list)
                                    simulation_type = simulation_type_list{i_simulation_type};
                                    for i_solver = 1:length(solver_list)
                                        solver = solver_list{i_solver};
                                        for i_material = 1:length(material_list)
                                            material = material_list{i_material};
                                            for i_constraint_list = 1:length(constraint_list)
                                            constraint = constraint_list{i_constraint_list};
%                                              
%                                             expression = sprintf('simulation(''draw'',%i,''rerun_flag'',%i,''save_state'',%i,''dt'',%f,''T'',%f,''mesh_shape'',''%s'',''maxA'',%f,''simulation_type'',''%s'',''solver'',''%s'',''Y'',%f,''P'',%f,''rho'',%f,''a'',%f,''b'',%f,''material'',''%s'')',...
%                                                 draw,rerun_flag,save_state,dt,T,mesh_shape,maxA,simulation_type,solver,Y,P,rho,a,b,material);
%                                             eval(expression);
simulation('draw',draw,...
    'rerun_flag',rerun_flag,...
    'save_state',save_state,...
    'dt',dt,...
    'T',T,...
    'mesh_shape',mesh_shape,...
    'maxA',maxA,...
    'simulation_type',simulation_type,...
    'solver',solver,...
    'Y',Y,...
    'P',P,...
    'rho',rho,...
    'a',a,...
    'b',b,...
    'material',material,...
    'constraint',constraint,...
    'gravity',gravity)
%                                                 draw,rerun_flag,save_state,dt,T,mesh_shape,maxA,simulation_type,solver,Y,P,rho,a,b,material)
                                        end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end