draw = true;
rerun_flag = true;
save_state = true;
dt_list = [1/10000];
T_list = [2];
mesh_shape = 'small_bar';
simulation_type_list = {'CG'};
solver_list = {'BE'};
h_list = [0.5];
Y_list = [1e4];
P_list = [0.45];
rho_list = [1000];
a_list = [0];
b_list = [0.000];
material_list = {'neo-hookean'};
gravity = 'on';
axis_box = [-1 1 -1 1 -1 1]/4;

for i_h = 1:length(h_list)
    h = h_list(i_h);
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
                                            material_type = material_list{i_material};
%                                              
%                                             expression = sprintf('simulation(''draw'',%i,''rerun_flag'',%i,''save_state'',%i,''dt'',%f,''T'',%f,''mesh_shape'',''%s'',''h'',%f,''simulation_type'',''%s'',''solver'',''%s'',''Y'',%f,''P'',%f,''rho'',%f,''a'',%f,''b'',%f,''material'',''%s'')',...
%                                                 draw,rerun_flag,save_state,dt,T,mesh_shape,h,simulation_type,solver,Y,P,rho,a,b,material);
%                                             eval(expression);
simulation3D_motions_bar_polyfit_modH_line_search('draw',draw,...
    'rerun_flag',rerun_flag,...
    'save_state',save_state,...
    'dt',dt,...
    'T',T,...
    'mesh_shape',mesh_shape,...
    'h',h,...
    'simulation_type',simulation_type,...
    'solver',solver,...
    'Y',Y,...
    'P',P,...
    'rho',rho,...
    'a',a,...
    'b',b,...
    'material',material_type,...
    'gravity',gravity,...
    'axis_box',axis_box)
%                                                 draw,rerun_flag,save_state,dt,T,mesh_shape,h,simulation_type,solver,Y,P,rho,a,b,material)

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