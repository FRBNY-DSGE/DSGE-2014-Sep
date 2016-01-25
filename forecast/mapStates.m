function [A, states_names] = mapStates(mspec,nstate,nplotstate,all_states)

% Modified by MDC on 2/18/2014 to have extra output
%   "states_names" to facilitate plotting

% Modified by MDC on 3/31/2014 to load all states if nplotstate==nstate
%   which will be used to run forecast_mode_est_ant to forecast/save all states

A = zeros(nplotstate,nstate);
eval(['states',num2str(mspec)]);

if all_states
    A(1:nplotstate, 1:nplotstate) = eye(nplotstate); 
else
    switch mspec
        case {45,803,4501, 805, 825, 826, 827, 829, 830, 835, 8351, 8352, 836, 8361, 8362}
            %% Investment
            A(1,i_t) = 1;

            %% Consumption
            A(2,c_t) = 1;

            %% Output
            A(3,y_t) = 1;

            %% Net worth 
            A(4,:) = NaN;

            %% Leverage
            A(5,:) = NaN;
            
            %% Utilization
            A(6,u_t) = 1;
            
            %% qk_t
            A(7,:) = NaN;
            
            %% Capital
            A(8,k_t) = 1;
            
        case 50
            %% Investment
            A(1,i_t) = 1;

            %% Consumption
            A(2,c_t) = 1;

            %% Output
            A(3,y_t) = 1;
            
            %% Net worth 
            A(4,n_t) = 1;

            %% Leverage
            A(5,k_t) = 1;
            A(5,n_t) = -1;
            
            %% Utilization
            A(6,u_t) = 1;
            
            %% q_t
            A(7,qk_t) = 1;
            
            %% Capital
            A(8,k_t) = 1;

        case {51,510, 555, 556, 557, 5571, 558, 804, 514, 5143 5144, 906, 5131, 5155}
            states_names = cell(nplotstate, 1);

            %% Investment
            states_names{1} = 'Investment';
            A(1,i_t) = 1;

            %% Consumption
            states_names{2} = 'Consumption';
            A(2,c_t) = 1;

            %% Output
            states_names{3} = 'Output';
            A(3,y_t) = 1;

            %% Net worth 
            states_names{4} = 'Net worth';
            A(4,n_t) = 1;

            %% Leverage
            states_names{5} = 'Leverage';
            A(5,k_t) = 1;
            A(5,n_t) = -1;
            
            %% Utilization
            states_names{6} = 'Utilization';
            A(6,u_t) = 1;
            
            %% qk_t
            states_names{7} = 'qk_t';
            A(7,qk_t) = 1;
            
            %% Productivity
            states_names{8} = 'Productivity';
            %A(8,k_t) = 1;
            A(8,z_t) = 1;
            
            %% installed capital stock or level process
            if nplotstate == 9
                if any(mspec == [557 5571 558])
                    states_names{9} = 'Z-Lev';
                    A(9, zlev_t) = 1;   
                else
                    A(9,kbar_t) = 1;
                end
            end
            %% plus rk(t)
            if nplotstate == 10
                A(9,kbar_t) = 1;
                A(10,E_Rktil+10)  = 1;
            end

            if nplotstate == 12
                A(9,R_t)  = 1;
                A(10,pi_t) = 1;
                A(11,pi_t1)  = 1;
                A(12,pi_t2)  = 1;
            end
            
            if mspec == 5571
                A(10, y_t_n) = 1;
                A(11, r_t_n) = 1;
                A(12, zlev_t) = 1;
            end

        case {828, 8286, 8287}
            %% Output Gap
            A(1,y_f_t)=-1;
            A(1,y_t)=1;
        case {8281, 8282, 8283, 8284, 8285}
            %% Don't do anything

        case {904, 908, 909, 9043, 9044, 9045, 9046, 90451, 90452, 907, 910}

            %% Investment
            states_names{1} = 'Investment';
            A(1,i_t) = 1;
            %A(1,gamm_t) = 1;

            %% Consumption
            states_names{2} = 'Consumption';
            A(2,c_t) = 1;

            %% Output
            states_names{3} = 'Output';
            A(3,y_t) = 1;

            %% Net worth 
            states_names{4} = 'Net Worth';
            A(4,n_t) = 1;

            %% Leverage
            states_names{5} = 'Leverage';
            A(5,k_t) = 1;
            A(5,n_t) = -1;
            
            states_names{6} = 'Output Gap';
            A(6,y_t) = 1;
            A(6,y_f_t) = -1;
            
            %% qk_t
            states_names{7} = 'Marginal Cost';
            A(7,mc_t) = 1;
            
            %% Productivity
            states_names{8} = 'Productivity';
            %A(8,k_t) = 1;
            A(8,z_t) = 1;
            
            %% plus rk(t)
            states_names{9} = 'y_f_t';
            states_names{10} = 'rk(t)';
            A(9,y_f_t) = 1;
            A(10,E_Rktil+10)  = 1;
            
            %% installed capital stock or level process
            if nplotstate == 9
                states_names{9} = 'Capital Stock';
                A(9,mc_t) = 1;
            end


           % taylor rule variables
            if nplotstate == 11
                states_names{10} = 'y_f_t';
                A(11,y_f_t) = 1;
                A(13, pitil_t) = 1;
                A(14, Y_t) = 1;
                A(15, Ybar_t) = 1;
                A(16, pitil_t) = 1;
            end
            
            % %% natural rates: output and ffr
            % if nplotstate == 13
                % states_names{10} = 'Natural Output';
                % states_names{11} = 'Natural Rate';
                % states_names{12} = 'Nstate';
                % states_names{13} = 'Natural FFR';
                % A(10,y_f_t) = 1;
                % A(11,r_f_t) = 1;
                % A(12,nstate) = 1;
                % A(13,i_f_t) = 1;
            % end

        case {920, 921}

            states_names{1} = 'Output';
            A(1,y_t) = 1;

            states_names{2} = 'Flex Output';
            A(2,y_f_t) = 1;

            states_names{3} = 'Interest Rate';
            A(3,R_t) = 1;

            states_names{4} = 'Natural Rate';
            A(4,r_f_t) = 1;

            states_names{5} = 'Inflation';
            A(5,pi_t) = 1;

            states_names{6} = 'Annual Inflation';
            A(6,pi_a_t) = 1;
            
            states_names{7} = 'Output Gap';
            A(7,y_t) = 1;
            A(7,y_f_t) = -1;
            
            states_names{8} = 'Expected Inflation';
            A(8,E_pi) = 1;

            states_names{9} = 'Ex Ante Real Rate ';
            A(9,R_t) = 1;
            A(9,E_pi) = -1;

            states_names{10} = 'Long Run Inflation';
            A(10,pist_t) = 1;
           
            states_names{11} = 'Inflation (t-1)';
            A(11,pi_t1) = 1;
            
            states_names{12} = 'Inflation (t-2)';
            A(12,pi_t2) = 1;
            
            states_names{13} = 'Interest Rate (t-1)';
            A(13,R_t1) = 1;

            states_names{14} = 'Monetary Policy Shock';
            A(14,rm_t) = 1;

            %states_names{10} = 'Productivity';
            %A(10,z_t) = 1;

        otherwise, 
            error('"A" matrix not specified for this mspec')
    end
end
if ~exist('states_names')
    states_names = cell(arrayfun(@(x) ['state' num2str(x)], [1:nplotstate], 'UniformOutput',0)');
end
