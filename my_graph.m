%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim,figout)
            %% Plot production function.
            
            figure(1)
            
            plot(par.kgrid,sol.y)
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
            title('Production Function')
            
            fig_name = strcat(figout,'ypol.fig');
            savefig(fig_name)
            
            %% Plot capital policy function.
            
            figure(2)
            
            plot(par.kgrid,sol.k)
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$k_{t+1}$'},'Interpreter','latex') 
            title('Capital Policy Function')
            
            fig_name = strcat(figout,'kpol.fig');
            savefig(fig_name)
            
            %% Plot consumption policy function.
            
            figure(3)
            
            plot(par.kgrid,sol.c)
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function')
            
            fig_name = strcat(figout,'cpol.fig');
            savefig(fig_name)
            
            %% Plot investrment policy function.
            
            figure(4)
            
            plot(par.kgrid,sol.i)
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$i_{t}$'},'Interpreter','latex') 
            title('Investment Policy Function')
            
            fig_name = strcat(figout,'ipol.fig');
            savefig(fig_name)
            
            
            %% Plot government policy function.
            
            figure(5)
            
            plot(par.kgrid,sol.g)
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$g_{t}$'},'Interpreter','latex') 
            title('Government Policy Function')
            
            fig_name = strcat(figout,'npol.fig');
            savefig(fig_name)
            
            %% Plot value function.
            
            figure(6)
            
            plot(par.kgrid,sol.v)
                xlabel({'$k_{t}$'},'Interpreter','latex')
                ylabel({'$v_t(k_t,A_t)$'},'Interpreter','latex') 
            title('Value Function')

            fig_name = strcat(figout,'vfun.fig');
            savefig(fig_name)
            
            %% Plot simulated output.

            tgrid = linspace(1,par.T,par.T);

            figure(7)

            plot(tgrid,sim.ysim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$y^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Output')

            fig_name = strcat(figout,'ysim.fig');
            savefig(fig_name)

            %% Plot simulated capital choice.

            figure(8)

            plot(tgrid,sim.ksim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$k^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Capital Choice')

            fig_name = strcat(figout,'ksim.fig');
            savefig(fig_name)

            %% Plot simulated consumption.

            figure(9)

            plot(tgrid,sim.csim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$c^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Consumption')

            fig_name = strcat(figout,'csim.fig');
            savefig(fig_name)

            %% Plot simulated investment.

            figure(10)

            plot(tgrid,sim.isim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$i^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Investment')

            fig_name = strcat(figout,'isim.fig');
            savefig(fig_name)


            %% Plot simulated labor supply.

            figure(11)

            plot(tgrid,sim.gsim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$g^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Government Spending')

            fig_name = strcat(figout,'nsim.fig');
            savefig(fig_name)

            %% Plot simulated utility.

            figure(12)

            plot(tgrid,sim.usim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$u^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Utility')

            fig_name = strcat(figout,'usim.fig');
            savefig(fig_name)

            %% Plot simulated productivity.

            figure(13)

            plot(tgrid,sim.Asim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$A^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Productivity')

            fig_name = strcat(figout,'Asim.fig');
            savefig(fig_name)

        end

        function [] = compare_real_vs_simulated(par, sim)
            %data = readtable("vn_data.csv");
            data = readtable('vn_data.csv', 'VariableNamingRule', 'preserve');


            time_real = data.Year;
            Y_real = data.("GDP per capita");
            i_real = data.("Private Investment");
            c_real = data.("Private Consumption");

            time_sim = 1:par.T;
            Y_sim = sim.ysim;
            i_sim = sim.isim;
            c_sim = sim.csim;

            start_year = 1985;
            time_sim = (start_year:start_year + par.T - 1)';

            min_length = min([length(time_real), length(time_sim), length(Y_real), length(Y_sim)]);

            time_real = time_real(1:min_length);
            time_sim = time_sim(1:min_length);
            Y_real = Y_real(1:min_length);
            Y_sim = Y_sim(1:min_length);

            time_real = time_real(:);
            time_sim = time_sim(:);
            Y_real = Y_real(:);
            Y_sim = Y_sim(:);
            Y_real = log(Y_real);

            min_length2 = min([length(time_real), length(time_sim), length(i_real), length(i_sim)]);
            time_real = time_real(1:min_length2);
            time_sim = time_sim(1:min_length2);
            i_real = i_real(1:min_length2);
            i_sim = i_sim(1:min_length2);

            time_real = time_real(:);
            time_sim = time_sim(:);
            i_real = i_real(:);
            i_sim = i_sim(:);
            i_real = log(i_real + 1);

            min_length3 = min([length(time_real), length(time_sim), length(c_real), length(c_sim)]);
            time_real = time_real(1:min_length3);
            time_sim = time_sim(1:min_length3);
            c_real = c_real(1:min_length3);
            c_sim = c_sim(1:min_length3);

            time_real = time_real(:);
            time_sim = time_sim(:);
            c_real = c_real(:);
            c_sim = c_sim(:);
            c_real = log(c_real + 1);

            

            figure;
            plot(time_real, Y_real, 'r--', 'LineWidth', 2);
            hold on;
            plot(time_sim, Y_sim, 'b', 'LineWidth', 2);
            xlabel('Time'); ylabel('GDP');
            legend('Real Data', 'Simulated Data');
            title('GDP Comparison');
            grid on;

            figure;
            plot(time_real, i_real, 'r--', 'LineWidth', 2);
            hold on;
            plot(time_sim, i_sim, 'b', 'LineWidth', 2);
            xlabel('Time'); ylabel('Investment');
            legend('Real Data', 'Simulated Data');
            title('Investment Comparison');
            grid on;

            figure;
            plot(time_real, c_real, 'r--', 'LineWidth', 2);
            hold on;
            plot(time_sim, c_sim, 'b', 'LineWidth', 2);
            xlabel('Time'); ylabel('Consumption');
            legend('Real Data', 'Simulated Data');
            title('Consumption Comparison');
            grid on;

            valid_idx = (time_real >= 1995) & ~isnan(i_real) & ~isnan(c_real);
            i_real_stats = log(i_real(valid_idx));
            c_real_stats = log(c_real(valid_idx));

            stats = table();
            stats.Mean_Real = [mean(Y_real); mean(i_real_stats); mean(c_real_stats)];
            stats.Mean_Simulated = [mean(Y_sim); mean(i_sim); mean(c_sim)];
            stats.Std_Dev_Real = [std(Y_real); std(i_real_stats); std(c_real_stats)];
            stats.Std_Dev_Simulated = [std(Y_sim); std(i_sim); std(c_sim)];
            stats.Maximum_Real = [max(Y_real); max(i_real_stats); max(c_real_stats)];
            stats.Maximum_Simulated = [max(Y_sim); max(i_sim); max(c_sim)];
            stats.Minimum_Real = [min(Y_real); min(i_real_stats); min(c_real_stats)];
            stats.Minium_Simulated = [min(Y_sim); min(i_sim); min(c_sim)];
        
            % Display statistics
            var_names = {'GDP per capita', 'Investment', 'Consumption'};
            stats.Properties.RowNames = var_names;
            disp('--- Descriptive Statistics for Simulated vs Real Data ---');
            disp(stats);
        
        end
    end
end