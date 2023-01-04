%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: own_rep_check.m
% Author: Craig A. Chikis
% Date: 12/15/2021
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = own_rep_check()




	models = {@AA_QCU_environment_gen, @AkcigitAtes_environment_gen_1};
	rho_vec = power(1+linspace(5e-4, 0.065, 25),1/12)-1;
	g_vec = zeros(2, size(rho_vec, 2));
	for (ii = 1:size(g_vec, 1))

		m = models{ii}(); 

		m.alpha_fric = Inf;
		m.bridge_method = false;
		m.EIS = 1;
		m.r_guess = power(1.03,1/12)-1;
		m.two_adjust = 1;
		m.v_max_try = 1e4;

		m_vec = cell(1, size(g_vec, 2));
		for (jj = 1:size(g_vec, 2))
			m_iter = m; 
			m_iter.rho = rho_vec(jj); 
			m_vec{jj} = m_iter;
		end

		parfor (jj = 1:size(g_vec, 2))
			
			m_vec{jj} = workhorse_nested_robust(m_vec{jj});
			g_vec(ii, jj) = m_vec{jj}.growth_rate;


		end
	end
	save_this = [power(1+rho_vec,12)-1; power(1+g_vec,12)-1];
	save('Output/Store_Data/own_res_rep.mat', 'save_this')
end