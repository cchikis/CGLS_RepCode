%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: RDsales_quant.m
% Author: Craig A. Chikis
% Date: 07/25/2022D
% Note(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p50_approx, markup_array, bucket_ret] =  RDsales_quant(kappa, lambda, s_max, firm_distributions, nu_s, quant, subsidy, ...
                                                                  wage_share, B, gamma, investment, bucketing)


    G = @(x, B, gamma) (x./B).^(1./gamma);
    sales_divisor = sales_func(kappa, s_max, nu_s);
    RD_sales = wage_share*G(investment, B, gamma).*(1-subsidy)./sales_divisor;
    weightL = nu_s.^(1-kappa)./(nu_s.^(1-kappa)+1);
    weightF = 1./(nu_s.^(1-kappa) + 1);

    RD_salesL = RD_sales((s_max+1):end);
    RD_salesF = flip(RD_sales(1:(s_max+1)));
	if (~(kappa >= 9999))
		markup_array = [-s_max:s_max; 
						flip(RD_salesF), RD_salesL(2:end); 
						-s_max:s_max; 
						flip(abs(firm_distributions(2:end)))/2, abs(firm_distributions(1)), abs(firm_distributions(2:end)/2)];
	else
		markup_array = [-s_max:s_max; 
						flip(RD_salesF), RD_salesL(2:end); 
						-s_max:s_max; 
						zeros(1, s_max), abs(firm_distributions(1)), abs(firm_distributions(2:end))];

		markup_array(isinf(markup_array)) = nan;
	end
    [~, sort_idx] = sort(markup_array(1,:));

	markup_array = markup_array(:, sort_idx);
	markup_array(5,:) = cumsum(markup_array(4,:));

	cutpoints = [quant] - markup_array(5, :);

	p50_loc = [find(cutpoints(1, :) > 0 ,1,'last'), find(cutpoints(1, :) < 0, 1,'first')];

	if (length(p50_loc) < 2)
		p50_loc(2) = p50_loc(1);
	end
	


	p50_approx = ((markup_array(1, p50_loc(2)) - markup_array(1, p50_loc(1)))/...
					(markup_array(5, p50_loc(2)) - markup_array(5, p50_loc(1))))*(quant - markup_array(5,p50_loc(1))) + ...
				markup_array(1, p50_loc(1));

	
	
	if (isnan(p50_approx))
		p50_approx = markup_array(1, p50_loc(1));
	end

	sales_dist = markup_array;

	if (strcmp(bucketing, "quintile"))
		flipsalesdist = flip(sales_dist, 2);
		flipsalesdist(5, :) = cumsum(flipsalesdist(4, :));
		cutoff1 = flipsalesdist(1, find(cumsum(flipsalesdist(4, :)) >= 0.2, 1));
		q1_data = sales_dist(:, sales_dist(1, :) >= cutoff1);
		if (sum(q1_data(4, :)) > 0.2)
			diff = sum(q1_data(4, :)) - 0.2;
			q1_data(4, 1) = q1_data(4, 1) - diff;
		end
		cutoff2 = flipsalesdist(1, find(cumsum(flipsalesdist(4, :)) >= 0.4, 1));
		q2_data = sales_dist(:, sales_dist(1, :) >= cutoff2 & sales_dist(1, :) <= q1_data(1, 1));

		lb1 = q1_data(4, q1_data(1, :) == q2_data(1, end));
		if (isempty(lb1))
			lb1 = 0;
		end

		q2_data(4, end) = q2_data(4, end) - lb1; 
		if (sum(q2_data(4, :)) > 0.2)
			diff = sum(q2_data(4, :)) - 0.2;
			q2_data(4, 1) = q2_data(4, 1) - diff;
		end
		cutoff3 = flipsalesdist(1, find(cumsum(flipsalesdist(4, :)) >= 0.6, 1));
		q3_data = sales_dist(:, sales_dist(1, :) >= cutoff3 & sales_dist(1, :) <= q2_data(1, 1));

		lb2 = q2_data(4, q2_data(1, :) == q3_data(1, end));
		if (isempty(lb2))
			lb2 = 0;
		end
		lb1 = q1_data(4, q1_data(1, :) == q3_data(1, end));
		if (isempty(lb1))
			lb1 = 0;
		end

		q3_data(4, end) = q3_data(4, end) - lb2 - lb1; 
							
		if (sum(q3_data(4, :)) > 0.2)
			diff = sum(q3_data(4, :)) - 0.2;
			q3_data(4, 1) = q3_data(4, 1) - diff;
		end
		cutoff4 = flipsalesdist(1, find(cumsum(flipsalesdist(4, :)) >= 0.8, 1));
		q4_data = sales_dist(:, sales_dist(1, :) >= cutoff4 & sales_dist(1, :) <= q3_data(1, 1));

		lb3 = q3_data(4, q3_data(1, :) == q4_data(1, end));
		if (isempty(lb3))
			lb3 = 0;
		end

		lb2 = q2_data(4, q2_data(1, :) == q4_data(1 ,end)); 
		if (isempty(lb2))
			lb2 = 0;
		end
		lb1 = q1_data(4, q1_data(1, :) == q4_data(1, end));
		if (isempty(lb1))
			lb1 = 0;
		end

		q4_data(4, end) = q4_data(4, end) - lb3 - lb2 - lb1;
		if (sum(q4_data(4, :)) > 0.2)
			diff = sum(q4_data(4, :)) - 0.2;
			q4_data(4, 1) = q4_data(4, 1) - diff;
		end
		cutoff4 = flipsalesdist(1, find(cumsum(flipsalesdist(4, :)) >= 0.8, 1));
		q4_data = sales_dist(:, sales_dist(1, :) >= cutoff4 & sales_dist(1, :) <= q3_data(1, 1));

		lb3 = q3_data(4, q3_data(1, :) == q4_data(1, end));
		if (isempty(lb3))
			lb3 = 0;
		end

		lb2 = q2_data(4, q2_data(1, :) == q4_data(1 ,end)); 
		if (isempty(lb2))
			lb2 = 0;
		end
		lb1 = q1_data(4, q1_data(1, :) == q4_data(1, end));
		if (isempty(lb1))
			lb1 = 0;
		end

		q4_data(4, end) = q4_data(4, end) - lb3 - lb2 - lb1;
		if (sum(q4_data(4, :)) > 0.2)
			diff = sum(q4_data(4, :)) - 0.2;
			q4_data(4, 1) = q4_data(4, 1) - diff;
		end
		q5_data = sales_dist(:, sales_dist(1, :) <= q4_data(1, 1));

		lb4 = q4_data(4, q4_data(1, :) == q5_data(1, end));
		if (isempty(lb4))
			lb4 = 0;
		end

		lb3 = q3_data(4, q3_data(1, :) == q5_data(1, end));
		if (isempty(lb3))
			lb3 = 0;
		end

		lb2 = q2_data(4, q2_data(1, :) == q5_data(1 ,end)); 
		if (isempty(lb2))
			lb2 = 0;
		end
		lb1 = q1_data(4, q1_data(1, :) == q5_data(1, end));
		if (isempty(lb1))
			lb1 = 0;
		end

		q5_data(4, end) = q5_data(4, end) - lb4 - lb3 - lb2 - lb1;
		
		bucket_ret= {q1_data, q2_data, q3_data, q4_data, q5_data}; 
	
	end
		
		
		