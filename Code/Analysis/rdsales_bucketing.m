%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: rdsales_bucketing.m
% Author: Craig A. Chikis
% Date: 08/03/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_array = rdsales_bucketing(kappa, lambda, s_max, firm_distributions, nu_s, subsidy, wage_share, B, gamma, ...
                                          investment, bucketing)

    [~, bucket_array, bucket_ret] =  RDsales_quant(kappa, lambda, s_max, firm_distributions, nu_s, 0.5, subsidy, ...
            wage_share, B, gamma, investment, bucketing);


    mean_uncond = nansum(bucket_array(4, :).*bucket_array(2, :)); 
    sd_uncond = nanstd(bucket_array(2, :)); 

    [~, idx] = sort(bucket_array(2, :)); 
    bucket_array = bucket_array(:, idx); 
    bucket_array(5, :) = cumsum(bucket_array(4, :));
    try
        idx1 = find(bucket_array(5, :) > 0.25, 1);
        idx2 = find(bucket_array(5, :) <= 0.25, 1, 'last'); 

        p25_uncond = 1*(bucket_array(2, idx2) + ...
                (0.25 - bucket_array(5, idx2))*(bucket_array(2, idx1) - bucket_array(2, idx2))/(bucket_array(5, idx1) - bucket_array(5, idx2)));


        idx1 = find(bucket_array(5, :) > 0.5, 1);
        idx2 = find(bucket_array(5, :) <= 0.5, 1, 'last'); 

        p50_uncond = 1*(bucket_array(2, idx2) + ...
                (0.5 - bucket_array(5, idx2))*(bucket_array(2, idx1) - bucket_array(2, idx2))/(bucket_array(5, idx1) - bucket_array(5, idx2)));

        idx1 = find(bucket_array(5, :) > 0.75, 1);
        idx2 = find(bucket_array(5, :) <= 0.75, 1, 'last'); 

        p75_uncond = 1*(bucket_array(2, idx2) + ...
                (0.75 - bucket_array(5, idx2))*(bucket_array(2, idx1) - bucket_array(2, idx2))/(bucket_array(5, idx1) - bucket_array(5, idx2)));

    catch
        p25_uncond = bucket_array(2, find(bucket_array(5, :) > 0.25, 1));
        p50_uncond = bucket_array(2, find(bucket_array(5, :) > 0.5, 1));
        p75_uncond = bucket_array(2, find(bucket_array(5, :) > 0.75, 1));
    end




    q1_data = bucket_ret{1};
    q2_data = bucket_ret{2};
    q3_data = bucket_ret{3};
    q4_data = bucket_ret{4};
    q5_data = bucket_ret{5};


    q1_data(6, :) = q1_data(4, :)/sum(q1_data(4, :));
    [~, idx] = sort(q1_data(2, :)); 
    q1_data = q1_data(:, idx);
    q1_data(7, :) = cumsum(q1_data(6, :));

    q1_mean = nansum(q1_data(6,:).*q1_data(2, :))*100;
    q1_std = 100*nanstd(q1_data(2, :), q1_data(6, :));
    try
        idx1 = find(q1_data(7, :) > 0.25, 1);
        idx2 = find(q1_data(7, :) <= 0.25, 1, 'last'); 

        q1_p25 = 100*(q1_data(2, idx2) + (0.25 - q1_data(7, idx2))*(q1_data(2, idx1) - q1_data(2, idx2))/(q1_data(7, idx1) - q1_data(7, idx2)));


        idx1 = find(q1_data(7, :) > 0.5, 1);
        idx2 = find(q1_data(7, :) <= 0.5, 1, 'last'); 

        q1_p50 = 100*(q1_data(2, idx2) + (0.5 - q1_data(7, idx2))*(q1_data(2, idx1) - q1_data(2, idx2))/(q1_data(7, idx1) - q1_data(7, idx2)));

        idx1 = find(q1_data(7, :) > 0.75, 1);
        idx2 = find(q1_data(7, :) <= 0.75, 1, 'last'); 

        q1_p75 = 100*(q1_data(2, idx2) + (0.75 - q1_data(7, idx2))*(q1_data(2, idx1) - q1_data(2, idx2))/(q1_data(7, idx1) - q1_data(7, idx2)));

    catch
        q1_p25 = 100*q1_data(2, find(q1_data(7, :) >= 0.25, 1));
        q1_p50 = 100*q1_data(2, find(q1_data(7, :) >= 0.5, 1));
        q1_p75 = 100*q1_data(2, find(q1_data(7, :) >= 0.75, 1));
    end



    q2_data(6, :) = q2_data(4, :)/sum(q2_data(4, :));
    [~, idx] = sort(q2_data(2, :)); 
    q2_data = q2_data(:, idx);
    q2_data(7, :) = cumsum(q2_data(6, :));

    q2_mean = nansum(q2_data(6,:).*q2_data(2, :))*100;
    q2_std = 100*nanstd(q2_data(2, :), q2_data(6, :));
    try
        idx1 = find(q2_data(7, :) > 0.25, 1);
        idx2 = find(q2_data(7, :) <= 0.25, 1, 'last'); 

        q2_p25 = 100*(q2_data(2, idx2) + (0.25 - q2_data(7, idx2))*(q2_data(2, idx1) - q2_data(2, idx2))/(q2_data(7, idx1) - q2_data(7, idx2)));


        idx1 = find(q2_data(7, :) > 0.5, 1);
        idx2 = find(q2_data(7, :) <= 0.5, 1, 'last'); 

        q2_p50 = 100*(q2_data(2, idx2) + (0.5 - q2_data(7, idx2))*(q2_data(2, idx1) - q2_data(2, idx2))/(q2_data(7, idx1) - q2_data(7, idx2)));

        idx1 = find(q2_data(7, :) > 0.75, 1);
        idx2 = find(q2_data(7, :) <= 0.75, 1, 'last'); 

        q2_p75 = 100*(q2_data(2, idx2) + (0.75 - q2_data(7, idx2))*(q2_data(2, idx1) - q2_data(2, idx2))/(q2_data(7, idx1) - q2_data(7, idx2)));

    catch
        q2_p25 = 100*q2_data(2, find(q2_data(7, :) >= 0.25, 1));
        q2_p50 = 100*q2_data(2, find(q2_data(7, :) >= 0.5, 1));
        q2_p75 = 100*q2_data(2, find(q2_data(7, :) >= 0.75, 1));
    end



    q3_data(6, :) = q3_data(4, :)/sum(q3_data(4, :));
    [~, idx] = sort(q3_data(2, :)); 
    q3_data = q3_data(:, idx);
    q3_data(7, :) = cumsum(q3_data(6, :));

    q3_mean = nansum(q3_data(6,:).*q3_data(2, :))*100;
    q3_std = 100*nanstd(q3_data(2, :), q3_data(6, :));
    try
        idx1 = find(q3_data(7, :) > 0.25, 1);
        idx2 = find(q3_data(7, :) <= 0.25, 1, 'last'); 

        q3_p25 = 100*(q3_data(2, idx2) + (0.25 - q3_data(7, idx2))*(q3_data(2, idx1) - q3_data(2, idx2))/(q3_data(7, idx1) - q3_data(7, idx2)));


        idx1 = find(q3_data(7, :) > 0.5, 1);
        idx2 = find(q3_data(7, :) <= 0.5, 1, 'last'); 

        q3_p50 = 100*(q3_data(2, idx2) + (0.5 - q3_data(7, idx2))*(q3_data(2, idx1) - q3_data(2, idx2))/(q3_data(7, idx1) - q3_data(7, idx2)));

        idx1 = find(q3_data(7, :) > 0.75, 1);
        idx2 = find(q3_data(7, :) <= 0.75, 1, 'last'); 

        q3_p75 = 100*(q3_data(2, idx2) + (0.75 - q3_data(7, idx2))*(q3_data(2, idx1) - q3_data(2, idx2))/(q3_data(7, idx1) - q3_data(7, idx2)));

    catch
        q3_p25 = 100*q3_data(2, find(q3_data(7, :) >= 0.25, 1));
        q3_p50 = 100*q3_data(2, find(q3_data(7, :) >= 0.5, 1));
        q3_p75 = 100*q3_data(2, find(q3_data(7, :) >= 0.75, 1));
    end

    q4_data(6, :) = q4_data(4, :)/sum(q4_data(4, :));
    [~, idx] = sort(q4_data(2, :)); 
    q4_data = q4_data(:, idx);
    q4_data(7, :) = cumsum(q4_data(6, :));

    q4_mean = nansum(q4_data(6,:).*q4_data(2, :))*100;
    q4_std = 100*nanstd(q4_data(2, :), q4_data(6, :));
    try
        idx1 = find(q4_data(7, :) > 0.25, 1);
        idx2 = find(q4_data(7, :) <= 0.25, 1, 'last'); 

        q4_p25 = 100*(q4_data(2, idx2) + (0.25 - q4_data(7, idx2))*(q4_data(2, idx1) - q4_data(2, idx2))/(q4_data(7, idx1) - q4_data(7, idx2)));


        idx1 = find(q4_data(7, :) > 0.5, 1);
        idx2 = find(q4_data(7, :) <= 0.5, 1, 'last'); 

        q4_p50 = 100*(q4_data(2, idx2) + (0.5 - q4_data(7, idx2))*(q4_data(2, idx1) - q4_data(2, idx2))/(q4_data(7, idx1) - q4_data(7, idx2)));

        idx1 = find(q4_data(7, :) > 0.75, 1);
        idx2 = find(q4_data(7, :) <= 0.75, 1, 'last'); 

        q4_p75 = 100*(q4_data(2, idx2) + (0.75 - q4_data(7, idx2))*(q4_data(2, idx1) - q4_data(2, idx2))/(q4_data(7, idx1) - q4_data(7, idx2)));

    catch
        q4_p25 = 100*q4_data(2, find(q4_data(7, :) >= 0.25, 1));
        q4_p50 = 100*q4_data(2, find(q4_data(7, :) >= 0.5, 1));
        q4_p75 = 100*q4_data(2, find(q4_data(7, :) >= 0.75, 1));
    end


    q5_data(6, :) = q5_data(4, :)/sum(q5_data(4, :));
    [~, idx] = sort(q5_data(2, :)); 
    q5_data = q5_data(:, idx);
    q5_data(7, :) = cumsum(q5_data(6, :));

    q5_mean = nansum(q5_data(6,:).*q5_data(2, :))*100;
    q5_std = nanstd(q5_data(2, :), q5_data(6, :))*100;
    try
        idx1 = find(q5_data(7, :) > 0.25, 1);
        idx2 = find(q5_data(7, :) <= 0.25, 1, 'last'); 

        q5_p25 = 100*(q5_data(2, idx2) + (0.25 - q5_data(7, idx2))*(q5_data(2, idx1) - q5_data(2, idx2))/(q5_data(7, idx1) - q5_data(7, idx2)));


        idx1 = find(q5_data(7, :) > 0.5, 1);
        idx2 = find(q5_data(7, :) <= 0.5, 1, 'last'); 

        q5_p50 = 100*(q5_data(2, idx2) + (0.5 - q5_data(7, idx2))*(q5_data(2, idx1) - q5_data(2, idx2))/(q5_data(7, idx1) - q5_data(7, idx2)));

        idx1 = find(q5_data(7, :) > 0.75, 1);
        idx2 = find(q5_data(7, :) <= 0.75, 1, 'last'); 

        q5_p75 = 100*(q5_data(2, idx2) + (0.75 - q5_data(7, idx2))*(q5_data(2, idx1) - q5_data(2, idx2))/(q5_data(7, idx1) - q5_data(7, idx2)));

    catch
        q5_p25 = 100*q5_data(2, find(q5_data(7, :) >= 0.25, 1));
        q5_p50 = 100*q5_data(2, find(q5_data(7, :) >= 0.5, 1));
        q5_p75 = 100*q5_data(2, find(q5_data(7, :) >= 0.75, 1));
    end

    output_array = [q1_mean, q1_std, q1_p25, q1_p50, q1_p75;
                    q2_mean, q2_std, q2_p25, q2_p50, q2_p75;
                    q3_mean, q3_std, q3_p25, q3_p50, q3_p75;
                    q4_mean, q4_std, q4_p25, q4_p50, q4_p75;
                    q5_mean, q5_std, q5_p25, q5_p50, q5_p75;
                    100*[mean_uncond, sd_uncond, p25_uncond, p50_uncond, p75_uncond]]; 

end