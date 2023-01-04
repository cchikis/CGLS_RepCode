close all; clear all; n=50; xinit = zeros(1,2*n); xinit(n+1) = 1; flatpi=1;
sig=12; lamb = 1.21; rvec = [6:-0.1:0.1,0.05,0.02];
[pivec,~,~,~,~,lL,lF,l0] =compute_pi_fast(sig,lamb,n); 
pivec(n+1+flatpi:end) = pivec(n+1+flatpi); pivec(1:n+1-flatpi)=pivec(n+1-flatpi);
c=33.3569^2;
pi=pivec*c;
kap=3.9345;
pishr=pivec(n+1:-1:1)+pivec(n+1:end);

% The above code is identical to Lines 3-10 of LMS's script calibration_EMA_submit.m,
% except that compute_pi_fast.m here returns three additional scalar outputs:
% the production cost of a leader in a market with a productivity gap of 1 step (lL), 
% the production cost of a follower in such a market (lF), and 
% the production cost of a tied firm (l0).  See compute_pi_fast.m for details.  

%%%%%%%%%%%%%%%%%%%%%%%%
% Excercises with a discount rate decline from 3% to 1%
%%%%%%%%%%%%%%%%%%%%%%%%

% Solving for the BGP
r=3; dr=2;
xinit = zeros(1,2*n); xinit(n+1) = 1;
[xvec, muvec, ~, g,~] = gen_compute_eqm(lamb,pi,1,kap,r,xinit);
[xvec2, muvec2, ~, g2,~] = gen_compute_eqm(lamb,pi,1,kap,r-dr,xinit);
invcost=[xvec.^2/c,0]; invcost2=[xvec2.^2/c,0]; pishrvec=(pishr-invcost(n+1:end)-invcost(n+1:-1:1)); pishrvec2=(pishr-invcost2(n+1:end)-invcost2(n+1:-1:1));
LIst = pishrvec*muvec';
must = muvec*(0:n)';

% Lines 21-27 are identical to Lines 237-242 of LMS's main script (except
% that BGP productivity growth for r = 1 is returned as g2 in Line 24 of
% our script, so that we can use g2 later).  

LI_terminal = pishrvec2*muvec2';  % Profit share in the new BGP. Compare with Line 39 (profit share in initial BGP) of LMS's script.

eps=0.001; % As in LMS's code, a year is divided into periods of length epsilon years
horizon_in_years = 200;  % Calculate dynamics for this number of years after the decline in the discount rate
T=horizon_in_years / eps;  % Number of periods over which to calculate transition dynamics

%%%%%%%%%%%% Obtaining transition dynamics

% When PDA_leader = 0, sc18_transition.m calculates transition dynamics
% under the assumption that the Spillover from Followers PDA holds.  
% See Definition 2 of our paper, in Online Appendix B.  
% Under the Spillover from Followers PDA, a follower's productivity index 
% increases by 1 when the follower obtains an investment success or when
% technology diffusion occurs.  

% Obtain transition dynamics, fixing the time-scale error in LMS's code and
% using the Spillover from Follower PDA
PDA_leader = 0; 
percent_adjust = 0.01; % When percent_adjust = 0.01, sc18_transition calculates the evolution of the market state distribution fixing the time-scale error in LMS's code
[transg_fix100, transg_comp_fix100, transLI_fix100, ~,mu_fix100,ginit] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

% transg_fix100 (line 52) is the transition path for productivity growth, correcting the time-scale error in LMS but without including composition effects
% transg_comp_fix100 is the transition path for productivity growth, correcting the time-scale error and including composition effects
% The time-scale error in LMS's code is described in our paper's Section 2
% The composition effects omitted in LMS's code are described in our Section 3

% Next, obtain transition dynamics, using LMS transition dynamics for the market state
% distribution (i.e., not fixing time-scale error) and using the Spillover
% from Follower PDA
PDA_leader = 0; 
percent_adjust = 1; % When percent_adjust = 1, sc18_transition calculates the evolution of the market state distribution as in LMS's code, without fixing the time-scale error
[transg, ~, transLI, ~,mu] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

% Next, obtain transition dynamics, fixing the time-scale error in LMS's code,
% and using the Spillover from *Leader* PDA.  See Definition 3 of our paper, in Online Appendix B.
% Under the Spillover from Leader PDA, a leader's productivity index 
% increases by 1 when the leader obtain an investment success.  
% Results for the Spillover from Leader PDA are used only in Online Appendix B.
PDA_leader = 1;
percent_adjust = 0.01;
[~, transg_comp_fix100_PDAL, transLI_fix100_PDAL] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

%%%%%%%%%%%%%%%%%%%%%%%%
% Excercises with a discount rate decline from 3.59% to 0.33%
% (output stored as ``..._large'')
%%%%%%%%%%%%%%%%%%%%%%%%

r=3.59; dr=(3.59-0.33);
xinit = zeros(1,2*n); xinit(n+1) = 1;

% Solving for the BGP
[xvec_large, muvec_large, ~, g_large,~] = gen_compute_eqm(lamb,pi,1,kap,r,xinit);
[xvec2_large, muvec2_large, ~, g2_large,~] = gen_compute_eqm(lamb,pi,1,kap,r-dr,xvec_large);
invcost_large=[xvec_large.^2/c,0]; invcost2_large=[xvec2_large.^2/c,0]; 
pishrvec_large=(pishr-invcost_large(n+1:end)-invcost_large(n+1:-1:1));
pishrvec2_large=(pishr-invcost2_large(n+1:end)-invcost2_large(n+1:-1:1));
LIst_large = pishrvec_large*muvec_large';
must_large = muvec_large*(0:n)';
eps=0.001;
horizon_in_years = 200;
T=horizon_in_years / eps;  

% Obtain transition dynamics, fixing the time-scale error in LMS's code and
% using the Spillover from Follower PDA, for r = 3.59% to r = 0.33% exercise

PDA_leader = 0;
percent_adjust = 0.01;
[transg_fix100_large, transg_comp_fix100_large, transLI_fix100_large, ~,~,ginit_large] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2_large,percent_adjust,lL,lF,l0,sig,PDA_leader);

% Obtain transition dynamics, using LMS transition dynamics for state
% distribution (i.e., not fixing time-scale error) and using the Spillover
% from Follower PDA
PDA_leader = 0;
percent_adjust = 1;
[transg_large, ~, transLI_large, ~,~,~] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

[status,msg,msgID] = mkdir('figures_comment'); % Figures will be saved in figures_comment folder; create folder if missing

%%% Figure 1 of our main text

pre_horizon_in_years = horizon_in_years / 4;
Tst = pre_horizon_in_years / eps;
try close(fig1);catch;end
fig1=figure(1);

set(gcf, 'PaperUnits', 'inches');
    x_width=6.5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %   

xaxis = (1-Tst):T;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_fix100],'--r', 'LineWidth', 2); hold on;
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    ylim([0.75 1.5]); 
    xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    title({'Productivity growth'; 'over 800 quarters'}, 'Interpreter', 'latex');
    l1 = legend({"LMS Figure 8"; newline+ "Correcting" + newline + "time-scale" + newline + "error"}, ...
                 'Interpreter', 'latex', 'Location', 'northeast');
    set(l1, 'box', 'off')

zoom = 100;
Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (1-Tst):(T/zoom);
subplot(1,2,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg(1:T/zoom)],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_fix100(1:T/zoom)],'--r', 'LineWidth', 2); hold on;
xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); ylim([0.75 1.5]); %pbaspect([1.5 1 1]); ax.YGrid = 'on';
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    ylim([0.75 1.5]); 
    xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]);
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    title({'Productivity growth'; 'over 8 quarters'}, 'Interpreter', 'latex')

    saveas(gcf, "figures_comment/lms1.eps", 'epsc');  
   
%%% Figure 2 of our main text

try close(fig2);catch;end
Tst = pre_horizon_in_years / eps;
fig2=figure(2); 
set(gcf, 'PaperUnits', 'inches');
x_width=6.5;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
xaxis = (1-Tst):T;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,100*[LIst*ones(1,Tst),transLI_fix100],'--r', 'LineWidth', 2); hold on;
    yline(100*LI_terminal, 'LineStyle', ':', 'color', [0.4, 0.4, 0.4], 'LineWidth', 0.85);
    xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]);
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    ylim([-20, 20]); 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    title({'Profit share'; 'over 800 quarters'}, 'Interpreter', 'latex')

zoom = 100;
Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (1-Tst):(T/zoom);
subplot(1,2,2); plot(xaxis / (1/eps) * 4,100*[LIst*ones(1,Tst),transLI_fix100(1:T/zoom)],'--','LineWidth',2,'Color','r'); hold on;
    p3 = yline(100*LI_terminal, 'LineStyle', ':', 'color', [0.4, 0.4, 0.4], 'LineWidth', 0.85);
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    ylim([-20, 20]); 
    xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    title({'Profit share'; 'over 8 quarters'}, 'Interpreter', 'latex')

    saveas(gcf, "figures_comment/lms2.eps", 'epsc');

%%% Figure 3 of our main text

zoom = 10;
try close(fig3);catch;end
Tst = (pre_horizon_in_years / zoom) / eps;
fig3=figure(3); 
set(gcf, 'PaperUnits', 'inches'); 
x_width=6.5;
y_width=5.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

xaxis = (1-Tst):(T/zoom);

subplot(2,1,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg(1:T/zoom)],'-k', 'LineWidth', 1.75); hold on;
subplot(2,1,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_fix100(1:T/zoom)],'--r', 'LineWidth', 2); hold on;
subplot(2,1,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100(1:T/zoom)],':b', 'LineWidth', 2); hold on;
yline(ginit, 'LineStyle', ':', 'color', [0.4, 0.4, 0.4], 'LineWidth', 1.05)

xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]);   
ylim([0.65 2.3]); 
    title("Following a decline in the discount rate from 3\% to 1\%", 'Interpreter', 'latex')
    ax = gca;  
    ax.YGrid = 'off';
    box off; 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    l1 = legend({"Using LMS code"; newline+ "Correcting" + newline + "time-scale" + newline + "error only"; ...
                 newline + "Including composition"+ newline + "effects"}, ...
                'Interpreter', 'latex', 'Location', 'northeast');
    set(l1, 'box', 'off')

subplot(2,1,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_large(1:T/zoom)],'-k', 'LineWidth', 1.75); hold on;
subplot(2,1,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_fix100_large(1:T/zoom)],'--r', 'LineWidth', 2); hold on;
subplot(2,1,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100_large(1:T/zoom)],':b', 'LineWidth', 2); hold on;
yline(ginit, 'LineStyle', ':', 'color', [0.4, 0.4, 0.4], 'LineWidth', 1.05)
xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); 
    ylim([0.65 5]); 
    title("Following a decline in the discount rate from 3.59\% to 0.33\%", 'Interpreter' ,'latex')
    ax = gca;  
    ax.YGrid = 'off';
    box off; 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');

    saveas(gcf, "figures_comment/lms3.eps", 'epsc');

%%% Figure IA.3

zoom = 1;
try close(fig4);catch;end

Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (1-Tst):(T/zoom);
fig4=figure(4); 
set(gcf, 'PaperUnits', 'inches');
    x_width=6.5;
    y_width=3.5;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
xaxis = (1-Tst):(T/zoom);
plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100(1:T/zoom)], ':b', 'LineWidth', 2); hold on;
plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100_PDAL(1:T/zoom)],'LineStyle', '-.', 'color', [0, 0.5, 0], 'LineWidth', 2); hold on;
xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); ylim([0.75 4]); legend('Spillover from Followers PDA','Spillover from Leaders PDA');
xlabel('Quarters since $r$ declined','FontSize',12,'fontweight','normal','interpreter','latex'); 
    ylim([0.7 3.5]); 
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    l1 = legend({"Spillover from Followers PDA"; "Spillover from Leaders PDA"}, ...
                 'Interpreter', 'latex', 'Location', 'northeast');
    set(l1, 'box', 'off')

    saveas(gcf, "figures_comment/lms_IA3.eps", 'epsc');
   
  
%%%% Figure IA.1

try close(fig_IA1);catch;end
fig_IA1=figure(5); 
Tst = pre_horizon_in_years / eps;
xaxis = (1-Tst):T;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu_fix100],'--r', 'LineWidth', 2); hold on;
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    ylim([2 6.5]); 
    xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    %ylabel('Percent', 'Interpreter','latex');
    title({'Average leader-follower distance'; 'over 800 quarters'}, 'Interpreter', 'latex');
    l1 = legend({"LMS Figure 8"; newline+ "Correcting" + newline + "time-scale" + newline + "error"}, ...
                 'Interpreter', 'latex', 'Location', 'southeast');
    set(l1, 'box', 'off')
zoom = 100;
Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (1-Tst):(T/zoom);
subplot(1,2,2); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu(1:T/zoom)],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,2); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu_fix100(1:T/zoom)],'--r', 'LineWidth', 2); hold on;
xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); ylim([0.75 1.5]); %pbaspect([1.5 1 1]); ax.YGrid = 'on';
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    ylim([2 6.5]); 
    xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]);
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    %ylabel('Percent', 'Interpreter','latex');
    title({'Average leader-follower distance'; 'over 8 quarters'}, 'Interpreter', 'latex')

    saveas(gcf, "figures_comment/lms_IA1.eps", 'epsc');

%%% Figure IA.2

try close(fig_IA2);catch;end
fig_IA2=figure(6); 

set(gcf, 'PaperUnits', 'inches');
x_width=7.75;
y_width=3.75;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

subplot(1,3,1); plot(0:20,pivec(n+1:n+1+20),'-k', 'LineWidth', 1.75); hold on;
    plot(0:20,pivec(n+1:-1:n+1-20),'LineStyle', '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlabel('State, $s$', 'Interpreter', 'latex')
    title({"Flow profits"}, ...
           'Interpreter', 'latex')
    xlim([0, 20])
    l1 = legend('Leader', 'Follower', 'Interpreter', 'latex', 'Location', 'east');
    set(l1, 'box', 'off')

subplot(1,3,2);
plot(0:20,invcost(n+1:n+1+20),'-k', 'LineWidth', 1.75); hold on;
plot(0:20,invcost(n+1:-1:n+1-20),'LineStyle', '-', 'LineWidth', 2, 'color', [0.5, 0.5, 0.5]); 
plot(0:20,invcost2(n+1:n+1+20),'--k', 'LineWidth', 1.75); hold on;
    plot(0:20,invcost2(n+1:-1:n+1-20),'LineStyle', '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2); hold on;
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlabel('State, $s$', 'Interpreter', 'latex')
    title('Investment cost', 'Interpreter', 'latex')
    xlim([0, 20])
    l1=  legend(strcat("Leader, $r = 3\%$ "),...
           strcat("Follower, $r = 3\%$ "), ...
           strcat("Leader, $r = 1\%$ "), ...
           strcat("Follower, $r = 1\%$ "), ...
           'Interpreter', 'latex', 'Location', 'northeast');
    set(l1, 'box', 'off')

subplot(1,3,3); plot(0:20,muvec(1:20+1), '-k', 'LineWidth', 1.75); hold on;
    plot(0:20, muvec2(1:20+1), 'color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'LineStyle', '--');
    xlabel('State, $s$', 'Interpreter', 'latex')
    xlim([0, 20])
    title({'Stationary distribution'; 'of distance, $\{\mu_s\}_{s=0}^{\bar{s}}$'}, 'Interpreter', 'latex')
    l1 = legend(strcat("$r = 3\%$ "),...% num2str(100*(power(1+rho_vec(1),12)-1)), "\%"), ...
           strcat("$r = 1\%$ "), ...%num2str(100*(power(1+rho_vec(2),12)-1)), "\%"), ...
           'Interpreter', 'latex');
    set(l1, 'box', 'off')

    saveas(gcf, "figures_comment/lms_IA2.eps", 'epsc');

%%% Figure IA.4

zoom = 1;
try close(fig_IA3);catch;end
Tst = (pre_horizon_in_years / zoom) / eps;
fig_IA3=figure(7); 
set(gcf, 'PaperUnits', 'inches'); 
x_width=6.5;
y_width=5.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

xaxis = (1-Tst):(T/zoom);

subplot(2,1,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg(1:T/zoom)],'-k', 'LineWidth', 1.75); hold on;
subplot(2,1,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_fix100(1:T/zoom)],'--r', 'LineWidth', 2); hold on;
subplot(2,1,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100(1:T/zoom)],':b', 'LineWidth', 2); hold on;
yline(ginit, 'LineStyle', ':', 'color', [0.4, 0.4, 0.4], 'LineWidth', 1.05)
ylim([0.7 2]); 
xlim([-200, 800]);
title("Following a decline in the discount rate from 3\% to 1\%", 'Interpreter', 'latex')
ax = gca;  
ax.YGrid = 'on';
grid off;
box off; 
xlabel('Quarters since $r$ declined','Interpreter','latex'); 
ylabel('Percent', 'Interpreter','latex');
l1 = legend({"Using LMS code"; newline+ "Correcting" + newline + "time-scale" + newline + "error only"; ...
             newline + "Including composition"+ newline + "effects"}, ...
            'Interpreter', 'latex', 'Location', 'northeast');
set(l1, 'box', 'off')

subplot(2,1,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_large(1:T/zoom)],'-k', 'LineWidth', 1.75); hold on;
subplot(2,1,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_fix100_large(1:T/zoom)],'--r', 'LineWidth', 2); hold on;
subplot(2,1,2); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100_large(1:T/zoom)],':b', 'LineWidth', 2); hold on;
yline(ginit, 'LineStyle', ':', 'color', [0.4, 0.4, 0.4], 'LineWidth', 1.05)
ylim([0.7 5]); 
xlim([-200, 800]);
title("Following a decline in the discount rate from 3.59\% to 0.33\%", 'Interpreter' ,'latex')
ax = gca;
grid off; 
box off; 
xlabel('Quarters since $r$ declined','Interpreter','latex'); 
ylabel('Percent', 'Interpreter','latex');

    saveas(gcf, "figures_comment/lms_IA4.eps", 'epsc');
    
    
