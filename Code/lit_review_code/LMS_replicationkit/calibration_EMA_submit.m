function [] = calibration_EMA_submit()
    addpath(genpath('./export'));
    close all; clear all; n=50; xinit = zeros(1,2*n); xinit(n+1) = 1; flatpi=1;
    sig=12; lamb = 1.21; rvec = [6:-0.1:0.1,0.05,0.02];
    pivec=compute_pi_fast(sig,lamb,n);
    pivec(n+1+flatpi:end) = pivec(n+1+flatpi); pivec(1:n+1-flatpi)=pivec(n+1-flatpi);
    c=33.3569^2;
    pi=pivec*c;
    kap=3.9345;
    pishr=pivec(n+1:-1:1)+pivec(n+1:end);

    rstar = find(abs(0.3 - rvec) == min(abs(0.3 - rvec)), 1);
    x_coll = zeros(length(rvec), 2*n+1);
    mu_coll = zeros(length(rvec), n+1);
    % tic;
    for i=1:length(rvec)
        [xvec, muvec, ~, gvec(i),~,flag] = gen_compute_eqm(lamb,pi,1,kap,rvec(i),xinit);
        x_coll(i, :) = [xvec, 0];
        mu_coll(i, :) = muvec;
        if flag>0;xinit=xvec;end
        invcost=[xvec.^2/c,0];
        avgmu(i)=(0:n)*muvec';
        LI(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    end%;
    % toc;

    LMS_bgp.investment = x_coll(rstar, :);
    LMS_bgp.firm_distributions = mu_coll(rstar, :);
    LMS_bgp.rho = rvec(rstar);
    vtmp =  zeros(1,2*n); for j=1:2*n; vtmp(j+1) = vtmp(j)+LMS_bgp.investment(j); end;
    LMS_bgp.v = vtmp;
    LMS_bgp.rvec = rvec;
    LMS_bgp.gvec = gvec;

    save('../../../Output/Store_Data/lms_export.mat', 'LMS_bgp');
    close all
    % black = [0 0 0]; grey = [.7 .7 .7];

    % %% Figure 5A: growth
    % figure; ax=axes; plot(rvec,gvec,'-','LineWidth',3,'Color',black); xtickformat(ax, 'percentage'); ytickformat(ax, 'percentage'); 
    % ylim([0.68 1.12]); pbaspect([3 1 1]); ax.YGrid = 'on';
    % box off; set(ax,'FontSize', 20, 'ytick', [0.7 0.8 0.9 1 1.1]); 
    % xl=xlabel('growth-adjusted interest rate','FontSize',28,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 1.85; xl.Position(2)=xl.Position(2)+0.042;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % xline(3.6,'linestyle','--','linewidth',2)
    % xline(0.3,'linestyle','--','linewidth',2)
    % annotation('textarrow',[0.65 0.6],[0.54 0.62],'String','1984-2000','fontsize',20)
    % annotation('textarrow',[0.225 0.175],[0.53 0.435],'String','2001-2016','fontsize',20)
    % %export_fig('figures/fig1.pdf', '-pdf','-transparent');

    % %% Figure 5B: profit share
    % figure; ax=axes; plot(rvec,LI,'-','LineWidth',3,'Color',black); xtickformat(ax, 'percentage');
    % ylim([0.135 0.175]); pbaspect([3 1 1]); ax.YGrid = 'on'
    % box off; set(gca,'FontSize', 20, 'ytick', [0.14 0.15 0.16 0.17]); 
    % xl=xlabel('growth-adjusted interest rate','FontSize',28,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 1.85; xl.Position(2)=xl.Position(2)+0.004;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % xline(3.6,'linestyle','--','linewidth',2)
    % xline(0.3,'linestyle','--','linewidth',2)
    % %export_fig('figures/fig2.pdf', '-pdf','-transparent');

    % %% Figure 5C: average distance
    % figure; ax=axes; plot(rvec,avgmu,'-','LineWidth',3,'Color',black); xtickformat(ax, 'percentage');
    % %ylim([0.135 0.175]); 
    % pbaspect([3 1 1]); ax.YGrid = 'on'
    % box off; set(gca,'FontSize', 20, 'ytick', [0 5 10 15 20]); 
    % xl=xlabel('growth-adjusted interest rate','FontSize',28,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 1.85; xl.Position(2)=xl.Position(2)+1.7;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % xline(3.6,'linestyle','--','linewidth',2)
    % xline(0.3,'linestyle','--','linewidth',2)
    % %export_fig('figures/fig3.pdf', '-pdf','-transparent');

    % %% Figure 6A: value function at 4% and 2% interest rate
    % [x1,mu1]=gen_compute_eqm(lamb,pi,1,kap,3,xinit);
    % [x2,mu2]=gen_compute_eqm(lamb,pi,1,kap,1,xinit);
    % v = zeros(1,2*n); for j=1:2*n; v(j+1) = v(j)+x1(j); end; vL1 = v(n+1:2*n+1); vF1 = v(n+1:-1:1);
    % v = zeros(1,2*n); for j=1:2*n; v(j+1) = v(j)+x2(j); end; vL2 = v(n+1:2*n+1); vF2 = v(n+1:-1:1);

    % figure; hold off; s=12;
    % plot(0:s, vL1(1:s+1),'-','LineWidth',3,'Color',black); hold on; pbaspect([2.5 1 1]);
    % plot(0:s, vF1(1:s+1),'-','LineWidth',3,'Color',grey);
    % plot(0:s, vL2(1:s+1),'--','LineWidth',3,'Color',black);
    % plot(0:s, vF2(1:s+1),'--','LineWidth',3,'Color',grey);
    % h1=legend('$v_s$ (leader), $\hat{r}=4\%$','$v_{-s}$ (follower), $\hat{r}=4\%$', '$v_s$ (leader), $\hat{r}=2\%$', '$v_{-s}$ (follower), $\hat{r}=2\%$');
    % set(h1,'FontSize',22,'interpreter','latex','Location','east');
    % %h1.Position(1)=h1.Position(1)+0.05;
    % %title(['Value functions'],'interpreter','latex','FontSize',20);
    % xlabel('state $s$','interpreter','latex','FontSize',25);
    % set(gca, 'FontSize', 25,'ytick',[]);
    % ylim([0,135]); box off; legend boxoff;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % %export_fig('figures/fig4.pdf', '-pdf','-transparent');

    % %% Figure 6B: investment by state
    % figure; hold off; s=12;
    % plot(0:s, x1(n+1:n+s+1),'-','LineWidth',3,'Color',black); hold on; pbaspect([2.5 1 1]);
    % plot(0:s, x1(n+1:-1:n+1-s),'-','LineWidth',3,'Color',grey);
    % plot(0:s, x2(n+1:n+s+1),'--','LineWidth',3,'Color',black);
    % plot(0:s, x2(n+1:-1:n+1-s),'--','LineWidth',3,'Color',grey);
    % h1=legend('$\eta_s$ (leader), $\hat{r}=4\%$','$\eta_{-s}$ (follower), $\hat{r}=4\%$', '$\eta_s$ (leader), $\hat{r}=2\%$', '$\eta_{-s}$ (follower), $\hat{r}=2\%$');
    % %title(['\hspace{30mm} Top Panel',10,'Investment by productivity gap: leader ($\eta_s$) and follower ($\eta_{-s}$)'],'interpreter','latex','FontSize',20);
    % xlabel('state $s$','interpreter','latex','FontSize',25); box off; legend boxoff;
    % set(gca, 'FontSize', 25,'ytick',[]);
    % set(h1,'FontSize',22,'interpreter','latex','Location','northeast');
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % %export_fig('figures/fig5.pdf', '-pdf','-transparent');

    % %% Figure 6C: distribution of productivity gap
    % figure; hold off; s=12;
    % plot(0:s,mu1(1:s+1),'-','LineWidth',3,'Color',black); hold on; pbaspect([2.5 1 1]);
    % plot(0:s,mu2(1:s+1),'--','LineWidth',3,'Color',black);
    % h1=legend('$\hat{r}=4\%$','$\hat{r}=2\%$');
    % %title('Steady-state distribution of productivity gap: $\mu_s$','interpreter','latex','FontSize',25);
    % xlabel('state $s$','interpreter','latex','FontSize',25);
    % set(gca,'ytick',[])
    % %ylabel('$\mu_s$','interpreter','latex','fontsize',30);
    % set(h1,'FontSize',22,'interpreter','latex','Location','northeast');
    % set(gca, 'FontSize', 25);
    % axis([0 s 0 0.25]); box off; legend boxoff;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % %export_fig('figures/fig6.pdf', '-pdf','-transparent');

    % %% Policy counterfactuals for Figure 7 (section 6.1)
    % xinit = zeros(1,2*n); xinit(n+1) = 1;
    % picf1=pivec;
    % picf1(n+1+flatpi:end) = picf1(n+1+flatpi)*0.9;
    % c=33.3569^2;
    % pi=picf1*c;
    % kap=3.9345;
    % pishrcf1=picf1(n+1:-1:1)+picf1(n+1:end);
    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec_cf1(i),~,flag] = gen_compute_eqm(lamb,pi,1,kap,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     avgmu_cf1(i)=(0:n)*muvec';
    %     LI_cf1(i)=(pishrcf1-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % xinit = zeros(1,2*n); xinit(n+1) = 1;
    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec_cf2(i),~,flag] = gen_compute_eqm(lamb,pivec*c,1,kap*1.1,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     avgmu_cf2(i)=(0:n)*muvec';
    %     LI_cf2(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % xinit = zeros(1,2*n); xinit(n+1) = 1;
    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec_cf3(i),~,flag] = gen_compute_eqm(lamb,pi/0.9,1.1,kap,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     avgmu_cf3(i)=(0:n)*muvec';
    %     LI_cf3(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % %% Figure 7A: growth
    % figure; ax=axes; plot(rvec+gvec,gvec,'-','LineWidth',5,'Color',black); hold on; xtickformat(ax, 'percentage'); ytickformat(ax, 'percentage'); 
    % ylim([0.68 1.2]); pbaspect([1.5 1 1]); ax.YGrid = 'on';
    % xlim([0.85 6]);
    % box off; set(ax,'FontSize', 25, 'ytick', [0.7 0.8 0.9 1 1.1 1.2],'xtick',[1 2 3 4 5 6]); 
    % xl=xlabel('interest rate','FontSize',35,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 2; xl.Position(2)=xl.Position(2)+0.03;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % plot(rvec+gvec_cf1,gvec_cf1,'-.','markersize',10,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(rvec+gvec_cf3,gvec_cf3,':','markersize',10,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(rvec+gvec_cf2,gvec_cf2,'-*','markersize',8,'LineWidth',4,'Color',[1 1 1]*0.5);
    % %xline(4.69,'linestyle','--','linewidth',2)
    % %xline(1.09,'linestyle','--','linewidth',2)
    % h1=legend('baseline','10\% tax on leader profits', '10\% subsidy to follower investment', '10\% higher $\kappa$');
    % set(h1,'FontSize',35,'Location','southeast','interpreter','latex');
    % %export_fig('figures/fig_cf1.pdf', '-pdf','-transparent');

    % %% Figure 7B: LI
    % figure; ax=axes; plot(rvec+gvec,LI,'-','LineWidth',5,'Color',black); hold on; xtickformat(ax, 'percentage');
    % ylim([0.125 0.175]); pbaspect([1.5 1 1]); ax.YGrid = 'on';
    % xlim([0.85 6]);
    % box off; set(gca,'FontSize', 25, 'ytick', [0.13 0.14 0.15 0.16 0.17],'xtick',[1 2 3 4 5 6]); 
    % xl=xlabel('interest rate','FontSize',35,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 2; xl.Position(2)=xl.Position(2)+0.003;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % plot(rvec+gvec_cf1,LI_cf1,'-.','markersize',8,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(rvec+gvec_cf3,LI_cf3,':','markersize',8,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(rvec+gvec_cf2,LI_cf2,'-*','markersize',8,'LineWidth',4,'Color',[1 1 1]*0.5);
    % %xline(4.69,'linestyle','--','linewidth',2)
    % %xline(1.09,'linestyle','--','linewidth',2)
    % h1=legend('baseline','10\% tax on leader profits', '10\% subsidy to follower investment', '10\% higher $\kappa$');
    % set(h1,'FontSize',35,'Location','northeast','interpreter','latex');
    % %export_fig('figures/fig_cf2.pdf', '-pdf','-transparent');

    % %% Policy counterfactuals for Figure 8: state dependent kappa (section 6.1)

    % xinit = zeros(1,2*n); xinit(n+1) = 1;
    % pivec=compute_pi_fast(sig,lamb,n); 
    % pivec(n+1+flatpi:end) = pivec(n+1+flatpi); pivec(1:n+1-flatpi)=pivec(n+1-flatpi);
    % c=33.3569^2; pi=pivec*c;
    % rvec = [6:-0.1:0.1,0.05,0.02];

    % len=5;cap=1.5;kapvec1=kap*[linspace(cap,1,len),ones(1,n-len)];
    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec_ref1(i),~,flag] = gen_compute_eqm_state_dep_kappa(lamb,pi,1,kapvec1,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     LI_ref1(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % len=5;cap=2;kapvec2=kap*[linspace(cap,1,len),ones(1,n-len)];
    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec_ref2(i),~,flag] = gen_compute_eqm_state_dep_kappa(lamb,pi,1,kapvec2,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     LI_ref2(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % len=10;cap=2;kapvec3=kap*[linspace(cap,1,len),ones(1,n-len)];
    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec_ref3(i),~,flag] = gen_compute_eqm_state_dep_kappa(lamb,pi,1,kapvec3,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     LI_ref3(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % %%
    % % Figure 8B: growth
    % figure; ax=axes; plot(rvec+gvec,gvec,'-','LineWidth',5,'Color',black); hold on; xtickformat(ax, 'percentage'); ytickformat(ax, 'percentage'); 
    % ylim([0.68 1.65]); pbaspect([1.5 1 1]); ax.YGrid = 'on';
    % xlim([0.85 6]);
    % box off; set(ax,'FontSize', 25, 'ytick', [0.8 1 1.2 1.4 1.6],'xtick',[1 2 3 4 5 6]); 
    % xl=xlabel('interest rate','FontSize',35,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 2; xl.Position(2)=xl.Position(2)+0.03*(0.97/0.52);
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % plot(rvec+gvec_ref1,gvec_ref1,'-.','markersize',10,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(rvec+gvec_ref2,gvec_ref2,':','markersize',10,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(rvec+gvec_ref3,gvec_ref3,'-*','markersize',8,'LineWidth',4,'Color',[1 1 1]*0.5);
    % h1=legend('baseline','state-dependent $\kappa$: spec 1','state-dependent $\kappa$: spec 2','state-dependent $\kappa$: spec 3');
    % set(h1,'FontSize',35,'Location','southeast','interpreter','latex');
    % %export_fig('figures/fig_ref1.pdf', '-pdf','-transparent');

    % %%
    % % Figure 8A: kap_s
    % plotvec = 1:20;
    % figure; ax=axes; plot(plotvec,kap*ones(1,plotvec(end)),'-','LineWidth',5,'Color',black); hold on;
    % pbaspect([1.5 1 1]); ax.YGrid = 'on';
    % ylim([3.5 8]);
    % box off; set(ax,'FontSize', 25, 'ytick', [4 5 6 7 8],'xtick',[1 5 10 15 20]); 
    % xl=xlabel('state','FontSize',35,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 9; xl.Position(2)=xl.Position(2)+0.3;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % plot(plotvec,kapvec1(plotvec),'-.','LineWidth',4,'Color',[1 1 1]*0.5); hold on;
    % plot(plotvec,kapvec2(plotvec),':','markersize',10,'LineWidth',4,'Color',[1 1 1]*0.5);
    % plot(plotvec,kapvec3(plotvec),'-*','markersize',10,'LineWidth',4,'Color',[1 1 1]*0.5);
    % h1=legend('baseline','state-dependent $\kappa$: spec 1','state-dependent $\kappa$: spec 2','state-dependent $\kappa$: spec 3');
    % set(h1,'FontSize',35,'Location','northeast','interpreter','latex');
    % %export_fig('figures/fig_ref3.pdf', '-pdf','-transparent');



    % %% transitional dynamics
    % r=3; dr=2;
    % [xvec, muvec, ~, g,~] = gen_compute_eqm(lamb,pi,1,kap,r,xinit);
    % [xvec2, muvec2, ~, ~,~] = gen_compute_eqm(lamb,pi,1,kap,r-dr,xinit);
    % invcost=[xvec.^2/c,0]; invcost2=[xvec2.^2/c,0]; pishrvec=(pishr-invcost(n+1:end)-invcost(n+1:-1:1)); pishrvec2=(pishr-invcost2(n+1:end)-invcost2(n+1:-1:1));
    % LIst = pishrvec*muvec';
    % must = muvec*(0:n)';
    % T=2000; eps=0.001; %lambda=1.2; kappa=0.05;
    % [transg, transLI, vLovF,mu,ginit,vLovFinit,vL1,vL2,vF1,vF2] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2);

    % %% impulse response: growth (narrow, paper version)
    % try close(fig1);catch;end
    % fig1=figure(1); ax=axes; Tst=500; xaxis = (1-Tst):T;
    % plot(xaxis /250,[ginit*ones(1,Tst),transg],'-','LineWidth',3,'Color',black); hold on;
    % ylim([0.75 1.5]); pbaspect([1.5 1 1]); ax.YGrid = 'on';
    % box off; set(ax,'FontSize', 20, 'ytick', [0.8 1 1.2 1.4 1.6]); 
    % set(ax,'FontSize', 20, 'xtick', [-2 0 2 4 6 8]); 
    % xl=xlabel('Quarters since $r$ declined','FontSize',28,'fontweight','normal','interpreter','latex'); 
    % xl.Position(1) = xl.Position(1) + 3.3; xl.Position(2)=xl.Position(2)+0.03;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % %exportgraphics(gcf,'figures/fig_ipr_growth_narrow.pdf', 'backgroundcolor','none');

    % %% impulse response: distance(narrow, paper version)
    % try close(fig1);catch;end
    % fig1=figure(1); ax=axes; Tst=500; xaxis = (1-Tst):T;
    % plot(xaxis /250,[must*ones(1,Tst),mu],'-','LineWidth',3,'Color',black); hold on;
    % %axis([1-Tst T gst*0.6 transg(1)*1.1])
    % %xtickformat(ax, 'percentage'); ytickformat(ax, 'percentage'); 
    % ylim([1 7]); pbaspect([1.5 1 1]); ax.YGrid = 'on';
    % box off; set(ax,'FontSize', 20, 'ytick', [2 4 6]); 
    % set(ax,'FontSize', 20, 'xtick', [-2 0 2 4 6 8]); 
    % xl=xlabel('Quarters since $r$ declined','FontSize',28,'fontweight','normal','interpreter','latex'); 
    % xl.Position(1) = xl.Position(1) + 3.3; xl.Position(2)=xl.Position(2)+0.22;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % %exportgraphics(gcf,'figures/fig_ipr_dist_narrow.pdf', 'backgroundcolor','none');


    % %% Figure 3: illustration of g against r
    % close all; clear all; n=50; xinit = zeros(1,2*n); xinit(n+1) = 1; flatpi=1;
    % sig=12; lamb = 1.21; rvec = [10:-0.1:0.1,0.05,0.02];
    % pivec = 1*(1-lamb.^(-1*(1:n)));pivec(flatpi:end) = pivec(flatpi);pivec=[zeros(1,n+1),pivec];c=33.3569^2;
    % pi=pivec*c;
    % kap=3.9345;
    % pishr=pivec(n+1:-1:1)+pivec(n+1:end);

    % tic;for i=1:length(rvec)
    %     [xvec, muvec, ~, gvec(i),~,flag] = gen_compute_eqm(lamb,pi,1,kap,rvec(i),xinit);
    %     if flag>0;xinit=xvec;end
    %     invcost=[xvec.^2/c,0];
    %     avgmu(i)=(0:n)*muvec';
    %     LI(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
    % end;toc;

    % black = [0 0 0]; grey = [.7 .7 .7];

    % figure; ax=axes; plot(rvec,gvec,'-','LineWidth',3,'Color',black); xtickformat(ax, 'percentage'); ytickformat(ax, 'percentage'); 
    % ylim([0.68 1.05]); pbaspect([2 1 1]); ax.YGrid = 'on'
    % box off; set(ax,'FontSize', 20, 'ytick', [],'xtick',[]); 
    % xl=xlabel('interest rate $r$','interpreter','latex','FontSize',28,'fontweight','normal'); xl.Position(1) = xl.Position(1) + 4;
    % yl=ylabel('Productivity growth $g$', 'interpreter','latex','fontsize',28);yl.Position(2)=yl.Position(2)+0.07;
    % set(gcf,'OuterPosition', [100, 800, 1200, 1200]);
    % annotation('textarrow',[0.225 0.145],[0.37 0.385],'String','$\kappa\cdot\ln\lambda$','interpreter','latex','fontsize',30)
    % %export_fig('figures/fig0.pdf', '-pdf','-transparent');


end

