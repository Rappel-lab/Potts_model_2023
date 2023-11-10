%% Parameter Search Plots for Cancer cells in ECM %%
clear

% local vs global remodeling (symmetric % change)
% par1_list_FIXED = [0, 5, 12.5, 25, 37.5, 45, 50];       % global remodeling
% par1_list_FIXED = [0, 5, 10, 12.5, 25, 37.5, 40, 45, 50];
% par1_list_FIXED = [5,45];
% par2_list_FIXED = [0]; % local remodeling
% par1name = 'rho';
% par2name = 'KD1';
% repeatNN = 10;
% run_counter = 0;


% polarization - fold vs rate
par1_list_FIXED = [8, 16, 32, 64, 128, 256, 512, 1024, 2028];    % rate
par2_list_FIXED = 2.^(-9:1);    % fold
par1name = 'rate';
par2name = 'fold';
repeatNN = 1;
run_counter = 0;

% par1_list_FIXED = 10.^linspace(-1.5,1,14);
% par2_list_FIXED = 10.^linspace(-0.5,2,14);

% recent %
% par1_list = [1.1];        % spheroid
% par2_list = [1.2];
% par1_list = [0.45];      % network
% par2_list = [41];
% par1_list_FIXED = [0.45, 1.1];      % network
% par2_list_FIXED = [41, 1.2];

% par1name = 'rho';
% par2name = 'KDplus';
% repeatN_idx = '';
% repeatNN = 10;
% run_counter = 0;

save_2cell_png = 0; % save a png of the last 2 cell snapshot

[circularity_ALL, v1_avg_ALL, v2_avg_ALL, angv1_avg_ALL, angv2_avg_ALL, per_lin1_total_ALL, ...
per_lin2_total_ALL, per_ang1_total_ALL, per_ang2_total_ALL, ncell_ALL, cell2ECMratio_ALL, P1_ALL, P2_ALL, ... 
pro_len1_ALL, pro_len2_ALL, pro_maxlen1_ALL, pro_maxlen2_ALL, pro_lifetime1_ALL, pro_lifetime2_ALL, pro_count1_ALL, pro_count2_ALL] ...
= deal(nan([length(par1_list_FIXED),length(par2_list_FIXED),repeatNN],'double'));

%% Load data
for iii = 1:length(par1_list_FIXED)
    
    for jjj = 1:length(par2_list_FIXED)
%         if iii ~= jjj       %@%@%@    for individual case
%             continue        %@%@%@
%         end                 %@%@%@
        for rrr = 1:repeatNN
        
        clearvars -except iii jjj rrr par1* par2* repeatN* run_counter *_ALL save_2cell_png
        
%         par1_list = 10.^linspace(-1.5,1,14); %@%@%@
%         par2_list = 10.^linspace(-0.5,2,14); %@%@%@

        par1 = par1_list_FIXED(iii);
        par2 = par2_list_FIXED(jjj);
        load(sprintf('sim_%s%0.2e_%s%0.2e_%s_data.mat',par1name,par1,par2name,par2,num2str(rrr)))
        
        if exist('exitnote','var')
            continue
        end
        
        if length(circularity) > 1
            circularity_ALL(iii,jjj,rrr) = (4 * pi * sum(allAreas))./(sum(allPerimeters) .^ 2);
        else
            circularity_ALL(iii,jjj,rrr) = circularity;
        end
        v1_avg_ALL(iii,jjj,rrr) = v1_avg;
        v2_avg_ALL(iii,jjj,rrr) = v2_avg;
        angv1_avg_ALL(iii,jjj,rrr) = angv1_avg;
        angv2_avg_ALL(iii,jjj,rrr) = angv2_avg;
        per_lin1_total_ALL(iii,jjj,rrr) = per_lin1_total;
        per_lin2_total_ALL(iii,jjj,rrr) = per_lin2_total;
        per_ang1_total_ALL(iii,jjj,rrr) = per_ang1_total;
        per_ang2_total_ALL(iii,jjj,rrr) = per_ang2_total;
        ncell_ALL(iii,jjj,rrr) = p.ncell;
        %%
        % Protrusion %
        ang_thres = 2.5;
        len_thres = 4;
        time_thes = 0; % min
        Pangle_all = atan(Py_all./Px_all);
        Pangle_all(1,:) = [];
        P_all = sqrt(Px_all.^2 + Py_all.^2);

        if ~isempty(t_points)
            Pangle_all_2cell = Pangle_all(t_points(1):t_points(end),1:2);
            P_all_2cell = P_all(t_points(1):t_points(end),1:2);
        else
            Pangle_all_2cell = [nan, nan];
            P_all_2cell = [nan, nan];
        end
        [protrusion_lifetime1, protrusion_len1, protrusion_P_avg1] = GetProtrusion(Pangle_all_2cell(:,1),P_all_2cell(:,1),ang_thres);
        [protrusion_lifetime2, protrusion_len2, protrusion_P_avg2] = GetProtrusion(Pangle_all_2cell(:,2),P_all_2cell(:,2),ang_thres);

        valid1 = (protrusion_len1 >= len_thres) & (protrusion_lifetime1*dt/60 >= time_thes);
        valid2 = (protrusion_len2 >= len_thres) & (protrusion_lifetime2*dt/60 >= time_thes);

        protrusion_lifetime1(~valid1) = []; protrusion_len1(~valid1) = []; protrusion_P_avg1(~valid1) = [];
        protrusion_lifetime2(~valid2) = []; protrusion_len2(~valid2) = []; protrusion_P_avg2(~valid2) = [];

        protrusion_lifetime1_avg = mean(protrusion_lifetime1);
        protrusion_len1_avg = mean(protrusion_len1);
        protrusion_P_avg1_avg = mean(protrusion_P_avg1);
        protrusion_lifetime2_avg = mean(protrusion_lifetime2);
        protrusion_len2_avg = mean(protrusion_len2);
        protrusion_P_avg2_avg = mean(protrusion_P_avg2);
        %
        pro_len1_ALL(iii,jjj,rrr) = protrusion_len1_avg;
        pro_len2_ALL(iii,jjj,rrr) = protrusion_len2_avg;
        pro_lifetime1_ALL(iii,jjj,rrr) = protrusion_lifetime1_avg;
        pro_lifetime2_ALL(iii,jjj,rrr) = protrusion_lifetime2_avg;
        
        if ~isempty(protrusion_len1)
            pro_maxlen1_ALL(iii,jjj,rrr) = max(protrusion_len1);
            pro_count1_ALL(iii,jjj,rrr) = length(protrusion_len1);
        else
            pro_maxlen1_ALL(iii,jjj,rrr) = NaN;
            pro_count1_ALL(iii,jjj,rrr) = NaN;
        end
        
        if ~isempty(protrusion_len2)
            pro_maxlen2_ALL(iii,jjj,rrr) = max(protrusion_len2);
            pro_count2_ALL(iii,jjj,rrr) = length(protrusion_len2);
        else
            pro_maxlen2_ALL(iii,jjj,rrr) = NaN;
            pro_count2_ALL(iii,jjj,rrr) = NaN;
        end
        %%
        P1_ALL(iii,jjj,rrr) = P(1);
        if length(P) > 1
            P2_ALL(iii,jjj,rrr) = P(2);
        end
        
        if sum(nsig_2cell(:)) ~= 0
            cell2ECMratio_ALL(iii,jjj,rrr) = sum(nsig_2cell(:) > 0)/sum(S_2cell(:) == 0);
        else
            cell2ECMratio_ALL(iii,jjj,rrr) = sum(nsig(:) > 0)/sum(Snew(:) == 0);
        end

        if save_2cell_png
                FIG_2cell = figure('Position',[10 10 subplotN*512 512]);
                colormap turbo
            if sum(nsig_2cell(:)) ~= 0
                imagesc(nsig_2cell*floor(87/2)+int8((nsig_2cell>0)*30)+int8(S_2cell/Smax*10),[0 128])
                hold on 
                quiver([cmx_2cell],[cmy_2cell],[Px_2cell]*Pscale,[Py_2cell]*Pscale,'off','color','white','LineWidth',1)
                hold off
            else
                imagesc(nsig*floor(87/p.ncell)+int8((nsig>0)*30)+int8(Snew/Smax*10),[0 128])
                hold on 
                quiver([cmx],[cmy],[Px]*Pscale,[Py]*Pscale,'off','color','white','LineWidth',1)
                hold off
            end
            saveas(FIG_2cell,sprintf('sim_%s%0.2e_%s%0.2e_%s_2cell.png',par1name,par1,par2name,par2,num2str(rrr)))
            close(FIG_2cell)
        end
        run_counter = run_counter + 1;
        fprintf('%s:%0.2e, %s:%0.2e, repeat:%s [%d/%d]\n',par1name,par1,par2name,par2,num2str(rrr),run_counter,length(par1_list_FIXED)*length(par2_list_FIXED)*repeatN);
        end
    end
end

%% PART A
%% Average over repeats (do not average if for whisket plot)

success_run = repeatNN - sum(isnan(v1_avg_ALL),[3]);

circularity_ALL = nanmean(circularity_ALL,3);
v1_avg_ALL = nanmean(v1_avg_ALL,3);
v2_avg_ALL = nanmean(v2_avg_ALL,3);
angv1_avg_ALL = nanmean(angv1_avg_ALL,3);
angv2_avg_ALL = nanmean(angv2_avg_ALL,3);
per_lin1_total_ALL = nanmean(per_lin1_total_ALL,3);
per_lin2_total_ALL = nanmean(per_lin2_total_ALL,3);
per_ang1_total_ALL = nanmean(per_ang1_total_ALL,3);
per_ang2_total_ALL = nanmean(per_ang2_total_ALL,3);
ncell_ALL = nanmean(ncell_ALL,3);
cell2ECMratio_ALL = nanmean(cell2ECMratio_ALL,3);
        
%% plot name conversion
par1plotname = [];
par2plotname = [];

if strcmp(par1name,'rho')
    par1plotname = 'global degradation';
end

if strcmp(par2name,'KDplus')
    par2plotname = 'local degradation';
end

%% Label change
xtnew = [0,1,2];
xtlbl = 10.^xtnew;

ytnew = [-1,0,1];
ytlbl = 10.^ytnew;

%% Find boundary
bp = [];    % boundary parameter
x_newmin = -0.6;
x_newmax = 2.2;
y_newmin = -1.7;
y_newmax = 1.2;

% one cell %
[X, Y] = meshgrid(log10(par2_list_FIXED),log10(par1_list_FIXED));
onecell = ncell_ALL == 1;
bp.xvalid_onecell = X(onecell);
bp.yvalid_onecell = Y(onecell);
bp.xvalid_onecell(bp.xvalid_onecell == min(bp.xvalid_onecell)) = x_newmin;
bp.yvalid_onecell(bp.yvalid_onecell == min(bp.yvalid_onecell)) = y_newmin;
bp.k_onecell = boundary(bp.xvalid_onecell,bp.yvalid_onecell,0.01);

% rotational cell %
angv_avg = (angv1_avg_ALL + angv2_avg_ALL)/2;
angv_thres = nanmean(angv_avg(:));
rotcell = angv_avg >= angv_thres;
bp.xvalid_rotcell = X(rotcell);
bp.yvalid_rotcell = Y(rotcell);
bp.xvalid_rotcell(bp.xvalid_rotcell == min(bp.xvalid_rotcell)) = x_newmin;
bp.yvalid_rotcell(bp.yvalid_rotcell == max(bp.yvalid_rotcell)) = y_newmax;
bp.k_rotcell = boundary(bp.xvalid_rotcell,bp.yvalid_rotcell,0.01);

% unconstrained cell %
freecell = cell2ECMratio_ALL <= 0.8;
bp.xvalid_freecell = X(freecell);
bp.yvalid_freecell = Y(freecell);
bp.xvalid_freecell(bp.xvalid_freecell == min(bp.xvalid_freecell)) = x_newmin;
bp.xvalid_freecell(bp.xvalid_freecell == max(bp.xvalid_freecell)) = x_newmax;
bp.yvalid_freecell(bp.yvalid_freecell == max(bp.yvalid_freecell)) = y_newmax;
bp.k_freecell = boundary(bp.xvalid_freecell,bp.yvalid_freecell,0.5);

% adjust boundary slightly %
onecell_rotcell_ybound = (max(bp.yvalid_onecell) + min(bp.yvalid_rotcell))/2;
bp.yvalid_onecell(bp.yvalid_onecell == max(bp.yvalid_onecell)) = onecell_rotcell_ybound;
bp.yvalid_rotcell(bp.yvalid_rotcell == min(bp.yvalid_rotcell)) = onecell_rotcell_ybound;

%% Make plots
bpstatus = 0;

f_circularity = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), circularity_ALL, [0 1])
hold on
colorbar
title('Circularity')
cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)

saveas(f_circularity,'circularity_nobp.png')

% xtnew = linspace(min(par2_list), max(par2_list), 7);                      % New 'XTick' Values
% xtlbl = 10.^round(linspace(log10(min(par2_list)), log10(max(par2_list)), 7),2);                        % New 'XTickLabel' Vector
% set(gca, 'XTick', xtnew, 'XTickLabel', xtlbl)                               % Label Ticks
%%
f_v = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), (v1_avg_ALL + v2_avg_ALL)/2)
hold on
c = colorbar;
c.Label.String = '\mum / h';
title('Linear Velocity')
cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)

saveas(f_v,'linear_velocity.png')
%%
f_angv = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), angv_avg)
hold on
c = colorbar;
c.Label.String = 'rad / h';
title('Anuglar Velocity')
cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)

saveas(f_angv,'anuglar_velocity.png')

%%
f_per_lin = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), (per_lin1_total_ALL + per_lin2_total_ALL)/2, [0 1])
hold on
c = colorbar;
title('Linear Persistence Ratio')
cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)

saveas(f_per_lin,'linear_persistence_ratio.png')
%%
f_per_ang = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), (per_ang1_total_ALL + per_ang2_total_ALL)/2, [0 1])
hold on
c = colorbar;
title('Angular Persistence Ratio')
cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)

saveas(f_per_ang,'anugular_persistence_ratio.png')
%%
f_ncell = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), ncell_ALL)
hold on
c = colorbar;
title('Number of Cells at 60 Hours')

cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)
saveas(f_ncell,'ncell_60h.png')


%%
f_cell2ECM = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), cell2ECMratio_ALL)
hold on
c = colorbar;
title('Cell to ECM Volume Ratio')

cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)
saveas(f_cell2ECM,'cell2ECM_ratio.png')

%%
f_suc = figure;
imagesc(log10(par2_list_FIXED), log10(par1_list_FIXED), success_run)
hold on
c = colorbar;
title('Number of successful run included')

cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)
saveas(f_cell2ECM,'success_run.png')

%% PART B
%% Whisket Plot
groupname = ['Spheroid','Network'];
xname = {'Rotational', 'Invasive'};
markersz = 30;
Ndata = repeatNN;

% local vs global remodeling (symmetric % change)
spheroid_iii = 2;
spheroid_jjj = 1;

% previous
% spheroid_iii = 9;
% spheroid_jjj = 4;
% 
% spheroid_iii = 2;
% spheroid_jjj = 2;

spheroid_v1 = v1_avg_ALL(spheroid_iii,spheroid_jjj,:);
spheroid_v2 = v2_avg_ALL(spheroid_iii,spheroid_jjj,:);
spheroid_v = [spheroid_v1(:);spheroid_v2(:)];
spheroid_angv = reshape([angv1_avg_ALL(spheroid_iii,spheroid_jjj,:); angv2_avg_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_perlin = reshape([per_lin1_total_ALL(spheroid_iii,spheroid_jjj,:); per_lin2_total_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_perang = reshape([per_ang1_total_ALL(spheroid_iii,spheroid_jjj,:); per_ang2_total_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_P = reshape([P1_ALL(spheroid_iii,spheroid_jjj,:); P2_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_prolen = reshape([pro_len1_ALL(spheroid_iii,spheroid_jjj,:); pro_len2_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_promaxlen = reshape([pro_maxlen1_ALL(spheroid_iii,spheroid_jjj,:); pro_maxlen2_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_prolifetime = reshape([pro_lifetime1_ALL(spheroid_iii,spheroid_jjj,:); pro_lifetime2_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);
spheroid_procount = reshape([pro_count1_ALL(spheroid_iii,spheroid_jjj,:); pro_count2_ALL(spheroid_iii,spheroid_jjj,:)],[Ndata*2 1]);

% local vs global remodeling (symmetric % change)
network_iii = 1;
network_jjj = 1;

% previous
% network_iii = 7;
% network_jjj = 12;
% 
% network_iii = 1;
% network_jjj = 1;

network_v1 = v1_avg_ALL(network_iii,network_jjj,:);
network_v2 = v2_avg_ALL(network_iii,network_jjj,:);
network_v = [network_v1(:);network_v2(:)];
network_angv = reshape([angv1_avg_ALL(network_iii,network_jjj,:); angv2_avg_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_perlin = reshape([per_lin1_total_ALL(network_iii,network_jjj,:); per_lin2_total_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_perang = reshape([per_ang1_total_ALL(network_iii,network_jjj,:); per_ang2_total_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_P = reshape([P1_ALL(network_iii,network_jjj,:); P2_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_prolen = reshape([pro_len1_ALL(network_iii,network_jjj,:); pro_len2_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_promaxlen = reshape([pro_maxlen1_ALL(network_iii,network_jjj,:); pro_maxlen2_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_prolifetime = reshape([pro_lifetime1_ALL(network_iii,network_jjj,:); pro_lifetime2_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);
network_procount = reshape([pro_count1_ALL(network_iii,network_jjj,:); pro_count2_ALL(network_iii,network_jjj,:)],[Ndata*2 1]);

xplot = repmat(1:2,size(network_v,1),1);
%% Linear speed
f_linear = figure('Position',[50 50 200 300]);
hold on
ylabel('Linear Speed (\mum/h)')

swarmchart(xplot(:,1),spheroid_v,markersz,'filled')
swarmchart(xplot(:,2),network_v,markersz,'filled')
box = boxplot([spheroid_v, network_v],'symbol','','Colors','k','Labels',xname);
saveas(f_linear,'box_linearV.png')
%% Angular speed
f_ang = figure('Position',[50 50 200 300]);
hold on
ylabel('Angular Speed (rad/h)')

swarmchart(xplot(:,1),spheroid_angv,markersz,'filled')
swarmchart(xplot(:,2),network_angv,markersz,'filled')
box = boxplot([spheroid_angv, network_angv],'symbol','','Colors','k','Labels',xname);
saveas(f_ang,'box_angV.png')

%% Linear Persistence
f_perlin = figure('Position',[50 50 200 300]);
hold on
ylabel('Cell Persistence Ratio (Linear)')
ylim([0 1])

swarmchart(xplot(:,1),spheroid_perlin,markersz,'filled')
swarmchart(xplot(:,2),network_perlin,markersz,'filled')
box = boxplot([spheroid_perlin, network_perlin],'symbol','','Colors','k','Labels',xname);
saveas(f_perlin,'box_perlin.png')

%% Angular Persistence
f_perang = figure('Position',[50 50 200 300]);
hold on
ylabel('Cell Persistence Ratio (Anuglar)')
ylim([0 1])

swarmchart(xplot(:,1),spheroid_perang,markersz,'filled')
swarmchart(xplot(:,2),network_perang,markersz,'filled')
box = boxplot([spheroid_perang, network_perang],'symbol','','Colors','k','Labels',xname);
saveas(f_perang,'box_perang.png')

%% Polarization
f_P = figure('Position',[50 50 200 300]);
hold on
ylabel('Polarization vector')
ylim([0 1])

swarmchart(xplot(:,1),spheroid_P,markersz,'filled')
swarmchart(xplot(:,2),network_P,markersz,'filled')
box = boxplot([spheroid_P, network_P],'symbol','','Colors','k','Labels',xname);
saveas(f_P,'box_P.png')

%% Protrusion length
f_prolen = figure('Position',[50 50 200 300]);
hold on
ylabel('Protrusion Length (A.U.)')
ylim([0 inf])

swarmchart(xplot(:,1),spheroid_prolen,markersz,'filled')
swarmchart(xplot(:,2),network_prolen,markersz,'filled')
box = boxplot([spheroid_prolen, network_prolen],'symbol','','Colors','k','Labels',xname);
saveas(f_prolen,'box_prolen.png')

%% Protrusion max length
f_promaxlen = figure('Position',[50 50 200 300]);
hold on
ylabel('Max Protrusion Length (A.U.)')
ylim([0 inf])

swarmchart(xplot(:,1),spheroid_promaxlen,markersz,'filled')
swarmchart(xplot(:,2),network_promaxlen,markersz,'filled')
box = boxplot([spheroid_promaxlen, network_promaxlen],'symbol','','Colors','k','Labels',xname);
saveas(f_promaxlen,'box_promaxlen.png')

%% Protrusion lifetime
f_prolifetime = figure('Position',[50 50 200 300]);
hold on
ylabel('Protrusion Lifetime (min)')
ylim([0 inf])

swarmchart(xplot(:,1),spheroid_prolifetime*dt/60,markersz,'filled')
swarmchart(xplot(:,2),network_prolifetime*dt/60,markersz,'filled')
box = boxplot([spheroid_prolifetime*dt/60, network_prolifetime*dt/60],'symbol','','Colors','k','Labels',xname);
saveas(f_prolifetime,'box_prolifetime.png')

%% Protrusion rate
f_procount = figure('Position',[50 50 200 300]);
hold on
ylabel('Protrusion Rate (protrusions/hr)')
ylim([0 inf])

swarmchart(xplot(:,1),spheroid_procount/10,markersz,'filled')
swarmchart(xplot(:,2),network_procount/10,markersz,'filled')
box = boxplot([spheroid_procount/10, network_procount/10],'symbol','','Colors','k','Labels',xname);
saveas(f_procount,'box_procount.png')
%% PART C
%% Speed vs rate for polarization
markersz = 30;
scale = 10^5;
gamma64_idx = 4;
gamma2_idx = 9;
%%
xplot = par1_list_FIXED/scale; 
f_speed = figure('Position',[50 50 200 300]);
hold on
ylabel('Speed  (\mum/h)')
ylim([0 inf])

scatter(xplot, v1_avg_ALL(:,gamma2_idx), markersz, 'r')
scatter(xplot, v2_avg_ALL(:,gamma2_idx), markersz, 'r')

scatter(xplot, v1_avg_ALL(:,gamma64_idx), markersz, 'b')
scatter(xplot, v2_avg_ALL(:,gamma64_idx), markersz, 'b')

set(gca,'Xscale','log')
xlabel('rate k+')

legend('\gamma = 2', '', '\gamma =64', 'Position',[0.4 0.85 0.01 0.01])
saveas(f_speed,'speed_vs_rate.png')
%%
% %% Rename
% % Get all text files in the current folder
% files = dir('*_data.mat');
% rename_arr = [6, 7, 8, 9, 10];
% %%
% % Loop through each file 
% for id = 1:length(files)
%     % Get the file name 
%     [~, f,ext] = fileparts(files(id).name);
%     
%     ridx = strfind(f,'_data') - 1;
%     if str2double(f(ridx)) > 5
%         continue
%     end
%     newridx = num2str(rename_arr(str2double(f(ridx))));
%     f(ridx) = [];
%     f = insertBefore(f,'_data',newridx);
%     rename = strcat(f,ext); 
%     movefile(files(id).name, rename); 
% end