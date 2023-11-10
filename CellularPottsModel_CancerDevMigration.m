%% Cellular Potts Model for Cancer Development and Migration in Dense ECM
clear 

repeatN = 3; % number of simultation to repeat
run_counter = 0;
rng(11693); % set random seed for reproducibility
%%
for repeatN_idx = 1:repeatN

clearvars -except par1* par2* repeatN* run_counter

% rng(93); % set random seed for reproducibility
parsave = 0;
draweveryN = 100; % draw figure every N frame, 0 = do not draw
savemovie = 1; % 1 for saving avi
everyN = 100; % save movie every N frame
framerate = 10; % how many frame per second (fps)
subplotN = 1; % number of plot to make
maketiff = 0; % 1 for saving tiff from avi
% computeforce = 0; % 1 for computing force
Pscale = 10; % scale of quiver plot length for P


dx = 1; % um
dt = 20; % sec
tN = round((168*3600)/dt);
% tN = round((90*3600)/dt); %@%@% 168
MeshN = 128/dx;
stop_flipping = tN; % stop doing MCS after the Nth frame
exp_D_steadystate = 0; % on or off
[x, y] = meshgrid(1:MeshN, 1:MeshN);
scale = 10^5;

% q.J = 0; % cAMP degradation %%%
% q.alpha = 0; % cAMP basal generation
% q.thres = 0.3954/7.5; % spike

q.rho = 50/scale; % spike %%% used to be 1 [1.1, 2023/05/11]

% q.D = 0; % spike
% DifS = 0; %%% used to be 8
% DifD = 1; %%% um^2/s = 10^-8 cm^2/s
DifD = 0.0005/100; %%% um^2/s = 0.5 * 10^-9 cm^2/s


% q.A = 5*0.0178; % log sensing
% q.K = 1e-5; % log sensing
% q.Kg = 0.85*17.68; % G
% q.Kr = 4; % G
% q.a = 0.1863; % middle branch G
% q.tau = 0.3; % R
% q.Kv = 10000/scale; % P, production
% q.Kv = par1*par2/scale; % P, production [640, 2023/05/11], [2016, 2023/06/02]
q.Kv = 256/scale; % P, production [640, 2023/05/11], [2016, 2023/06/02]
% q.Kp = par1/scale; % P, degradation [10, 2023/05/11], [15.8, 2023/06/02]
q.Kp = q.Kv*0.0625; % P, degradation [10, 2023/05/11], [15.8, 2023/06/02]
q.KD1 = 0/scale; % D, local production [1.2, 2023/05/11]
% q.KD1 = (50-par1)*10/scale; % D, local production [1.2, 2023/05/11]
q.KD2 = 9.6/scale; % D, degradation // half-life 2hour ~= 1e-4 s^-1 [9.6]
q.Kdeg = 4444/scale; % D on S, s^-1 // literature ~= 16h^-1 [444]
q.Sthres = 0.1; % threshold of S to be treated as 0
% sigma = 0.6;

% par1_list = [0.1, 1, 10, 100];
% par2_list = [0.1, 0.5, 1, 2, 10];
% par1name = 'rate';
% par2name = 'fold';

p.ncell = 1; % number of cell, maximum 127 for int8
p.area = round(10*10/dx^2); % cell area
% p.areamin = p.area*0.2; % minimum initial cell area
p.area_divide = p.area*0.85; % threshold area that cell start dividing
p.growth = p.area*0.35/(18*3600/dt)*2; % growth in cell size per time step
p.growthmulti = 0.5;
p.xlam  = 300;  % Strength of Area constraint
p.perim = p.area*1.2; %sqrt(4*pi()*p.area)*1.27; % perimeter of cell
% p.perimstr = 0; % strength of perimeter constraint
% p.len = 5; % cell length along long axis/P
% p.lenstr = 10; % strength of length constraint
p.con   = -0;  % Strength of polarization
p.mmt   = -400;  % Strength of self momemtum
p.ecm   = 25000; % Strength of ECM
% p.spring = 8; % Strength of spring, -ve means it wants to extend
p.temp  =  10;  % Temperature
p.rseed = rand(); % Random seed, only used initially
% p.dx = 2; % micrometers (not actually used in the MCS functions)
%           % rather, it is used in getting velocity calculation
%           % and in the calculation of the effective grid size for diffusion
p.PottsIterations = (1/4);
% p.stng = 0;   % force strength tp signal
% p.rThresh = 0.2; % force threshold
p.maxcells = 86;
p.divide_t = (18*3600/dt)/2;
%%
scale1 = 5;
jbg = 300;
p.jma(1,1) = 0; % medium-medium
p.jma(1,2) = 400; % cell-medium
p.jma(2,1) = 400; % cell-medium
p.jma(2,2) = 300; % cell-cell

scale2 = 1;
p.jma(1,3) = 0; % medium-ECM
p.jma(3,1) = 0; % medium-ECM
p.jma(2,3) = 300; % cell-ECM
p.jma(3,2) = 300; % cell-ECM
p.jma(3,3) = 0; % ECM-ECM

rho_jma_k = p.jma(2,3); % modifier for q.rho to depend on Jce, neutral effect if sets as same to jce

% DcAMP = 400; % cAMP diffusion microns^2/min
% q.J=10;      % cAMP degradation rate (units: 1/min)
% 
% q.a = 0.1863;      % middle branch of cubic (Vasiev model)
% q.tau = 2.0;       % Relaxation variable's timescale (Vasiev's model)
% q.kr = 10.0;       % Strength of repressor on activator (Vasiev)
% q.kg=0.85*17.68;   % Amplitude of cubic (Vasiev calls this excitability)
% q.thresh=0.3954;   % When excitability exceeds this value, secrete cAMP
% 
% q.alph0=0.1*q.J;      % Baseline cAMP production / cell area
% q.D=1000*0.1*q.J;     % Amount of cAMP released / excited cell area (Not diffusion!!)
% q.A=5*0.0178;         % a in the logrithmic module I(S)
% q.K=1e-5;             % K in the logrithmic module I(S)
% q.dx=p.dx/sqrt(DcAMP);% effective grid size encodes the diffusion
% q.dt=0.002;           % time step (min), should be <(q.dx)^2/2
%%
% Set up ECM %
Smax = 1;
fiber_len = 5;
density = 0.75;
[~,dS] = deal(zeros([MeshN MeshN],'single')); % Activator
S = ones([MeshN MeshN],'single')*Smax ; % as ECM, homogeneous
% S = normrnd(Smax,sqrt(Smax/10),[MeshN MeshN]); % as ECM, guassian distribution

% S = makeECM([MeshN MeshN], fiber_len, density)*Smax; % as ECM, heterogeneous fiber
% pores_map = bwconncomp(1 - S, 4);
% pores_size_all = cellfun(@length, pores_map.PixelIdxList);
% fprintf('ECM density (fraction): %0.3f, Fiber length: %0.0f pixel,  Pore size: %0.1f +/- %0.1f pixel^2 \n', sum(S == 1, 'all')/MeshN^2, fiber_len, mean(pores_size_all), std(pores_size_all))
% figure('Position', [0, 0, 800, 350])
% subplot(1,2,1)
% imagesc(S)
% subplot(1,2,2)
% imagesc(labelmatrix(pores_map))
%%

% Additional ECM structure
% S = makeStructure(S,'hole',20,x,y);
% S = makeStructure(S,'stripe',5,x,y);

[D,dD] = deal(zeros([MeshN MeshN],'single')); % Degradation Enzyme

% [g,r,s,dg,dr] = deal(zeros([1 p.ncell],'single')); % concentration at ith cell
[nsig,ntyp,nsig_2cell] = deal(zeros([MeshN MeshN],'int8')); % cell label, cell type
[nar,nl,cmx,cmy,Px,Py] = deal(zeros([1 p.ncell],'single')); % cell area, center of mass x y, polarization
%%
% Pott Initialization %
ecc = 1;
cellside = round(sqrt(p.ncell));
initarea = 0.375*p.area;
lx = round(initarea^(1/2))*ecc; ly = round(initarea/lx);
rr = sqrt(initarea/pi());
space = 0;
label = datasample(1:p.ncell,p.ncell,'Replace',false); % assign random label nsig

% Make Movie %
if savemovie
    v = VideoWriter(['run' num2str(repeatN_idx) '_mov'],'MPEG-4');
    v.Quality = 99;
    v.FrameRate = framerate;
    open(v);
end

% % Square like %
% for i = 1:p.ncell
%     x0 = round(mod(i-1,cellside)*(lx+space) + MeshN/2 - lx*cellside/2);
%     y0 = round((floor((i-1)/cellside))*(ly+space) + MeshN/2 - ly*cellside/2);
%     nsig(y0+(1:ly), x0+(1:lx)) = label(i);
%     ntyp(y0+(1:ly), x0+(1:lx)) = 1;
% end

% Rounded cell %
for i = 1:p.ncell
    x0 = round(mod(i-1,cellside)*(rr+space) + MeshN/2);
    y0 = round((floor((i-1)/cellside))*(rr+space) + MeshN/2);
    
    nsig(sqrt((x - x0).^2 + (y - y0).^2) < rr) = label(i);
    ntyp(sqrt((x - x0).^2 + (y - y0).^2) < rr) = 1;
end

% % Strip like %
% for i = 1:p.ncell
%     x0 = round((i-1)*(lx+space) + MeshN/2 - lx*p.ncell/2);
%     y0 = round(MeshN/2);
%     nsig(y0+(1:ly), x0+(1:lx)) = label(i);
%     ntyp(y0+(1:ly), x0+(1:lx)) = 1;
% end


for i = 1:p.ncell
    nar(i) = sum((nsig == i),[1 2]);
    [row, col] = find(nsig == i);
        for j = 1:length(row)
            tmp = sum([nsig(row(j)+1,col(j)),nsig(row(j)-1,col(j)), nsig(row(j),col(j)+1), nsig(row(j),col(j)-1)] ~= nsig(row(j),col(j))); % count edge that is not the same cell
            nl(i) = nl(i) + tmp;
        end
    cmx(i) = sum(x .* (nsig == i), [1 2])/nar(i);
    cmy(i) = sum(y .* (nsig == i), [1 2])/nar(i);
end


% remove ECM where cells are initially plated at
S(ntyp == 1) = 0;

% assign ECM type
ntyp(S > 0) = 2;

% [x, y] = meshgrid(1:MeshN, 1:MeshN);
% for i = 1:p.ncell
%     tmp = ceil(rand(1,2)*(MeshN/2-sqrt(p.area))) + MeshN/4;
%     lx = 10;
%     ly = p.area/lx;
%     nsig(tmp(2):(tmp(2)+ly-1),tmp(1):(tmp(1)+lx-1)) = i;
%     ntyp(tmp(2):(tmp(2)+ly-1),tmp(1):(tmp(1)+lx-1)) = 1;
% end
% 
% for i = 1:p.ncell
%     nar(i) = sum((nsig == i),[1 2]);
%     cmx(i) = sum(x .* (nsig == i), [1 2])/nar(i);
%     cmy(i) = sum(y .* (nsig == i), [1 2])/nar(i);
% end
% 
% toosmall = find(nar < p.areamin);
% while ~isempty(toosmall)
%     for i = toosmall
%         % first remove the cells
%         if nar(i) ~= 0
%             tmp = (nsig == i);
%             nsig(tmp) = 0; 
%             ntyp(tmp) = 0;
%         end
%         
%         % then add them back
%         tmp = ceil(rand(1,2)*(MeshN/2-sqrt(p.area))) + MeshN/4;
%         lx = 10;
%         ly = p.area/lx;
%         nsig(tmp(2):(tmp(2)+ly-1),tmp(1):(tmp(1)+lx-1)) = i;
%         ntyp(tmp(2):(tmp(2)+ly-1),tmp(1):(tmp(1)+lx-1)) = 1;
%     end
%     
%     for i = 1:p.ncell
%         nar(i) = sum((nsig == i),[1 2]);
%         cmx(i) = sum(x .* (nsig == i), [1 2])/nar(i);
%         cmy(i) = sum(y .* (nsig == i), [1 2])/nar(i);
%     end
%     
%     toosmall = find(nar < 10);
% end

% Velocity Noise from CM position
Amp = 1;
cmxOld = cmx - Amp*(rand([1 p.ncell])-0.5); cmyOld = cmy - Amp*(rand([1 p.ncell])-0.5);
% [cmx_all, cmy_all] = deal(zeros(tN/100,p.ncell));

% Eight neighbours update
nbxupd = [1,1,0,-1,-1,-1,0,1];
nbyupd = [0,1,1,1,0,-1,-1,-1];

% Four neighbous update
% nbxupd = [1,1,0,0,-1,-1,0,0];
% nbyupd = [0,0,1,1,0,0,-1,-1];

% Random Noise
Amp = 0;
% S = 2*Amp*(rand([MeshN MeshN]));
% g = 0.2*Amp*(rand([1 p.ncell]));
% r = 0.2*Amp*(rand([1 p.ncell]));

% CM for whole pancake
CM = zeros(2,tN);
CM(1,1) = sum(x.*(nsig>0),[1 2])/sum(nar);
CM(1,2) = sum(y.*(nsig>0),[1 2])/sum(nar);

F = cell(2,tN);
[Vx, Vy, CMx_all, CMy_all, Px_all, Py_all] = deal(zeros(tN,p.maxcells,'single'));

% Target area list
targetarea = nar; % initialize target area

if exp_D_steadystate
    exp_D_Sy0t = zeros(tN,length(x0:MeshN));
    exp_D_Sr = zeros(tN,1);
end
%%
FIG = figure('Position',[10 10 subplotN*512 512]);
colormap turbo
ic = 0;

for i = 2:tN
    
%     if (p.ncell >= 25) && (mod(i,50)==0)
%         disp([nar(25), targetarea(25)]) %#%#%
%     end

    % slower growth and division %
    p.growth = p.growth*(p.growthmulti)^(1/tN);
    p.divide_t = p.divide_t*(1/p.growthmulti)^(1/tN);
    
    velx=(cmx-cmxOld)/dt; vely=(cmy-cmyOld)/dt; vel=(velx.^2 + vely.^2).^(1/2);
    
    % Divide %
%     if p.divide_t && ~mod(i,p.divide_t) && (p.ncell*2<p.maxcells)    % fixed time, all cells
%         % variable to double: 
%         % cmx, cmy, cmxOld, cmyOld, s, r, g, velx, vely, vel
%         outcell = divide_vars({cmx, cmy, nar, nl, Px, Py, s, r, g, velx, vely, vel});
%         [cmx, cmy, nar, nl, Px, Py, s, r, g, velx, vely, vel] = outcell{:};
%         
%         for k = 1:p.ncell
%             dauid = k + p.ncell; % daughter id
%             
%             ang = regionprops(nsig == k,'Orientation').Orientation;
%             ang_plus = ang + 90; % +ve minor axis
%             ang_minus = ang - 90; % -ve minor axis
%             
%             [row, col] = find(nsig == k);
%             cellang = atan2d((y - cmy(k)),(x - cmx(k)));
%             cell1 = (cellang <= ang_plus) & (cellang > ang_minus) & (nsig == k);
%             cell2 = ~cell1 & (nsig == k);
%             nsig(cell1) = k;
%             nsig(cell2) = dauid;
%             
%             nar(k) = sum((nsig == k),[1 2]);
%             nar(dauid) = sum((nsig == dauid),[1 2]);
%             nl(k)= regionprops(nsig == k,'Perimeter').Perimeter;
%             nl(dauid)= regionprops(nsig == dauid,'Perimeter').Perimeter;
%             cmx(k) = sum(x .* (nsig == k), [1 2])/nar(k);
%             cmy(k) = sum(y .* (nsig == k), [1 2])/nar(k);
%             cmx(dauid) = sum(x .* (nsig == dauid), [1 2])/nar(dauid);
%             cmy(dauid) = sum(y .* (nsig == dauid), [1 2])/nar(dauid);
% 
%         end

%     if p.divide_t && (rand() < 1/p.divide_t*p.ncell) && (p.ncell<p.maxcells)    % probabilistic, one cell
    if p.divide_t && (p.ncell<p.maxcells)
        validcells = find(nar >= p.area_divide);
        for j = 1:length(validcells)
            cellid = validcells(j);
            if (rand() < 1/p.divide_t)
                if p.ncell == 2
                    nsig_2cell = nsig;  % store end of 2 cell stage for circularity calculation
                    cmx_2cell = cmx;
                    cmy_2cell = cmy;
                    Px_2cell = Px;
                    Py_2cell = Py;
                    S_2cell = Snew;
                    i_2cell = i;
                end
                dauid = p.ncell + 1; % daughter id

                ang = regionprops(nsig == cellid,'Orientation').Orientation; % angle between x-axis and major axis
                ang_plus = ang + 90; % +ve minor axis
                ang_minus = ang - 90; % -ve minor axis

        %         [row, col] = find(nsig == cellid);
                cellang = -1*(atan2d((y - cmy(cellid)),(x - cmx(cellid))));
                cell1 = (cellang <= ang_plus) & (cellang > ang_minus) & (nsig == cellid);
                cell2 = ~cell1 & (nsig == cellid);
                
%                 if abs(sum(cell1, [1 2]) - sum(cell2, [1 2])) > nar(cellid)/4 % abort if division is too asymetric
%                     continue
%                 end
                
                nsig(cell1) = cellid;
                nsig(cell2) = dauid;
                
                outcell = divide_vars_single({cmx, cmy, nar, nl, Px, Py, velx, vely, vel}, cellid);
                [cmx, cmy, nar, nl, Px, Py, velx, vely, vel] = outcell{:};

                nar(cellid) = sum((nsig == cellid),[1 2]);
                nar(dauid) = sum((nsig == dauid),[1 2]);
                nl(cellid)= regionprops(nsig == cellid,'Perimeter').Perimeter;
                nl(dauid)= regionprops(nsig == dauid,'Perimeter').Perimeter;
                cmx(cellid) = sum(x .* (nsig == cellid), [1 2])/nar(cellid);
                cmy(cellid) = sum(y .* (nsig == cellid), [1 2])/nar(cellid);
                cmx(dauid) = sum(x .* (nsig == dauid), [1 2])/nar(dauid);
                cmy(dauid) = sum(y .* (nsig == dauid), [1 2])/nar(dauid);

                targetarea(cellid) = nar(cellid);
                targetarea(dauid) = nar(dauid);

                p.ncell = p.ncell + 1;
            end
        end
    end
    
%     if (mean(nar) >= p.area_divide) && p.divide_t && (rand() < 1/p.divide_t*p.ncell) && (p.ncell<p.maxcells)
% %     if (mean(nar) >= p.area_divide) && p.divide_t && ~mod(i,p.divide_t) && (p.ncell<p.maxcells)    % fixed time, one cell
%         cellid = randi([1,p.ncell]);
%         dauid = p.ncell + 1; % daughter id
%         
%         outcell = divide_vars_single({cmx, cmy, nar, nl, Px, Py, s, r, g, velx, vely, vel}, cellid);
%         [cmx, cmy, nar, nl, Px, Py, s, r, g, velx, vely, vel] = outcell{:};
%         
%         ang = regionprops(nsig == cellid,'Orientation').Orientation; % angle between x-axis and major axis
%         ang_plus = ang + 90; % +ve minor axis
%         ang_minus = ang - 90; % -ve minor axis
% 
% %         [row, col] = find(nsig == cellid);
%         cellang = -1*(atan2d((y - cmy(cellid)),(x - cmx(cellid))));
%         cell1 = (cellang <= ang_plus) & (cellang > ang_minus) & (nsig == cellid);
%         cell2 = ~cell1 & (nsig == cellid);
%         nsig(cell1) = cellid;
%         nsig(cell2) = dauid;
% 
%         nar(cellid) = sum((nsig == cellid),[1 2]);
%         nar(dauid) = sum((nsig == dauid),[1 2]);
%         nl(cellid)= regionprops(nsig == cellid,'Perimeter').Perimeter;
%         nl(dauid)= regionprops(nsig == dauid,'Perimeter').Perimeter;
%         cmx(cellid) = sum(x .* (nsig == cellid), [1 2])/nar(cellid);
%         cmy(cellid) = sum(y .* (nsig == cellid), [1 2])/nar(cellid);
%         cmx(dauid) = sum(x .* (nsig == dauid), [1 2])/nar(dauid);
%         cmy(dauid) = sum(y .* (nsig == dauid), [1 2])/nar(dauid);
%         
%         targetarea = nar;
%         
%         p.ncell = p.ncell + 1;
%     end
    
    targetarea = targetarea + p.growth;
    targetarea = min(targetarea,p.area);
    
    cmxOld=cmx; cmyOld=cmy;
%     Sx=circshift(S,[0,-1])-S;
%     Sy=circshift(S,[-1,0])-S;

%     % Singal S %
%     for xx = 2:(MeshN-1)
%         for yy = 2:(MeshN-1)
%             dS(xx,yy) = DifS*(1/dx^2)*(S(xx+1,yy) + S(xx-1,yy) + S(xx,yy+1) + S(xx,yy-1) - 4*S(xx,yy)); % diffusion       
%         end
%     end
%     
%     % PBC for S %
%     for yy = 2:(MeshN-1)
%        dS(1,yy) = DifS*(1/dx^2)*(S(2,yy) + S(MeshN,yy) + S(1,yy+1) + S(1,yy-1) - 4*S(1,yy));
%        dS(MeshN,yy) = DifS*(1/dx^2)*(S(1,yy) + S(MeshN-1,yy) + S(MeshN,yy+1) + S(MeshN,yy-1) - 4*S(MeshN,yy));
%     end
%     
%     for xx = 2:(MeshN-1)
%        dS(xx,1) = DifS*(1/dx^2)*(S(xx+1,1) + S(xx-1,1) + S(xx,2) + S(xx,MeshN) - 4*S(xx,1));
%        dS(xx,MeshN) = DifS*(1/dx^2)*(S(xx+1,MeshN) + S(xx-1,MeshN) + S(xx,1) + S(xx,MeshN-1) - 4*S(xx,MeshN));
%     end
    
%     dS(1,1) = DifS*(1/dx^2)*(S(2,1) + S(MeshN,1) + S(1,2) + S(1,MeshN) - 4*S(1,1));
%     dS(MeshN,1) = DifS*(1/dx^2)*(S(1,1) + S(MeshN-1,1) + S(MeshN,2) + S(MeshN,MeshN) - 4*S(MeshN,1));
%     dS(1,MeshN) = DifS*(1/dx^2)*(S(2,MeshN) + S(MeshN,MeshN) + S(1,1) + S(1,MeshN-1) - 4*S(1,MeshN));
%     dS(MeshN,MeshN) = DifS*(1/dx^2)*(S(1,MeshN) + S(MeshN-1,MeshN) + S(MeshN,1) + S(MeshN,MeshN-1) - 4*S(MeshN,MeshN));

    
    % Degradation D %
%     for xx = 2:(MeshN-1)
%         for yy = 2:(MeshN-1)
%             dD(xx,yy) = DifD*(1/dx^2)*(D(xx+1,yy) + D(xx-1,yy) + D(xx,yy+1) + D(xx,yy-1) - 4*D(xx,yy)); % diffusion       
%         end
%     end
    
    dD = DifD*(1/dx^2)*(circshift(D,[1 0]) + circshift(D,[-1 0]) + circshift(D,[0 1]) + circshift(D,[0 -1]) - 4*D);
    
    % PBC for D %
    for yy = 2:(MeshN-1)
       dD(1,yy) = DifD*(1/dx^2)*(D(2,yy) + D(MeshN,yy) + D(1,yy+1) + D(1,yy-1) - 4*D(1,yy));
       dD(MeshN,yy) = DifD*(1/dx^2)*(D(1,yy) + D(MeshN-1,yy) + D(MeshN,yy+1) + D(MeshN,yy-1) - 4*D(MeshN,yy));
    end
    
    for xx = 2:(MeshN-1)
       dD(xx,1) = DifD*(1/dx^2)*(D(xx+1,1) + D(xx-1,1) + D(xx,2) + D(xx,MeshN) - 4*D(xx,1));
       dD(xx,MeshN) = DifD*(1/dx^2)*(D(xx+1,MeshN) + D(xx-1,MeshN) + D(xx,1) + D(xx,MeshN-1) - 4*D(xx,MeshN));
    end
    
    dD(1,1) = DifD*(1/dx^2)*(D(2,1) + D(MeshN,1) + D(1,2) + D(1,MeshN) - 4*D(1,1));
    dD(MeshN,1) = DifD*(1/dx^2)*(D(1,1) + D(MeshN-1,1) + D(MeshN,2) + D(MeshN,MeshN) - 4*D(MeshN,1));
    dD(1,MeshN) = DifD*(1/dx^2)*(D(2,MeshN) + D(MeshN,MeshN) + D(1,1) + D(1,MeshN-1) - 4*D(1,MeshN));
    dD(MeshN,MeshN) = DifD*(1/dx^2)*(D(1,MeshN) + D(MeshN-1,MeshN) + D(MeshN,1) + D(MeshN,MeshN-1) - 4*D(MeshN,MeshN));
   
    
%     % No flux at boundary for S %
%     
%     for yy = 2:(MeshN-1)
%        dS(1,yy) = DcAMP*(1/dx^2)*(S(2,yy) + S(1,yy) + S(1,yy+1) + S(1,yy-1) - 4*S(1,yy));
%        dS(MeshN,yy) = DcAMP*(1/dx^2)*(S(MeshN,yy) + S(MeshN-1,yy) + S(MeshN,yy+1) + S(MeshN,yy-1) - 4*S(MeshN,yy));
%     end
%     
%     for xx = 2:(MeshN-1)
%        dS(xx,1) = DcAMP*(1/dx^2)*(S(xx+1,1) + S(xx-1,1) + S(xx,2) + S(xx,1) - 4*S(xx,1));
%        dS(xx,MeshN) = DcAMP*(1/dx^2)*(S(xx+1,MeshN) + S(xx-1,MeshN) + S(xx,MeshN) + S(xx,MeshN-1) - 4*S(xx,MeshN));
%     end
%     
%     dS(1,1) = DcAMP*(1/dx^2)*(S(2,1) + S(1,1) + S(1,2) + S(1,1) - 4*S(1,1));
%     dS(MeshN,1) = DcAMP*(1/dx^2)*(S(MeshN,1) + S(MeshN-1,1) + S(MeshN,2) + S(MeshN,1) - 4*S(MeshN,1));
%     dS(1,MeshN) = DcAMP*(1/dx^2)*(S(2,MeshN) + S(1,MeshN) + S(1,MeshN) + S(1,MeshN-1) - 4*S(1,MeshN));
%     dS(MeshN,MeshN) = DcAMP*(1/dx^2)*(S(MeshN,MeshN) + S(MeshN-1,MeshN) + S(MeshN,MeshN) + S(MeshN,MeshN-1) - 4*S(MeshN,MeshN));
    
%     dS =  dS + q.alpha - q.J*S;
    
%     nsig = nsig + 1;
%     g = [0 g];
%     dS = dS + (g(nsig) > q.thres)*q.rho*q.D;
%     cmidx = round(cmy)+(round(cmx)-1)*MeshN;
%     dS = dS - (nsig > 1)*q.rho.*S - q.Kdeg.*D.*S; %%% ecm degradation at whole cell, old
%     dS = dS - q.Kdeg.*D.*S; %%% ecm degradation at whole cell, new
    dS = - q.Kdeg.*D.*S;
%     dS(nsig == 2) = dS(nsig == 2) - q.rho.*S(nsig == 2);
%     dS(cmidx) = dS(cmidx) - q.rho.*S(cmidx); %%% ecm degradation at center of mass
%     dS(cmidx(1)) = dS(cmidx(1)) - q.rho.*S(cmidx(1)); %%% ecm degradation at center of mass of first cell
%     g(1) = [];
%     nsig = nsig -1;
    
    if sum(isnan(dS(:)))
        pause
    end
    
    tipmap = getPolarizedTip(nsig, ntyp, Px, Py, cmx, cmy, x, y, nbxupd, nbyupd, MeshN);
    dD = dD + (nsig > 0)*q.rho*rho_jma_k/p.jma(2,3) + q.KD1*tipmap - q.KD2*D - q.Kdeg*D.*S;
    
%     if i < 5000
%         if(mod(i,100) < 5)
%             for j=1:1
%                 x=round((199*rand+1));
%                 y=round((199*rand+1));
%                 S(x,y)=S(x,y)+50000*rand;
%             end
%         end
%     end
%     if i < 5000
%     if (mod(i,100) == 0)
%         for j=1:2
%             x=round((100*rand+50));
%             y=round((100*rand+50));
%             S(x,y)=S(x,y)+15000;
%         end
%     end
%     end
    
    Snew = S + dS*dt; % + normrnd(0,q,MeshN,MeshN)*sqrt(dt);
    Dnew = D + dD*dt;
    
%     % Total s in cell %
%     for j = 1:p.ncell
%         s(j) = sum(S(nsig == j));
%     end
%     % Input I, logarithm sensing %
%     I = q.A*log10(1 + s/q.K./nar);
%     % Activator g in cell %
%     dg = -q.Kg*g.*(g-1).*(g-q.a) - q.Kr*r + I;
%     gnew = g + dg*dt;
%     % Repressor r in cell %
%     dr = (g-r)/q.tau;
%     rnew = r + dr*dt;
    
    S = Snew; 
%     g = gnew; r = rnew;
    D = Dnew;
    
    % Polarization P of cell %
    
    % Simple Euler %
%     P = (Px.^2 + Py.^2).^(1/2);
%     dPx = -q.Kp*Px + q.Kv*velx./(vel+1e-10).*(1-P);
%     dPy = -q.Kp*Py + q.Kv*vely./(vel+1e-10).*(1-P);
%     Px = Px + dPx*dt;
%     Py = Py + dPy*dt;
    
    for ttt = 1:100
        P = (Px.^2 + Py.^2).^(1/2);
        dPx = -q.Kp*Px + q.Kv*velx./(vel+1e-10).*(1-P);
        dPy = -q.Kp*Py + q.Kv*vely./(vel+1e-10).*(1-P);
        Px = Px + dPx*dt/100;
        Py = Py + dPy*dt/100;
    end
    
    % RK4 %
%     dP = @(P1,P2,v1,v,q) -q.Kp*P1 + q.Kv*v1./(v+1e-10).*(1-(P1.^2 + P2.^2).^(1/2));
%     
%     Px_k1 = dP(Px,Py,velx,vel,q);
%     Py_k1 = dP(Py,Px,vely,vel,q);
%     
%     Px_k2 = dP(Px + Px_k1*dt/2,Py + Py_k1*dt/2,velx,vel,q);
%     Py_k2 = dP(Py + Py_k1*dt/2,Px + Px_k1*dt/2,vely,vel,q);
%     
%     Px_k3 = dP(Px + Px_k2*dt/2,Py + Py_k2*dt/2,velx,vel,q);
%     Py_k3 = dP(Py + Py_k2*dt/2,Px + Px_k2*dt/2,vely,vel,q);
% 
%     Px_k4 = dP(Px + Px_k3*dt,Py + Py_k3*dt,velx,vel,q);
%     Py_k4 = dP(Py + Py_k3*dt,Px + Px_k3*dt,vely,vel,q);
% 
%     Px = Px + dt*(Px_k1 + 2*Px_k2 + 2*Px_k3 + Px_k4)/6;
%     Py = Py + dt*(Py_k1 + 2*Py_k2 + 2*Py_k3 + Py_k4)/6;
    

    
    % Initial Condition %
%     % half plane wave at middle --> spiral
%     if abs(i-900) < 300
%         S((-MeshN/4:MeshN/4)+MeshN/4+1,(-1)+MeshN/2) = 40;
%         R((-MeshN/4:MeshN/4)+MeshN/4+1,(0:1)+MeshN/2) = 10;
%     end
    
%     % plane wave from right --> plane wave
%     if abs(i-900) < 300
%         S(1:MeshN, end-5) = 40;
%         S(1:MeshN, end-4:end) = 5;
%     end
    %%
    % Pott
    % Find cell overlapping boundaries
%     iover = zeros(p.ncell,p.ncell,'single');
%     for xx = 1:MeshN-1
%         for yy = 1:MeshN-1
%             if (nsig(yy+1,xx) ~= nsig(yy,xx)) && (nsig(yy+1,xx)*nsig(yy,xx) ~= 0)
%                 iover(nsig(yy,xx), nsig(yy+1,xx)) = iover(nsig(yy,xx), nsig(yy+1,xx)) + 1;
%                 iover(nsig(yy+1,xx), nsig(yy,xx)) = iover(nsig(yy+1,xx), nsig(yy,xx)) + 1;
%             end
%             if (nsig(yy,xx+1) ~= nsig(yy,xx)) && (nsig(yy,xx)*nsig(yy,xx+1) ~= 0)
%                 iover(nsig(yy,xx), nsig(yy,xx+1)) = iover(nsig(yy,xx), nsig(yy,xx+1)) + 1;
%                 iover(nsig(yy,xx+1), nsig(yy,xx)) = iover(nsig(yy,xx+1), nsig(yy,xx)) + 1;
%             end
%         end
%     end
    
    % Now calculate the force
    [conx, cony] = deal(zeros(1,p.ncell,'single'));
% 	for ii = 1:p.ncell
% 		forcex = 0;
% 		forcey = 0;
% 		for jj = 1:p.ncell
% 			forcex = forcex + iover(ii,jj)*velx(jj);
% 			forcey = forcey + iover(ii,jj)*vely(jj);
%         end
% 		force = forcex*forcex + forcey*forcey;
% 		if force ~= 0
% 			conx(ii) = forcex*p.con/sqrt(force);
% 			cony(ii) = forcey*p.con/sqrt(force);
%         end
%     end
    
%     r = interp2(1:MeshN,1:MeshN,R,cmx,cmy); % for the moment being
%     conx = conx - p.stng*interp2(1:MeshN,1:MeshN,Sx,cmx,cmy) .* (r<p.rThresh);
%     cony = cony - p.stng*interp2(1:MeshN,1:MeshN,Sy,cmx,cmy) .* (r<p.rThresh);
%%
    RAN = ceil(rand([2 round(MeshN^2*p.PottsIterations)])*(MeshN-4))+2; % wall
%     RAN = ceil(rand([2 MeshN^2*p.PottsIterations])*(MeshN)); % PBC
    RAN8 = ceil(rand([1 round(MeshN^2*p.PottsIterations)])*8);
    
    if i <= stop_flipping
        for ii = 1:round(MeshN^2*p.PottsIterations)
           % Randomly pick a site to flip       
           xran = RAN(1,ii);
           yran = RAN(2,ii);

    %        [YBD, XBD] = find(abs((circshift(nsig,[0,1])-nsig)) + abs((circshift(nsig,[0,-1])-nsig)) + abs((circshift(nsig,[1,0])-nsig)) + abs((circshift(nsig,[-1,0])-nsig)));
    %        ranBD = ceil(rand()*length(XBD));
    %        xran = XBD(ranBD);
    %        yran = YBD(ranBD);
    %        if xran < 3 || xran > 198 || yran < 3 || yran > 198
    %            continue
    %        end

           % Randomly pick a neighbour site
           nbidx = RAN8(ii);
           xnb = xran + nbxupd(nbidx);
    %        if xnb == (MeshN + 1)
    %            xnb = 1;
    %        else if xnb == 0
    %                xnb = MeshN;
    %            end
    %        end
           ynb = yran + nbyupd(nbidx);
    %        if ynb == (MeshN + 1)
    %            ynb = 1;
    %        else if ynb == 0
    %                ynb = MeshN;
    %            end
    %        end

    %        if (nsig(yran,xran) ~= nsig(ynb,xnb))
           if (nsig(yran,xran) ~= nsig(ynb,xnb)) && (ntyp(ynb,xnb)~=2)

               % Connectivity constraint
    %             ran_connect = nsig(((xran+nbxupd)-1)*MeshN+(yran+nbyupd)) == nsig(yran,xran);
    %             nb_connect = nsig(((xran+nbxupd)-1)*MeshN+(yran+nbyupd)) == nsig(ynb,xnb);
    %             if (sum((ran_connect - circshift(ran_connect,1)) == 1) ~= 1) && (nsig(yran,xran) ~=0)
    %                 continue
    %             end
    %             if (sum((nb_connect - circshift(nb_connect,1)) == 1) ~= 1) && (nsig(ynb,xnb) ~=0)
    %                 continue
    %             end

                dE = energies(p, xran, yran, xnb, ynb, MeshN, nsig, ntyp, nar, nl, cmx, cmy, velx, vely, Px, Py, conx, cony, nbxupd, nbyupd, S, x ,y, targetarea);

                prob = exp(-dE/p.temp);
                ran = rand();

                if (ran <= prob)
                    % update center of mass and area %
                    if (ntyp(yran,xran) == 1) 
                        term1 = nar(nsig(yran,xran))/(nar(nsig(yran,xran))-1)*cmx(nsig(yran,xran));
                        term2 = xran/(nar(nsig(yran,xran))-1);
                        cmx(nsig(yran,xran)) = term1 - term2;
                        term1 = nar(nsig(yran,xran))/(nar(nsig(yran,xran))-1)*cmy(nsig(yran,xran));
                        term2 = yran/(nar(nsig(yran,xran))-1);
                        cmy(nsig(yran,xran)) = term1 - term2;

                        nar(nsig(yran,xran)) = nar(nsig(yran,xran)) - 1;
    %                     nl(nsig(yran,xran)) = nl(nsig(yran,xran)) + (((nsig(yran+1,xran) ~= nsig(yran,xran)) + (nsig(yran-1,xran) ~= nsig(yran,xran)) + (nsig(yran,xran+1) ~= nsig(yran,xran)) + (nsig(yran,xran-1) ~= nsig(yran,xran))) - 2)*2;
                        nl(nsig(yran,xran)) = nl(nsig(yran,xran)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(yran,xran))-2)*(-2);
                    end

                    if (ntyp(ynb,xnb) == 1)
                        term1 = nar(nsig(ynb,xnb))/(nar(nsig(ynb,xnb))+1)*cmx(nsig(ynb,xnb));
                        term2 = xran/(nar(nsig(ynb,xnb))+1);
                        cmx(nsig(ynb,xnb)) = term1 + term2;
                        term1 = nar(nsig(ynb,xnb))/(nar(nsig(ynb,xnb))+1)*cmy(nsig(ynb,xnb));
                        term2 = yran/(nar(nsig(ynb,xnb))+1);
                        cmy(nsig(ynb,xnb)) = term1 + term2;

                        nar(nsig(ynb,xnb)) = nar(nsig(ynb,xnb)) + 1;
    %                     nl(nsig(ynb,xnb)) = nl(nsig(ynb,xnb)) - (((nsig(yran+1,xran) ~= nsig(ynb,xnb)) + (nsig(yran-1,xran) ~= nsig(ynb,xnb)) + (nsig(yran,xran+1) ~= nsig(ynb,xnb)) + (nsig(yran,xran-1) ~= nsig(ynb,xnb))) - 2)*2;
                        nl(nsig(ynb,xnb)) = nl(nsig(ynb,xnb)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(ynb,xnb))-2)*2;
                    end

                    % Flip it %
                    if ntyp(yran,xran) == 2
                        S(yran,xran) = 0;
                    end

                    nsig(yran,xran)=nsig(ynb,xnb);
                    ntyp(yran,xran)=ntyp(ynb,xnb);
    % 				ic = ic + 1;
                end
           end
        end
    end
    
    % Compute force from Hamiltonian %  
%     if computeforce
%         [Fx,Fy] = ComputeForce(p, MeshN, nsig, ntyp, nar, nl, cmx, cmy, velx, vely, Px, Py, conx, cony, nbxupd, nbyupd, S);
%         F{1,i} = Fx; F{2,i} = Fy;
%         CM(1,i) = sum(x.*(nsig>0),[1 2])/sum(nar);
%         CM(2,i) = sum(y.*(nsig>0),[1 2])/sum(nar);
%     end
    %%
    
    Vx(i,1:p.ncell) = velx*dt;
    Vy(i,1:p.ncell) = vely*dt;
    CMx_all(i,1:p.ncell) = cmx;
    CMy_all(i,1:p.ncell) = cmy;
    Px_all(i,1:p.ncell) = Px;
    Py_all(i,1:p.ncell) = Py;
    CM(1,i) = mean(cmx);
    CM(2,i) = mean(cmy);
    
    if exp_D_steadystate
        exp_D_Sy0t(i,:) = S(y0,x0:MeshN);
        exp_D_Sr(i) = find(S(y0,x0:MeshN) <= q.Sthres, 1, 'last');
    end
    
%     avg_cellpancake_R = sum(((cmx - CM(i,1)).^2 + (cmy - CM(i,2)).^2).^0.5)/p.ncell;
%     S_static(i,1) = S(round(CM(i,2)),round(CM(i,1) + avg_cellpancake_R));
%     s_cells(i,:) = s;
    
    %%
    if (mod(i,draweveryN) == 0) && (draweveryN ~= 0)
        
        if subplotN > 1
            subplot(1,2,1)
            imagesc(Snew,[0 Smax])
            imagesc(Dnew,[0 1])
            title(num2str(i))
        end
%         
%         if computeforce
%             subplot(1,2,1)
%     %         imagesc((Fx.^2+Fy.^2).^(1/2))
%             quiver(Fx,Fy,10)
%             title(num2str(i))
%             set(gca,'Ydir','reverse')
%         end
        
        if subplotN > 1
            subplot(1,subplotN,2)
            imagesc(nsig*floor(87/p.ncell)+int8((nsig>0)*30)+int8(Snew/Smax*10),[0 128])
            hold on 
%             quiver([1 0 cmx],[0 0 cmy],[1 0 Px],[0 0 Py],0.5,'color','white','LineWidth',1)
            quiver([cmx],[cmy],[Px]*Pscale,[Py]*Pscale,'off','color','white','LineWidth',1)
            hold off
        else
            imagesc(nsig*floor(87/p.ncell)+int8((nsig>0)*30)+int8(Snew/Smax*10),[0 128])
            hold on 
%             quiver([0 MeshN cmx],[0 MeshN cmy],[0 0 Px]*10,[-1 +1 Py]*10,'off','color','white','LineWidth',1)
            quiver([cmx],[cmy],[Px]*Pscale,[Py]*Pscale,'off','color','white','LineWidth',1)
            hold off
        end
%         imagesc(nsig)
        title([num2str(i) ' , ' num2str(round(i*dt/60/60,1)) ' hours'])
        pause(0.01)
        
% %         cmx_all(i/50,:) = cmx;
% %         cmy_all(i/50,:) = cmy;
%         saveas(FIG,['C:\MATLAB_data\' 'Sim' num2str(i,['%0' num2str(length(num2str(tN))) 'd']) '.png'])
%         print(['C:\MATLAB_data\' 'Sim' num2str(i,['%0' num2str(length(num2str(tN))) 'd'])],'-dpng')
        if savemovie && (mod(i,everyN) == 0)
            drawnow; %Force the figure to render
%             frame = getframe; % get the latest frame plotted
%             frame = uint8(nsig*(127/p.ncell));
            frame = getframe(FIG); %Convert the figure to a movie frame of figure
            writeVideo(v, frame); %Write the frame to the movie file
        end
    end
    
%     Vx(i,:) = velx*dt;
%     Vy(i,:) = vely*dt;
%     CMx_all(i,:) = cmx;
%     CMy_all(i,:) = cmy;
    S(S <= q.Sthres) = 0;
    ntyp((S <= q.Sthres) & (ntyp==2)) = 0;
    
end


%% END of main simulation --> Post-sim processes and calculations
if savemovie
	close(v);
end 

if maketiff
    if ~savemovie
        disp('Need to save avi movie for making tiff!')
    else
        obj = VideoReader('mov.avi');
        vid = read(obj);
        frames = obj.NumberOfFrames;
        for i = 1:frames
            if i == 1
                imwrite(vid(:,:,1,i),['simcells_n' num2str(p.ncell) '_C2.tif'])
            else
                imwrite(vid(:,:,1,i),['simcells_n' num2str(p.ncell) '_C2.tif'],'WriteMode','append')
            end
        end
    end
end

%% calculate circularity
if sum(nsig_2cell(:)) > 0
    mask = (nsig_2cell > 0);
else
    mask = (ntyp == 1);
end

% old method (less accurate, give >1 circularity)
% props = regionprops(mask, 'Area', 'Perimeter');
% allAreas = [props.Area];
% allPerimeters = [props.Perimeter];
% circularity = (4 * pi * allAreas)./(allPerimeters .^ 2);

% more accurate method using alphaShape
[Ip,Jp] = find(mask);
Shp = alphaShape(Ip, Jp, 2);
circularity = 4*pi*area(Shp)/perimeter(Shp)^2;
%% 2-cell stage velocity
twocell_t0 = find(CMx_all(:,2) > 0,1,'first');
skiphour = 5; % only start calculating the cell velocities after 'skiphour' hours since t0
totalhour = 10; % calculate velocity for 'totalhour' hours
twocell_ti = twocell_t0 + skiphour*3600/dt;
twocell_tf = twocell_ti + totalhour*3600/dt;

if twocell_tf > tN
    exitnote = 'two cell stage is less than 10 hours.';
    close(FIG)
    run_counter = run_counter + 1;
    continue
end

t_points = twocell_ti:(3600/dt):twocell_tf;

% Rectilinear Speed %
[v1,v2,vx1,vy1,vx2,vy2] = deal(zeros([length(t_points)-1 1],'double'));
for i = 1:length(t_points)-1
   vx1(i) = (CMx_all(t_points(i+1),1) - CMx_all(t_points(i),1))*dx/((t_points(i+1) - t_points(i))*dt/3600); % um/hr
   vy1(i) = (CMy_all(t_points(i+1),1) - CMy_all(t_points(i),1))*dx/((t_points(i+1) - t_points(i))*dt/3600);
   v1(i) = sqrt(vx1(i)^2 + vy1(i)^2);
   
   vx2(i) = (CMx_all(t_points(i+1),2) - CMx_all(t_points(i),2))*dx/((t_points(i+1) - t_points(i))*dt/3600); % um/hr
   vy2(i) = (CMy_all(t_points(i+1),2) - CMy_all(t_points(i),2))*dx/((t_points(i+1) - t_points(i))*dt/3600);
   v2(i) = sqrt(vx2(i)^2 + vy2(i)^2);
end
v1_avg = mean(v1);
v2_avg = mean(v2);

% Angular Speed %
[angv1, angv2] = deal(zeros([length(t_points)-1 1],'double'));
for i = 1:length(t_points)-1
    CMavg = (CM(:,t_points(i+1)) + CM(:,t_points(i)))/2;
    
    dCMy1f = CMy_all(t_points(i+1),1) - CMavg(2);
    dCMx1f = CMx_all(t_points(i+1),1) - CMavg(1);
    dCMy1i = CMy_all(t_points(i),1) - CMavg(2);
    dCMx1i = CMx_all(t_points(i),1) - CMavg(1);
    dang1 = angdiff(atan2d(dCMy1i,dCMx1i)*(pi()/180),atan2d(dCMy1f,dCMx1f)*(pi()/180));
    angv1(i) = dang1/((t_points(i+1) - t_points(i))*dt/3600); % rad/hr
    
    dCMy2f = CMy_all(t_points(i+1),2) - CMavg(2);
    dCMx2f = CMx_all(t_points(i+1),2) - CMavg(1);
    dCMy2i = CMy_all(t_points(i),2) - CMavg(2);
    dCMx2i = CMx_all(t_points(i),2) - CMavg(1);
    dang2 = angdiff(atan2d(dCMy2i,dCMx2i)*(pi()/180),atan2d(dCMy2f,dCMx2f)*(pi()/180));
%     dang2 = mod(atan2d(dCMy2f,dCMx2f) - atan2d(dCMy2i,dCMx2i),360);
    angv2(i) = dang2/((t_points(i+1) - t_points(i))*dt/3600); % rad/hr
%     angv2(i) = min(360-dang2,dang2)*(pi()/180)/((t_points(i+1) - t_points(i))*dt/3600); % rad/hr
end
angv1_avg = mean(abs(angv1));
angv2_avg = mean(abs(angv2));
%%
% Persistence %
[per_lin1,per_lin2,per_ang1,per_ang2, per_col_lin, per_col_ang] = deal(zeros([length(t_points)-2 1],'single'));

for i = 1:length(t_points)-2
    % cell %
    per_lin1(i) = (vx1(i+1)*vx1(i) + vy1(i+1)*vy1(i)) >= 0;
    per_lin2(i) = (vx2(i+1)*vx2(i) + vy2(i+1)*vy2(i)) >= 0;
    per_ang1(i) = abs(sign(angv1(i+1)) - sign(angv1(i)))/(-2) + 1;
    per_ang2(i) = abs(sign(angv2(i+1)) - sign(angv2(i)))/(-2) + 1;
    
    % collective %
    per_col_lin(i) = vx1(i)*vx2(i) + vy1(i)*vy2(i) >= 0;
    per_col_ang(i) = angv1(i)*angv2(i) >= 0;
end
per_lin1_total = sum(per_lin1)/length(per_lin1);
per_lin2_total = sum(per_lin2)/length(per_lin2);
per_ang1_total = sum(per_ang1)/length(per_ang1);
per_ang2_total = sum(per_ang2)/length(per_ang2);
per_col_lin_total = sum(per_col_lin)/length(per_col_lin);
per_col_ang_total = sum(per_col_ang)/length(per_col_ang);
%%
% Protrusion %
ang_thres = 10;
len_thres = 12;
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

valid1 = protrusion_len1 >= len_thres;
valid2 = protrusion_len2 >= len_thres;

protrusion_lifetime1(~valid1) = []; protrusion_len1(~valid1) = []; protrusion_P_avg1(~valid1) = [];
protrusion_lifetime2(~valid2) = []; protrusion_len2(~valid2) = []; protrusion_P_avg2(~valid2) = [];

protrusion_lifetime1_avg = mean(protrusion_lifetime1);
protrusion_len1_avg = mean(protrusion_len1);
protrusion_P_avg1_avg = mean(protrusion_P_avg1);
protrusion_lifetime2_avg = mean(protrusion_lifetime2);
protrusion_len2_avg = mean(protrusion_len2);
protrusion_P_avg2_avg = mean(protrusion_P_avg2);
%%
[P1, P2] = deal(zeros([length(t_points)-1 1],'double'));
% Px1 = Px_all(:,1);
% Py1 = Py_all(:,1);
% dPx1 = diff(Px_all(:,1));
% dPy1 = diff(Py_all(:,1));

% calculate steady state P given model parameter
% gamma = q.Kv/q.Kp;
% if gamma ~= 1
%     Pss = (1-sqrt(1/gamma))/(1-1/gamma);
% else
%     Pss = 0.5;
% end

for i = 1:(length(t_points)-1)
    P1(i) = mean(P_all(t_points(i):t_points(i+1),1),[1]);
    P2(i) = mean(P_all(t_points(i):t_points(i+1),2),[1]);
end
P1_avg = mean(P1);
P2_avg = mean(P2);
%% exp_D pl

if exp_D_steadystate
    plot_everyN = 100;
    figure
    plot((1:plot_everyN:tN)*dt/3600,exp_D_Sr(1:plot_everyN:tN))

    figure
    hold on
    for j = 1:plot_everyN:tN
        plot(x0:MeshN,exp_D_Sy0t(j,:))
    end
    hold off
end

close(FIG)

save(['run' num2str(repeatN_idx) '_simcells_n' num2str(p.ncell) '_data.mat'],'Vx','Vy','CMx_all','CMy_all','MeshN')

run_counter = run_counter + 1;
fprintf('repeat:%s [%d/%d]\n',num2str(repeatN_idx),run_counter,repeatN);
end

% %%
% % Calculate angular velocity %
% COMx_cluster = sum(CMx_all, 2)./sum(CMx_all~=0, 2);
% COMy_cluster = sum(CMy_all, 2)./sum(CMy_all~=0, 2);
% 
% t_2cells_start = find(sum(CMx_all~=0,[2])==2,1,'first');
% t_2cells_end = find(sum(CMx_all~=0,[2])==2,1,'last');
% 
% rho = ((CMx_all - COMx_cluster).^2 + (CMy_all - COMy_cluster).^2).^0.5;
% angv = ((CMx_all - COMy_cluster).*Vy - (CMy_all - COMy_cluster).*Vx)./rho.^2; % w = r cross v/r^2
% angv_avg_2cells = mean(angv(t_2cells_start:t_2cells_end,1:2),[2]); % rad/frame
% v_avg_2cells = mean((Vx(t_2cells_start:t_2cells_end,1:2).^2 +  Vy(t_2cells_start:t_2cells_end,1:2).^2).^0.5,[2]);
% %%
% figure
% for i = 1:tN
%     scatter(CMx_all(i,:),CMy_all(i,:))
%     xlim([0 MeshN])
%     ylim([0 MeshN])
%     pause(0.01)
% end

% %%
% figure
% hold on
% for i = 1:p.ncell
%     plot(1:size(cmx_all,1),cmx_all(:,i))
% end
% 
% figure
% hold on
% for i = 1:p.ncell
%     plot(1:size(cmy_all,1),cmy_all(:,i))
% end

% %%
% 
% test = [0,1,0,0,1,1,1,1];
% circshift(test,1);
% sum((test - circshift(test,1)) == 1)