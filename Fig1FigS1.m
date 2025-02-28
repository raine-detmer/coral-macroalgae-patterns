% README: code for making Figure 1 and Figure S1


%% set up

% colors for plots
Mcol = [0.4667 0.6745 0.1882]; % macroalgae
Ccol = [0.3020 0.7451 0.9333]; % coral
Hcol = [0.9294 0.6941 0.1255]; % herbivores

% make custom colormap for surface plots
vec = [100; 0]; % nodes at which to place the colors
raw = [Mcol; Ccol]; % gradient ranges from Mcol to Ccol
N = 256;
CMmap = interp1(vec, raw, linspace(100,0, N),'pchip');

% region of bistability
flow = 0.1674; % lower tipping point (calculated in Fig2.m)
fup = 0.1878; % upper tipping point

%% Briggs: ODE (non-spatial) bifurcation diagram

% make the non-spatial bifurcation diagram

% define the symbols (state variables)
syms Mi C H Mv

% define parameters
gTC = 0.1; % growth of coral over turf/free space
gamma = 0.4; % growth of macroalgae over coral vs. turf
gTI = 0.4; % growth of invulnerable macroalgae on turf
rM = 0.5; % production of vulnerable macroalgae by invulnerable macroalgae
gTV = 0.2; % growth of vulnerable macroalgae on turf
dv = 2; % vulnerable macroalgae mortality rate
omega = 2; % maturation rate of macroalgae from vulnerable to invulnerable stage
di = 0.4; % invulnerable macroalgae mortality rate
phiC = 0.01; % coral external recruitment rate
dC = 0.02; % coral mortality rate
phiM = 0.01; % macroalgae external recruitment rate

% herbivore parameters
rH = 0.2; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality rate
f = 0; % herbivore fishing pressure
phiH = 0.05; % external recruitment rate of herbivores

% set of fishing pressures at which to calculate equilibria 
fset = linspace(0.12, 0.2, 120); % 120 fishing pressures ranging from 0.12 to 0.2

% holding vectors for equilibrium values
Cstars = NaN(length(fset), 4); % coral
Mstars = NaN(length(fset), 4); % macroalgae (vuln + invuln)

% turn off warning about numerical stability
warning('off','symbolic:numeric:NumericalInstability')

% calculate the model equilibria for each fishing pressure in fset
for i = 1:length(fset)%for each element from 1 to length of fset
    
    fi = fset(i); % set fishing pressure to the ith element of fset

    % get the model eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the equilibrium values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just want positive and real equilibria 
    % store the values of the equilibrium C and M cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

% process results

% need to rearrange the equilibria to get smooth lines when plotting
% bistability region = 3 real, positive equilibria so 3rd column of Cstars
% is not NA in this region but is NaN outside of it
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region = last row where 3rd column of Cstars is not NaN (isnan = 0)
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region = first row where 3rd column of Cstars is not NaN

% use vertcat to concatenate vertical vectors
% for f on x axis: want to make sure the vectors containing each equilibria
% (stable upper, stable lower, and unstable) has the same length as fset
% but is NaN in the places where each of these equilibria don't exist
Cups = vertcat(Cstars(1:bstart-1, 1), Cstars(bstart:bend, 3), Cstars(bend+1:end, 4)); % upper (high C) stable equilibrium; need to make sure the length stays the same so concatenate with NaNs from other columns
Cmids = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:bend, 2), Cstars(bend+1:end, 4)); % middle (unstable) equilibrium
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 1)); % lower (low C) stable equilibrium

% repeat for macroalgae 
Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still (Mups =
% equilibrium macroalgal cover at the high-coral equilibria, Mlows =
% equilibrium macroalgal cover at the low-coral equilibria)


%% Briggs: PDE bifurcation diagram

% PDE parameters
diffs = [0.05,0.05,0.25, 0]; % diffusion rates of MI, C, H, and Mv
taxisM = 0; % herbivore taxis rate toward macroalgae
taxisC = -0.75; % herbivore taxis rate toward coral
taxisT = 0; % herbivore taxis rate toward turf/free space

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges

% space parameters
len = 400; % length of habitat
xset = linspace(-len/2,len/2,800); % spatial points at which to record model output

% time parameters
t_end = 3*50000; % number of timesteps
tset = linspace(0,t_end,2*2500); % timepoints at which to record model output


% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85; % coral cover at locations in initial coral patches
C0low = 0.05; % coral cover at locations in initial macroalgal patches
M0high = 0.85; % total macroalgal cover (MI + Mv) at locations in initial macroalgal patches
M0low = 0.05; % total macroalgal cover (MI + Mv) at locations in initial coral patches

% for icchoice = 3 (random initial conditions)
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4 (step-wise initial patches)
C0widths = round(length(xset)/16);  % step widths
initC = stepfun(C0widths, xset); % locations in xset with initial coral patches (where initial coral cover = C0high)

% for icchoice = 5 (sinusoidal initial conditions)
 ampC0 = (C0high-C0low)/2; % amplitude of coral patches
 ampM0 = (M0high-M0low)/2; % amplitude of macroalgal patches
 period0 = 0.4; % controls the angular frequency of the sine waves

% peak characteristics
pkthresh = 0.05; % min prominence that a peak has to have to count as a patch
dthresh = 0.25*len; % threshold distance from edge before a peak gets considered as a patch
b1 = xset(1) + dthresh; % lower boundary of landscape for peak consideration
b2 = xset(end)-dthresh; % upper boundary of landscape for peak consideration

% get the indeces of xset corresponding to these boundaries (will use these for intervals to 
% take spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1))); % approximate index = index corresponding to value of xset that is closest to b1
b2i = find(abs(xset-b2)==min(abs(xset-b2)));


%  values of fishing pressure at which to simulate the PDE model
fset21 = linspace(0.13, 0.19, 20); 

% holding arrays for storing model output (at each time point/spatial location) at each fishing pressure
Cruns = NaN(length(tset), length(xset),length(fset21)); % coral
Mruns = NaN(length(tset), length(xset),length(fset21)); % total macroalgae (MI + Mv)
Hruns = NaN(length(tset), length(xset),length(fset21)); % herbivores

% also record average abundance for each fishing pressure
Cmeans = NaN(length(fset21));
Mmeans = NaN(length(fset21));
Hmeans = NaN(length(fset21));

% record the maxes and mins of the coral and macroalgal peaks as well
Cmxs = NaN(length(fset21));
Cmns = NaN(length(fset21));
Mmxs = NaN(length(fset21));
Mmns = NaN(length(fset21));

summ10 = 1; % 1 = record metrics from 3 peaks closest to center, 0 = record all peaks
pksumm = NaN(6,3,length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

    tic % use tictoc to time simulations
    for i = 1:length(fset21) % for each element from 1 to length of fset21

    ftest = fset21(i); % set fishing pressure to the ith element of fset21

     % run PDE
     [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record full results
    Cruns(:, :,i) = solij(:,:,2);
    Mruns(:, :, i) = solij(:,:,1)+ solij(:,:,4);
    Hruns(:, :, i) = solij(:,:,3);

    % record spatial averages at final time point (averaged over locations
    % between xset(b1i) and xset(b2i))
    Cmeans(i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i) = mean(solij(end, b1i:b2i, 3));

    % record max and min coral and macroalgal cover (between xset(b1i) and
    % xset(b2i))
    Cmxs(i) = max(solij(end, b1i:b2i, 2));
    Cmns(i) = min(solij(end, b1i:b2i, 2));
    Mmxs(i) = max(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Mmns(i) = min(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));

     % calculate the peak summary metrics
      Cvalsijk = solij(end, :, 2); % coral cover at each spatial location at final timepoint 
      % (note don't need to subset b1i:b2i because peakfun does that)
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4); % total macroalgal cover at each spatial location at final timepoint
      % use peakfun to calculate peak summary metrics
      [npks0, npks, pklambdas,pkwidths,pkproms,pkheights] = peakfun(Cvalsijk,Mvalsijk,summ10,xset, pkthresh, b1, b2);

            % record these
            pksumm(6,1,i) = npks0/(b2-b1); % total peak density
            pksumm(5,1,i) = npks/(b2-b1); % density of peaks for which C > M
            pksumm(1,1:length(pklambdas),i) = pklambdas; % peak wavelengths
            pksumm(2,1:length(pkwidths),i) = pkwidths; % peak widths
            pksumm(3,1:length(pkproms),i) = pkproms; % peak prominances
            pksumm(4,1:length(pkheights),i) = pkheights; % peak absolute heights

   
    end % end for loop

    toc % 40 seconds
 

%% plot initial conditions 
% to check initialization
% figure(1)
% plot(xset, Cruns(1,:,1,1,1), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
% ylim([0 1])
% xlabel('Location','FontSize',22)
% ylabel('Prop. cover','FontSize',22)
% hold on 
% plot(xset, Mruns(1,:,1,1,1), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
% hold off
% legend('Coral','Macroalgae', 'Location','Northwest','NumColumns',2)
% lgd = legend;
% lgd.FontSize = 14;

%% get ensemble max/mins (from FigS6.m)
% load the ensemble max and min coral and macroalgal cover at each fishing
% pressure across a range of initial conditions (calculated in FigS6.m and
% stored in FigS6.mat)
load('code output/FigS6.mat','Mmxall','Mmnall', 'Cmxall', 'Cmnall')


%% plot bifurcation diagram for macroalgae (Fig. 1a) and coral (Fig. 1b)

fpts = [3, 7, 12, 15]; % elements of fset21 to highlight in the figure

fcol = [0.3020 0.2 0.9333]; % color for the fishing pressures

% polygon around region of bistability
pgon = polyshape([flow flow fup fup],[2 -1 -1 2]);


figure(1)
x0=10;
y0=10;
width=500;
height=950;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 1);
t.TileSpacing = 'compact';
nexttile
% start with bifurcation diagram for macroalgae
plot(pgon,'FaceColor','black','FaceAlpha',0.025) % polygon around bistability region
xlabel(t,'Fishing pressure (f)','FontSize',22) % t for shared label across tiles
ylabel('Equilibrium macroalgal cover','FontSize',22)
title('a) Macroalgae equilibria', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
ylim([0 0.85])
xlim([0.12, max(fset)])
hold on
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) 
text(0.1715, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% add the ensemble max and mins and plot them as a filled polygon
rng = find(Mmxall(:,1)-Mmnall(:,1)> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [Mmnall(rng,1)', fliplr(Mmxall(rng,1)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.15, 'EdgeColor', Mcol, 'EdgeAlpha', 0.5);
% add circle around f value for which spatial distribution is shown
plotj = fpts(1); % f = 0.136
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.002, Mmeans(plotj) + 0.04, 'f = 0.136', 'Color', fcol,'FontSize', 14)
plotj = fpts(2); % f = 0.149
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, Mmeans(plotj) + 0.04, 'f = 0.149', 'Color', fcol,'FontSize', 14)
plotj = fpts(3); % f = 0.165
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.009, Mmeans(plotj) + 0.04, 'f = 0.165', 'Color', fcol,'FontSize', 14)
plotj = fpts(4); % f = 0.174
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)+0.0013, Mmeans(plotj) - 0.025, 'f = 0.174', 'Color', fcol,'FontSize', 14)
hold off
hold on
% make legend
lg{1} = plot(nan, 'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5);
lg{2} = plot(nan, 'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5);
lg{3} = plot(nan, '.','MarkerSize',35,'Color', Mcol);
lg{4} = fill(nan, nan, Mcol, 'FaceAlpha',0.2, 'EdgeColor', Mcol, 'EdgeAlpha', 0.5);
legend([lg{:}],{'ODE, stable', 'ODE, unstable','PDE, means','PDE, range'}, 'Location', 'northwest')
hold off
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
% now repeat for coral
nexttile
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylabel('Equilibrium coral cover','FontSize',22)
title('b) Coral equilibria', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
ylim([0 0.85])
xlim([0.12, max(fset)])
hold on
plot(fset, Clows,'Color', Ccol, "LineStyle","-", 'LineWidth', 2.5) 
plot(fset, Cups,'Color', Ccol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Cmids,'Color', Ccol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Clows,'Color', Ccol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Cmeans, '.','MarkerSize',30,'Color', Ccol)
% add the ensemble max/mins
rng = find(Cmxall(:,1)-Cmnall(:,1)> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [Cmnall(rng,1)', fliplr(Cmxall(rng,1)')];
    fill(btwx, btwy, Ccol, 'FaceAlpha',0.15, 'EdgeColor', Ccol, 'EdgeAlpha', 0.5);
% add circle around f values for which the spatial distributions are shown
plotj = fpts(1); 
plot(fset21(plotj), Cmeans(plotj), '.','MarkerSize',30,'Color', Ccol)
plot(fset21(plotj), Cmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
plotj = fpts(2); 
plot(fset21(plotj), Cmeans(plotj), '.','MarkerSize',30,'Color', Ccol)
plot(fset21(plotj), Cmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
plotj = fpts(3); 
plot(fset21(plotj), Cmeans(plotj), '.','MarkerSize',30,'Color', Ccol)
plot(fset21(plotj), Cmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
plotj = fpts(4); 
plot(fset21(plotj), Cmeans(plotj), '.','MarkerSize',30,'Color', Ccol)
plot(fset21(plotj), Cmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
hold off

%% plot spatial distributions (Fig 1c)

% 4x1 panel plot
figure(2)
x0=10;
y0=10;
width=600;
height=550;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 1);
t.TileSpacing = 'compact';
nexttile
plotj = fpts(4); % f value
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Ccol)
text(-90, 1.5,'f = 0.174','FontSize',14, 'Color', fcol)
ylim([-0.01 1.75])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
title('c) Spatial distributions', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off
legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);
nexttile
plotj = fpts(3);% f value
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Ccol)
text(-90, 1.5,'f = 0.165','FontSize',14, 'Color', fcol)
ylim([-0.01 1.75])
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
hold off
hold on
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off
nexttile
plotj = fpts(2); % f value
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Ccol)
text(-90, 1.5,'f = 0.149','FontSize',14, 'Color', fcol)
ylim([-0.01 1.75])
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
hold off
hold on
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off
nexttile
plotj = fpts(1); % f value
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Ccol)
text(-90, 1.5,'f = 0.136','FontSize',14, 'Color', fcol)
ylim([-0.01 1.75])
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
hold off
hold on
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off

%% full space-time plots at each fishing pressure

% figure(3)
% x0=10;
% y0=10;
% width=1000;
% height=300;
% set(gcf,'position',[x0,y0,width,height])
% t=tiledlayout(1, 4);
% t.TileSpacing = 'compact';
% nexttile
% plotj = fpts(1);
% surf(Cruns(1:100, b1i:b2i,plotj),'FaceAlpha',1, 'EdgeColor','none')
% colormap(CMmap)
% title('a) f = 0.092','FontSize',16)
% xlabel(t,'Location','FontSize',18) % t for shared label
% ylabel(t,'Time','FontSize',18) % t for shared label
% view(0,90)
% nexttile
% plotj = fpts(2);
% surf(Cruns(1:100, b1i:b2i,plotj),'FaceAlpha',1, 'EdgeColor','none')
% title('b) f = 0.101','FontSize',16)
% view(0,90)
% colormap(CMmap)
% nexttile
% plotj = fpts(3);
% surf(Cruns(1:100, b1i:b2i,plotj),'FaceAlpha',1, 'EdgeColor','none')
% title('c) f = 0.109','FontSize',16)
% view(0,90)
% colormap(CMmap)
% nexttile
% plotj = fpts(4);
% surf(Cruns(1:100, b1i:b2i,plotj),'FaceAlpha',1, 'EdgeColor','none')
% title('d) f = 0.117','FontSize',16)
% view(0,90)
% colormap(CMmap)
% colorbar
% cb = colorbar;
% cbL = ylabel(cb,{'Coral'; 'cover'},'FontSize', 18);     
% set(cbL,'Rotation',0);


%% plot peak characteristics (Figure S1)
C1 = Ccol;

pgon = polyshape([flow flow fup fup],[30 -1 -1 30]); % polygon around region of bistability


figure(4)
x0=10;
y0=10;
width=400;
height=900;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 1);
t.TileSpacing = 'compact';
nexttile
% start with patch density
plot(fset21, squeeze(pksumm(5,1,:)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
hold on
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
hold off
xline([flow fup]) % bistability region
xlim([min(fset21) max(fset21)])
ylim([0 1.05*max(squeeze(pksumm(5,1,:)))])
ylabel({'Coral patch';'density'},'FontSize',22)
text(0.172, 0.09, 'Bistable', 'Color', 'black','FontSize', 16)
nexttile
% now patch widths
plot(fset21, squeeze(pksumm(2,1,:)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(fset21) max(fset21)])
ylim([0 26])
xline([flow fup])
ylabel({'Coral patch';'width'},'FontSize',22)
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
plot(fset21, squeeze(pksumm(2,2,:)),'Color',C1,"LineStyle","none",'Marker','.', 'LineWidth', 2.5)
plot(fset21, squeeze(pksumm(2,3,:)),'Color',C1,"LineStyle",":", 'LineWidth', 2.5)
hold off
nexttile
% now patch height (max coral cover in patch)
plot(fset21, squeeze(pksumm(4,1,:)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(fset21) max(fset21)])
ylim([0 0.85])
xline([flow fup]) 
ylabel({'Max coral';'cover in patch'},'FontSize',22)
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
plot(fset21, squeeze(pksumm(4,2,:)),'Color',C1,"LineStyle","none",'Marker','.', 'LineWidth', 2.5)
plot(fset21, squeeze(pksumm(4,3,:)),'Color',C1,"LineStyle",":", 'LineWidth', 2.5)
hold off
nexttile
% now mean coral cover
plot(fset21, Cmeans, 'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(fset21) max(fset21)])
ylim([0 0.85])
xline([flow fup]) 
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Mean';'coral cover'},'FontSize',22)
hold on
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
hold off


















