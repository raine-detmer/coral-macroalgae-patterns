% README: code for making Fig S22 (repeating figure 1 for scenario in which total herbivore
% abundance remains constant)


%% set up

% plot colors
Mcol = [0.4667 0.6745 0.1882]; % macroalgae
Ccol = [0.3020 0.7451 0.9333]; % coral
Hcol = [0.9294 0.6941 0.1255]; % herbivores

% region of bistability
flow = 0.1674; % lower tipping point (calculated in Fig2.m)
fup = 0.1878; % upper tipping point

%% Briggs: ODE (non-spatial) bifurcation diagram

% make the non-spatial bifurcation diagram

% define the symbols
syms Mi C H Mv

% define parameters
gTC = 0.1; 
gamma = 0.4; 
gTI = 0.4;
rM = 0.5; 
gTV = 0.2;
dv = 2; 
omega = 2; 
di = 0.4; 
phiC = 0.01; 
dC = 0.02;
phiM = 0.01; 

% herbivore parameters
rH = 0.2; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
f = 0; % herbivore fishing pressure
phiH = 0.05;

% set of fishing pressure values
fset = linspace(0.12, 0.24, 120);

% holding vectors for equilibrium values
Cstars = NaN(length(fset), 4); % coral
Mstars = NaN(length(fset), 4); % macroalgae (vuln + invuln)

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

for i = 1:length(fset)%for each fishing pressure
    
    fi = fset(i);
    H = ((rH-fi) + sqrt((rH-fi)^2 + 4*dH*phiH))/(2*dH);

    % get the eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    %eq3i = rH*H-dH*H*H-fi*H ==0; % no herbivore dynamics
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq4i],[Mi,C, Mv], [0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C and M cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

% process results

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region


% use vertcat to concatenate vertical vectors
% for f on x axis:
Cups = vertcat(Cstars(1:bstart-1, 1), Cstars(bstart:bend, 3), Cstars(bend+1:end, 4)); % need to make sure the length stays the same so concatenate with NaNs 
Cmids = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:bend, 2), Cstars(bend+1:end, 4)); % need to make sure the length stays the same so concatenate with NaNs 
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 1));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still (high and low
% coral eq)


%% Briggs: PDE bifurcation diagram

% PDE parameters
diffs = [0.05,0.05,0.25, 0]; % diffusion rates of MI, C, H, and Mv
taxisM = 0; 
taxisC = -0.25; % taxis rate toward coral
taxisT = 0;

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges

% space
len = 400;
xset = linspace(-len/2,len/2,800);

% time
t_end = 3*50000;
tset = linspace(0,t_end,2*2500); 


% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function


C0high = 0.85; % coral cover in initial coral patches
C0low = 0.05; % coral cover in initial macroalgal patches
M0high = 0.85; % total macroalgal cover in initial macroalgal patches
M0low = 0.05; % total macroalgal cover in initial coral patches

% for icchoice = 3
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4
C0widths = round(length(xset)/16);  % step widths
initC = stepfun(C0widths, xset); % elements of xset where coral is initially high


% for icchoice = 5
 ampC0 = (C0high-C0low)/2;
 ampM0 = (M0high-M0low)/2;
 period0 = 0.4;

% peak characteristics
pkthresh = 0.05; % min prominence that a peak has to have to count
dthresh = 0.25*len; % threshold distance from edge before a peak gets considered
b1 = xset(1) + dthresh; % lower boundary for peak consideration
b2 = xset(end)-dthresh; % upper boundary for peak consideration

% get the indeces of these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));


%  values of fishing pressure
fset21 = linspace(0.11, 0.21, 20); 

% holding arrays
Cruns = NaN(length(tset), length(xset),length(fset21));
Mruns = NaN(length(tset), length(xset),length(fset21));
Hruns = NaN(length(tset), length(xset),length(fset21));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(fset21));
Mmeans = NaN(length(fset21));
Hmeans = NaN(length(fset21));


summ10 = 1; % 1 = record metrics from 3 peaks closest to center of landscape, 0 = record all peaks

    tic
   for i = 1:length(fset21) % for each fishing pressure
  % for i = 16:20

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record full results
    Cruns(:, :,i) = solij(:,:,2);
    Mruns(:, :, i) = solij(:,:,1)+ solij(:,:,4);
    Hruns(:, :, i) = solij(:,:,3);

    % record spatial averages at final time point
    Cmeans(i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i) = mean(solij(end, b1i:b2i, 3));

   
    end 

    toc % 35 seconds, 168 with random IC
 



%% plot bifurcation diagram for macroalgae (Fig. 1a) and coral (Fig. 1b)

fpts = [2+2, 5+2, 10+2, 17+2]; % elements of fset21 to highlight in the figure


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
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22) % t for shared label
ylabel('Equilibrium macroalgal cover','FontSize',22)
title('a) Macroalgae equilibria', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
ylim([0 0.85])
xlim([min(fset), max(fset)])
hold on
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
text(0.1685, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% add circle around f value for which spatial distribution is shown
plotj = fpts(1); % f = 0.0921
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.002, Mmeans(plotj) + 0.035, 'f = 0.126', 'Color', fcol,'FontSize', 14)
plotj = fpts(2); % f = 0.1005
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.018, Mmeans(plotj) + 0.035, 'f = 0.142', 'Color', fcol,'FontSize', 14)
plotj = fpts(3); % f = 0.1089
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.018, Mmeans(plotj) + 0.03, 'f = 0.168', 'Color', fcol,'FontSize', 14)
plotj = fpts(4); % f = 0.1174
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
%text(fset21(plotj)+0.0013, Mmeans(plotj) - 0.025, 'f = 0.133', 'Color', fcol,'FontSize', 14)
text(fset21(plotj)-0.015, Mmeans(plotj) + 0.035, 'f = 0.205', 'Color', fcol,'FontSize', 14)
hold off
hold on
% make legend
lg{1} = plot(nan, 'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5);
lg{2} = plot(nan, 'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5);
lg{3} = plot(nan, '.','MarkerSize',35,'Color', Mcol);
%lg{4} = fill(nan, nan, Mcol, 'FaceAlpha',0.2, 'EdgeColor', Mcol, 'EdgeAlpha', 0.5);
legend([lg{1:3}],{'ODE, stable', 'ODE, unstable','PDE, means'}, 'Location', 'northwest')
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
xlim([min(fset), max(fset)])
hold on
plot(fset, Clows,'Color', Ccol, "LineStyle","-", 'LineWidth', 2.5) 
plot(fset, Cups,'Color', Ccol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Cmids,'Color', Ccol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Clows,'Color', Ccol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Cmeans, '.','MarkerSize',30,'Color', Ccol)
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
text(-90, 1.75,'f = 0.205','FontSize',14, 'Color', fcol)
ylim([-0.01 2])
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
text(-90, 1.75,'f = 0.168','FontSize',14, 'Color', fcol)
ylim([-0.01 2])
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
hold off
hold on
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off
nexttile
plotj = fpts(2); % f value
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Ccol)
text(-90, 1.75,'f = 0.142','FontSize',14, 'Color', fcol)
ylim([-0.01 2])
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
hold off
hold on
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off
nexttile
plotj = fpts(1); % f value
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Ccol)
text(-90, 1.75,'f = 0.126','FontSize',14, 'Color', fcol)
ylim([-0.01 2])
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Mcol)
hold off
hold on
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, plotj), 'LineWidth',2, 'Color', Hcol)
hold off
