% README: code for making Figure S4 and Figure S5


%% set up

% plot colors
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];
Hcol = [0.9294 0.6941 0.1255];

% colormap for space-time plots
vec = [100; 0]; % nodes at which to place the colors
raw = [Mcol; Ccol];
N = 256;
CMmap = interp1(vec, raw, linspace(100,0, N),'pchip');


flow = 0.1111; % lower tipping point (calculated in Fig2.m)
fup = 0.1229; % upper tipping point

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
rH = 0.2;%0.1; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
f = 0.08; % herbivore fishing pressure

% set of fishing values
fset = linspace(0.05, 0.15, 100);

% holding vector of eq values
Cstars = NaN(length(fset), 4);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 4);

Mistars = NaN(length(fset), 4);
Mvstars = NaN(length(fset), 4);

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

for i = 1:length(fset)%for each element of gset
    % get the fishing pressure
    fi = fset(i);

    % solve the equations
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

% process results

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% use vertcat to concatenate vertical vectors to get the upper, middle, and
% lower equilibria
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(4,)
Cmids = Cstars(:, 3);
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still

%% Briggs: PDE bifurcation diagram

% PDE parameters
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5;%0; % taxis rate toward coral
taxisT = 0;

diric = 0; % 0 = Neumann boundariess for constant habitat. 1 = Dirichlet boundaries for loss at the edges

% space
len = 400;
xset = linspace(-len/2,len/2,800);

% time
t_end = 3*50000;
tset = linspace(0,t_end,2*2500); 


% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85;
C0low = 0.05;%0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 3
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4
C0widths = round(length(xset)/16);  % step widths
initC = stepfun(C0widths, xset); 

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
fset21 = linspace(0.09, 0.13, 20); % for higher bistability region

% initial conditions: patch widths
wset = round([length(xset)/2, length(xset)/16, length(xset)/64]);

% holding arrays
Cruns = NaN(length(tset), length(b1i:b2i),length(fset21), length(wset));
Mruns = NaN(length(tset), length(b1i:b2i),length(fset21), length(wset));
Hruns = NaN(length(tset), length(b1i:b2i),length(fset21), length(wset));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(fset21), length(wset));
Mmeans = NaN(length(fset21), length(wset));
Hmeans = NaN(length(fset21), length(wset));

% record the maxes and mins as well
Cmxs = NaN(length(fset21), length(wset));
Cmns = NaN(length(fset21), length(wset));
Mmxs = NaN(length(fset21), length(wset));
Mmns = NaN(length(fset21), length(wset));

tic
for j = 1:length(wset) % for each initial condition

    C0widths = wset(j);
    initC = stepfun(C0widths, xset); 

    for i = 1:length(fset21) % for each fishing pressure
    %for i = 1:5

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(:, :,i,j) = solij(:,b1i:b2i,2);
    Mruns(:, :, i,j) = solij(:,b1i:b2i,1)+ solij(:,b1i:b2i,4);
    Hruns(:, :, i,j) = solij(:,b1i:b2i,3);

    % record spatial averages at final time point
    Cmeans(i,j) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i,j) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i,j) = mean(solij(end, b1i:b2i, 3));

    Cmxs(i,j) = max(solij(end, b1i:b2i, 2));
    Cmns(i,j) = min(solij(end, b1i:b2i, 2));
    Mmxs(i,j) = max(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Mmns(i,j) = min(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
   
    end 

    

end 

toc % 164 seconds




%% plot bifurcation diagrams (Fig. S4)
% couldn't make the figure tall enough when all together so plot each row
% separately

fpts = [6, 10]; % fishing pressures for which to show spatial distributions

fcol = [0.3020 0.2 0.9333];

pgon = polyshape([flow flow fup fup],[2 -1 -1 2]);


figure(1)
x0=10;
y0=10;
width=1200;
height=325;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 7); % (rows, columns)
t.TileSpacing = 'compact';

pp = 1;
MmeansP = Mmeans(:,pp);
MmxsP = Mmxs(:,pp);
MmnsP = Mmns(:,pp);
plotj = fpts(1);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax1 = nexttile(1,[1,1]); % nexttile(x, [r,c]) means put the upper left corner of the
% axes in tile x and then make the plot span r rows and c columns
plot(ax1,xset(b1i:b2i),HrunsP(1,:), 'LineWidth',2, 'Color', Hcol)
xticklabels([]) % get rid of the x axis tick labels
ylim([0 1.25])
hold on
plot(xset(b1i:b2i),CrunsP(1,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),MrunsP(1,:), 'LineWidth',2, 'Color', Mcol)
hold off
title('a) Initial patch width = 1/2 total space', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax2 = nexttile(2,[4,3]); % nexttile(x, [r,c]) means put the upper left corner of the
plot(ax2, pgon,'FaceColor','black','FaceAlpha',0.025)
xline([flow fup]) % bistability region
ylim([0 0.85])
xlim([0.08, 0.14])
hold on
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
text(0.112, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% add the max and mins from the other ICs
MmxsP = Mmxs(:,2);
MmnsP = Mmns(:,2);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 0.6);
MmxsP = Mmxs(:,3);
MmnsP = Mmns(:,3);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 0.6);
% NEWER UPDATE: add the maxes and mins
MmxsP = Mmxs(:,pp);
MmnsP = Mmns(:,pp);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 1, 'LineWidth', 1.5);
%plot the PDE means
% other ICs
MmeansP = Mmeans(:,2);
plot(fset21, MmeansP, 'o','MarkerSize',8,'Color', Mcol, 'LineWidth', 1.2)
MmeansP = Mmeans(:,3);
plot(fset21, MmeansP, 'o','MarkerSize',8,'Color', Mcol, 'LineWidth', 1.2)
% current IC
MmeansP = Mmeans(:,pp);
plot(fset21, MmeansP, '.','MarkerSize',30,'Color', Mcol)
% add circle around f value for which spatial distribution is shown
plotj = fpts(1); % f value = 0.1005
plot(fset21(plotj), MmeansP(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), MmeansP(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, MmeansP(plotj) + 0.03, 'f = 0.101', 'Color', fcol,'FontSize', 14)
plotj = fpts(2); % f value =  0.1089
plot(fset21(plotj), MmeansP(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), MmeansP(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, MmeansP(plotj) + 0.03, 'f = 0.109', 'Color', fcol,'FontSize', 14)
hold off
% add the legend
hold on
lg{1} = plot(nan, 'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5);
lg{2} = plot(nan, 'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5);
lg{3} = plot(nan, '.','MarkerSize',35,'Color', Mcol);
lg{4} = fill(nan, nan, Mcol, 'FaceAlpha',0.2, 'EdgeColor', Mcol, 'EdgeAlpha', 0.5);
legend([lg{:}],{'ODE, stable', 'ODE, unstable','PDE, means','PDE, range'}, 'Location', 'northwest')
hold off
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
title(lgd,{'Macroalgal equilibria'})
% nexttile: equilibrium spatial distributions
plotj = fpts(2);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax3 = nexttile(5,[2,3]); 
plot(ax3,xset(b1i:b2i),MrunsP(end,:), 'LineWidth',2, 'Color', Mcol)
ylim([0 1.6])
text(-90, 1.4,'f = 0.109','FontSize',14, 'Color', fcol)
text(-40, 1.4,'Coral','FontSize',14, 'Color', Ccol)
text(-10, 1.4,'Macroalgae','FontSize',14, 'Color', Mcol)
text(40, 1.4,'Herbivores','FontSize',14, 'Color', Hcol)
hold on
plot(xset(b1i:b2i),CrunsP(end,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP(end,:), 'LineWidth',2, 'Color', Hcol)
hold off
% nexttile
plotj = fpts(1);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax4 = nexttile(19,[2,3]); 
plot(ax4,xset(b1i:b2i),MrunsP(end,:), 'LineWidth',2, 'Color', Mcol)
ylim([0 1.6])
text(-90, 1.4,'f = 0.101','FontSize',14, 'Color', fcol)
hold on
plot(xset(b1i:b2i),CrunsP(end,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP(end,:), 'LineWidth',2, 'Color', Hcol)
hold off

figure(2)
x0=10;
y0=10;
width=1200;
height=325;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 7); % (rows, columns)
t.TileSpacing = 'compact';

pp = 2;
MmeansP = Mmeans(:,pp);
MmxsP = Mmxs(:,pp);
MmnsP = Mmns(:,pp);
plotj = fpts(1);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax1 = nexttile(1,[1,1]); % nexttile(x, [r,c]) means put the upper left corner of the
% axes in tile x and then make the plot span r rows and c columns
plot(ax1,xset(b1i:b2i),HrunsP(1,:), 'LineWidth',2, 'Color', Hcol)
xticklabels([]) % get rid of the x axis tick labels
ylim([0 1.25])
hold on
plot(xset(b1i:b2i),CrunsP(1,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),MrunsP(1,:), 'LineWidth',2, 'Color', Mcol)
hold off
title('b) Initial patch width = 1/16 total space', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax2 = nexttile(2,[4,3]); 
plot(ax2, pgon,'FaceColor','black','FaceAlpha',0.025)
xline([flow fup]) % bistability region
ylim([0 0.85])
xlim([0.08, 0.14])
hold on
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% add the max and mins from the other ICs
MmxsP = Mmxs(:,1);
MmnsP = Mmns(:,1);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 0.6);
MmxsP = Mmxs(:,3);
MmnsP = Mmns(:,3);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 0.6);
% add the maxes and mins from the current ICs
MmxsP = Mmxs(:,pp);
MmnsP = Mmns(:,pp);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 1, 'LineWidth', 1.5);
%plot the PDE means
% other ICs
MmeansP = Mmeans(:,1);
plot(fset21, MmeansP, 'o','MarkerSize',8,'Color', Mcol, 'LineWidth', 1.2)
MmeansP = Mmeans(:,3);
plot(fset21, MmeansP, 'o','MarkerSize',8,'Color', Mcol, 'LineWidth', 1.2)
% current IC
MmeansP = Mmeans(:,pp);
plot(fset21, MmeansP, '.','MarkerSize',30,'Color', Mcol)
% add circle around f value for which spatial distribution is shown
plotj = fpts(1); % f value = 0.1005
plot(fset21(plotj), MmeansP(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), MmeansP(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, MmeansP(plotj) + 0.03, 'f = 0.101', 'Color', fcol,'FontSize', 14)
plotj = fpts(2); % f value =  0.1089
plot(fset21(plotj), MmeansP(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), MmeansP(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, MmeansP(plotj) + 0.03, 'f = 0.109', 'Color', fcol,'FontSize', 14)
hold off

% nexttile: equilibrium spatial distributions
plotj = fpts(2);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax3 = nexttile(5,[2,3]); 
plot(ax3,xset(b1i:b2i),MrunsP(end,:), 'LineWidth',2, 'Color', Mcol)
ylim([0 1.6])
text(-90, 1.4,'f = 0.109','FontSize',14, 'Color', fcol)
hold on
plot(xset(b1i:b2i),CrunsP(end,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP(end,:), 'LineWidth',2, 'Color', Hcol)
hold off
% nexttile
plotj = fpts(1);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax4 = nexttile(19,[2,3]); 
plot(ax4,xset(b1i:b2i),MrunsP(end,:), 'LineWidth',2, 'Color', Mcol)
ylim([0 1.6])
text(-90, 1.4,'f = 0.101','FontSize',14, 'Color', fcol)
hold on
plot(xset(b1i:b2i),CrunsP(end,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP(end,:), 'LineWidth',2, 'Color', Hcol)
hold off


figure(3)
x0=10;
y0=10;
width=1200;
height=325;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 7); % (rows, columns)
t.TileSpacing = 'compact';

pp = 3;
MmeansP = Mmeans(:,pp);
MmxsP = Mmxs(:,pp);
MmnsP = Mmns(:,pp);
plotj = fpts(1);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax1 = nexttile(1,[1,1]); % nexttile(x, [r,c]) means put the upper left corner of the
% axes in tile x and then make the plot span r rows and c columns
plot(ax1,xset(b1i:b2i),HrunsP(1,:), 'LineWidth',2, 'Color', Hcol)
xticklabels([]) % get rid of the x axis tick labels
ylim([0 1.25])
hold on
plot(xset(b1i:b2i),CrunsP(1,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),MrunsP(1,:), 'LineWidth',2, 'Color', Mcol)
hold off
title('c) Initial patch width = 1/64 total space', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax2 = nexttile(2,[4,3]);
plot(ax2, pgon,'FaceColor','black','FaceAlpha',0.025)
xline([flow fup]) % bistability region
ylim([0 0.85])
xlim([0.08, 0.14])
hold on
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% add the max and mins from the other ICs
MmxsP = Mmxs(:,1);
MmnsP = Mmns(:,1);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 0.6);
MmxsP = Mmxs(:,2);
MmnsP = Mmns(:,2);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 0.6);
% add the maxes and mins from the current IC
MmxsP = Mmxs(:,pp);
MmnsP = Mmns(:,pp);
rng = find(MmxsP-MmnsP> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [MmnsP(rng)', fliplr(MmxsP(rng)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.08, 'EdgeColor', Mcol, 'EdgeAlpha', 1, 'LineWidth', 1.5);
%plot the PDE means
% other ICs
MmeansP = Mmeans(:,1);
plot(fset21, MmeansP, 'o','MarkerSize',8,'Color', Mcol, 'LineWidth', 1.2)
MmeansP = Mmeans(:,2);
plot(fset21, MmeansP, 'o','MarkerSize',8,'Color', Mcol, 'LineWidth', 1.2)
% current IC
MmeansP = Mmeans(:,pp);
plot(fset21, MmeansP, '.','MarkerSize',30,'Color', Mcol)
% add circle around f value for which spatial distribution is shown
plotj = fpts(1); % f value = 0.1005
plot(fset21(plotj), MmeansP(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), MmeansP(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, MmeansP(plotj) + 0.03, 'f = 0.101', 'Color', fcol,'FontSize', 14)
plotj = fpts(2); % f value =  0.1089
plot(fset21(plotj), MmeansP(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), MmeansP(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, MmeansP(plotj) + 0.03, 'f = 0.109', 'Color', fcol,'FontSize', 14)
hold off

% nexttile: equilibrium spatial distributions
plotj = fpts(2);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax3 = nexttile(5,[2,3]);
plot(ax3,xset(b1i:b2i),MrunsP(end,:), 'LineWidth',2, 'Color', Mcol)
ylim([0 1.6])
text(-90, 1.4,'f = 0.109','FontSize',14, 'Color', fcol)
hold on
plot(xset(b1i:b2i),CrunsP(end,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP(end,:), 'LineWidth',2, 'Color', Hcol)
hold off
% nexttile
plotj = fpts(1);
MrunsP = Mruns(:,:,plotj, pp);
CrunsP = Cruns(:,:,plotj, pp);
HrunsP = Hruns(:,:,plotj, pp);
ax4 = nexttile(19,[2,3]); 
plot(ax4,xset(b1i:b2i),MrunsP(end,:), 'LineWidth',2, 'Color', Mcol)
ylim([0 1.6])
text(-90, 1.4,'f = 0.101','FontSize',14, 'Color', fcol)
hold on
plot(xset(b1i:b2i),CrunsP(end,:), 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP(end,:), 'LineWidth',2, 'Color', Hcol)
hold off


%% surface plots


fpts = [6, 10]; % fishing pressures to show

tmx = 500; % number of timepoints to plot

tset([1, 100, 200, 300, 400, 500]); 
% 0    0.2971    0.5971    0.8972    1.1972    1.4973

x0=10;
y0=10;
width=800;
height=900;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(3, 2);
t.TileSpacing = 'compact';
nexttile
pp = 1;
plotj = fpts(1);
CrunsP = Cruns(1:tmx,:,plotj, pp);
surf(CrunsP,'FaceAlpha',1, 'EdgeColor','none')
colormap(CMmap)
view(0,90)
title({'f = 0.101'; 'Initial patch width = 1/2 total space'},'FontSize',14)
%xticklabels([-100 -50 0 50 100]) % make the xaxis labels match the distribution plots
xticklabels([ -100  -75  -50  -25    0   25   50   75  100]) % make the xaxis labels match the distribution plots
yticklabels([ 0*10^4    0.3*10^4   0.6*10^4    0.9*10^4    1.2*10^4    1.5*10^4]) % make the yaxis labels match actual timepoints
nexttile
plotj = fpts(2);
CrunsP = Cruns(1:tmx,:,plotj, pp);
surf(CrunsP,'FaceAlpha',1, 'EdgeColor','none')
colormap(CMmap)
view(0,90)
title({'f = 0.109'; 'Initial patch width = 1/2 total space'},'FontSize',14)
xticklabels([ -100  -75  -50  -25    0   25   50   75  100]) % make the xaxis labels match the distribution plots
yticklabels([ 0*10^4    0.3*10^4   0.6*10^4    0.9*10^4    1.2*10^4    1.5*10^4]) % make the yaxis labels match actual timepoints
nexttile
pp = 2;
plotj = fpts(1);
CrunsP = Cruns(1:tmx,:,plotj, pp);
surf(CrunsP,'FaceAlpha',1, 'EdgeColor','none')
view(0,90)
title({'Initial patch width = 1/16 total space'},'FontSize',14)
xticklabels([ -100  -75  -50  -25    0   25   50   75  100]) % make the xaxis labels match the distribution plots
yticklabels([ 0*10^4    0.3*10^4   0.6*10^4    0.9*10^4    1.2*10^4    1.5*10^4]) % make the yaxis labels match actual timepoints
nexttile
plotj = fpts(2);
CrunsP = Cruns(1:tmx,:,plotj, pp);
surf(CrunsP,'FaceAlpha',1, 'EdgeColor','none')
view(0,90)
title('Initial patch width = 1/16 total space','FontSize',14)
xticklabels([ -100  -75  -50  -25    0   25   50   75  100]) % make the xaxis labels match the distribution plots
yticklabels([ 0*10^4    0.3*10^4   0.6*10^4    0.9*10^4    1.2*10^4    1.5*10^4]) % make the yaxis labels match actual timepoints
colorbar
cb = colorbar;
cbL = ylabel(cb,{'Coral'; 'cover'},'FontSize', 14);     
set(cbL,'Rotation',0);
nexttile
pp = 3;
plotj = fpts(1);
CrunsP = Cruns(1:tmx,:,plotj, pp);
surf(CrunsP,'FaceAlpha',1, 'EdgeColor','none')
view(0,90)
title('Initial patch width = 1/64 total space','FontSize',14)
xticklabels([ -100  -75  -50  -25    0   25   50   75  100]) % make the xaxis labels match the distribution plots
yticklabels([ 0*10^4    0.3*10^4   0.6*10^4    0.9*10^4    1.2*10^4    1.5*10^4]) % make the yaxis labels match actual timepoints
nexttile
plotj = fpts(2);
CrunsP = Cruns(1:tmx,:,plotj, pp);
surf(CrunsP,'FaceAlpha',1, 'EdgeColor','none')
view(0,90)
title('Initial patch width = 1/64 total space','FontSize',14)
xticklabels([ -100  -75  -50  -25    0   25   50   75  100]) % make the xaxis labels match the distribution plots
yticklabels([ 0*10^4    0.3*10^4   0.6*10^4    0.9*10^4    1.2*10^4    1.5*10^4]) % make the yaxis labels match actual timepoints
xlabel(t, 'Location','FontSize',18) % t for shared label
ylabel(t, 'Time','FontSize',18) % t for shared label

