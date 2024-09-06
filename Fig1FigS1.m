% README: code for making Figure 1 and Figure S1


%% set up

% plot colors
Mcol = [0.4667 0.6745 0.1882]; % macroalgae
Ccol = [0.3020 0.7451 0.9333]; % coral
Hcol = [0.9294 0.6941 0.1255]; % herbivores

% make custom colormap for surface plots
vec = [100; 0]; % nodes at which to place the colors
raw = [Mcol; Ccol];
N = 256;
CMmap = interp1(vec, raw, linspace(100,0, N),'pchip');

% region of bistability
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

% set of fishing pressure values
fset = linspace(0.05, 0.15, 100);

% holding vectors for equilibrium values
Cstars = NaN(length(fset), 4); % coral
Mstars = NaN(length(fset), 4); % macroalgae (vuln + invuln)

%Mistars = NaN(length(fset), 4);
%Mvstars = NaN(length(fset), 4);

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

for i = 1:length(fset)%for each fishing pressure
    
    fi = fset(i);

    % get the eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
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
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3); % unstable eq
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still (high and low
% coral eq)


%% Briggs: PDE bifurcation diagram

% PDE parameters
diffs = [0.05,0.05,0.2, 0]; % diffusion rates of MI, C, H, and Mv
taxisM = 0; 
taxisC = -0.5; % taxis rate toward coral
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
fset21 = linspace(0.09, 0.13, 20); 

% holding arrays
Cruns = NaN(length(tset), length(xset),length(fset21));
Mruns = NaN(length(tset), length(xset),length(fset21));
Hruns = NaN(length(tset), length(xset),length(fset21));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(fset21));
Mmeans = NaN(length(fset21));
Hmeans = NaN(length(fset21));

% record the maxes and mins as well
Cmxs = NaN(length(fset21));
Cmns = NaN(length(fset21));
Mmxs = NaN(length(fset21));
Mmns = NaN(length(fset21));

summ10 = 1; % 1 = record peak summaries, 0 = record all peaks
pksumm = NaN(6,3,length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

    tic
    for i = 1:length(fset21) % for each fishing pressure

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(:, :,i) = solij(:,:,2);
    Mruns(:, :, i) = solij(:,:,1)+ solij(:,:,4);
    Hruns(:, :, i) = solij(:,:,3);

    % record spatial averages at final time point
    Cmeans(i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i) = mean(solij(end, b1i:b2i, 3));

    Cmxs(i) = max(solij(end, b1i:b2i, 2));
    Cmns(i) = min(solij(end, b1i:b2i, 2));
    Mmxs(i) = max(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Mmns(i) = min(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));

     % calculate the peak summary metrics
      Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
      [npks0, npks, pklambdas,pkwidths,pkproms,pkheights] = peakfun(Cvalsijk,Mvalsijk,summ10,xset, pkthresh, b1, b2);

            % record these
            pksumm(6,1,i) = npks0/(b2-b1);
            pksumm(5,1,i) = npks/(b2-b1);
            pksumm(1,1:length(pklambdas),i) = pklambdas;
            pksumm(2,1:length(pkwidths),i) = pkwidths;
            pksumm(3,1:length(pkproms),i) = pkproms;
            pksumm(4,1:length(pkheights),i) = pkheights;

   
    end 

    toc % 31 seconds
 

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
load('code output/FigS6.mat','Mmxall','Mmnall', 'Cmxall', 'Cmnall')


%% plot bifurcation diagram for macroalgae (Fig. 1a) and coral (Fig. 1b)

fpts = [2, 6, 10, 14]; % elements of fset21 to highlight in the figure

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
xlim([0.08, 0.14])
hold on
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
text(0.113, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% add the ensemble max and mins
rng = find(Mmxall(:,1)-Mmnall(:,1)> 0.01);
    btwx = [fset21(rng), fliplr(fset21(rng))];
    btwy = [Mmnall(rng,1)', fliplr(Mmxall(rng,1)')];
    fill(btwx, btwy, Mcol, 'FaceAlpha',0.15, 'EdgeColor', Mcol, 'EdgeAlpha', 0.5);
% add circle around f value for which spatial distribution is shown
plotj = fpts(1); % f = 0.0921
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.002, Mmeans(plotj) + 0.03, 'f = 0.092', 'Color', fcol,'FontSize', 14)
plotj = fpts(2); % f = 0.1005
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.008, Mmeans(plotj) + 0.03, 'f = 0.101', 'Color', fcol,'FontSize', 14)
plotj = fpts(3); % f = 0.1089
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)-0.0072, Mmeans(plotj) + 0.03, 'f = 0.109', 'Color', fcol,'FontSize', 14)
plotj = fpts(4); % f = 0.1174
plot(fset21(plotj), Mmeans(plotj), '.','MarkerSize',30,'Color', Mcol)
plot(fset21(plotj), Mmeans(plotj), 'o','MarkerSize',8,'Color', fcol, 'LineWidth',2.5)
text(fset21(plotj)+0.0013, Mmeans(plotj) - 0.025, 'f = 0.117', 'Color', fcol,'FontSize', 14)
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
xlim([0.08, 0.14])
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
text(-90, 1.5,'f = 0.117','FontSize',14, 'Color', fcol)
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
text(-90, 1.5,'f = 0.109','FontSize',14, 'Color', fcol)
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
text(-90, 1.5,'f = 0.101','FontSize',14, 'Color', fcol)
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
text(-90, 1.5,'f = 0.092','FontSize',14, 'Color', fcol)
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

pgon = polyshape([flow flow fup fup],[30 -1 -1 30]);


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
text(0.113, 0.09, 'Bistable', 'Color', 'black','FontSize', 16)
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


















