% README: code for making Figure S6 (and the max/min cover within the Busse
% balloon plotted in Fig. 1a-b)


%% set up
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];

flow = 0.1111; % lower tipping point (calculated in Fig2.m)
fup = 0.1229; % upper tipping point

%% parameter set up

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
rH = 0.2;% herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
f = 0.08; % herbivore fishing pressure

% PDE parameters
diffs = [0.05,0.05,0.2, 0]; % diffusion rates 
taxisM = 0; 
taxisC = -0.5; % taxis rate toward coral
taxisT = 0;

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges

% space
len = 400;
xset = linspace(-len/2,len/2,800);

% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85;
C0low = 0.05;%0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 4
C0widths = round(length(xset)/64);  % step widths
initC = stepfun(C0widths, xset); 

% for icchoice = 5
 ampC0 = (C0high-C0low)/2;
 ampM0 = (M0high-M0low)/2;
 period0 = 0.4;

% for icchoice = 3
rnsize = 1; % magnitude of random variation (0-1)

ftest = 0.1106; % just below tipping point

% time
t_end = 3*50000;% 6000 then 100000 then 500000
tset = linspace(0,t_end,2*2500); % 600 then 10000 then 50000

% peak characteristics
pkthresh = 0.05; % min prominence that a peak has to have to count
dthresh = 0.25*len; % threshold distance from edge before a peak gets considered
b1 = xset(1) + dthresh; % lower boundary for peak consideration
b2 = xset(end)-dthresh; % upper boundary for peak consideration
summ10 = 1; % 1 = record peak summaries, 0 = record all peaks

% get the indeces of these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));

%% step IC: vary step width

icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function


%  values of fishing pressure
fset21 = linspace(0.09, 0.13, 20); 

% low and high cover inside and outside patches
lowset = [0.05];
highset = [0.85];

% patch widths
wset = round([length(xset), length(xset)/2, length(xset)/4, length(xset)/8, length(xset)/16, length(xset)/32, length(xset)/64, length(xset)/96, length(xset)/128, 0]);

% first just do an example simulation to get a plot illustrating the
% initial conditions
C0widths = wset(5);
initC = stepfun(C0widths, xset); 
C0high = highset(1);
C0low = lowset(1);
M0high = highset(1);
M0low = lowset(1);

[solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record initial results
    CrunsST0 = solij(1,b1i:b2i,2);
    MrunsST0 = solij(1,b1i:b2i,1)+ solij(1,b1i:b2i,4);
    HrunsST0 = solij(1,b1i:b2i,3);


%% full runs

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));
Mruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));
Hruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(highset), length(wset));
Mmeans = NaN(length(fset21), length(highset), length(wset));
Hmeans = NaN(length(fset21), length(highset), length(wset));


tic
for k = 1:length(wset) % for each step width
   
    % get the initial conditions
    C0widths = wset(k);
    initC = stepfun(C0widths, xset); 

 for j = 1:length(highset) % for each step height

     C0high = highset(j);
     C0low = lowset(j);
     M0high = highset(j);
     M0low = lowset(j);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record final spatial distributions
    Cruns(1, :,i, j, k) = solij(end,:,2);
    Mruns(1, :, i, j, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, j, k) = solij(end,:,3);

    % record spatial averages at final time point
    Cmeans(i, j, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, j, k) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, j, k) = mean(solij(end, b1i:b2i, 3));

    end 
 end

end 

toc % took 421 seconds

% save these results
CrunsST = Cruns;
MrunsST = Mruns;
HrunsST = Hruns;

CmeansST = Cmeans;
MmeansST = Mmeans;
HmeansST = Hmeans;

%% sin IC: vary patch frequency

icchoice = 5; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

% frequencies
wset = linspace(0.1,1,10);
% period0*xi ; technically f is the frequency for sin(2pi*f*xi), so here
% period0 is the angular frequency, and the frequency would be period0/2pi

% first just get an example of the initial conditions
period0 = wset(2);
C0high = (highset(1) + lowset(1))/2;
M0high = (highset(1) + lowset(1))/2;
ampC0 = (highset(1)-lowset(1))/2;
ampM0 = (highset(1)-lowset(1))/2;

[solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    CrunsS0 = solij(1,b1i:b2i,2);
    MrunsS0 = solij(1,b1i:b2i,1)+ solij(1,b1i:b2i,4);
    HrunsS0 = solij(1,b1i:b2i,3);


%% full simulations
% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));
Mruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));
Hruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(highset), length(wset));
Mmeans = NaN(length(fset21), length(highset), length(wset));
Hmeans = NaN(length(fset21), length(highset), length(wset));


tic
for k = 1:length(wset) % for each step width

    period0 = wset(k);

 for j = 1:length(highset) % for each step height

    C0high = (highset(j) + lowset(j))/2;
    M0high = (highset(j) + lowset(j))/2;
    ampC0 = (highset(j)-lowset(j))/2;
    ampM0 = (highset(j)-lowset(j))/2;

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :,i, j, k) = solij(end,:,2);
    Mruns(1, :, i, j, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, j, k) = solij(end,:,3);

    % record spatial averages at final time point
    Cmeans(i, j, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, j, k) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, j, k) = mean(solij(end, b1i:b2i, 3));

    end 
 end

end 

toc % took 353 seconds

% save these results
CrunsS = Cruns;
MrunsS = Mruns;
HrunsS = Hruns;

CmeansS = Cmeans;
MmeansS = Mmeans;
HmeansS = Hmeans;


%% random

% initial conditions
icchoice = 3; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 5
C0high = 0.5;
C0low = 0.05;%0.05;
M0high = 0.5;
M0low = 0.05;

% get example initial conditions
rng(1) % random 1
[solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    CrunsR0 = solij(1,b1i:b2i,2);
    MrunsR0 = solij(1,b1i:b2i,1)+ solij(1,b1i:b2i,4);
    HrunsR0 = solij(1,b1i:b2i,3);

%% full simulations

% vary max cover and random replicate
highset = [0.5]; % Coral cover
highMset = [0.5]; % M cover

wset = [1, 1000, 10000]; % random replicate

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));
Mruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));
Hruns = NaN(1, length(xset),length(fset21), length(highset), length(wset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(highset), length(wset));
Mmeans = NaN(length(fset21), length(highset), length(wset));
Hmeans = NaN(length(fset21), length(highset), length(wset));


tic
for k = 1:length(wset) % for each replicate


 for j = 1:length(highset) % for each step height

    C0high = highset(j);
    M0high = highMset(j);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
     rng(wset(k)) % random k
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :,i, j, k) = solij(end,:,2);
    Mruns(1, :, i, j, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, j, k) = solij(end,:,3);

    % record spatial averages at final time point
    Cmeans(i, j, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, j, k) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, j, k) = mean(solij(end, b1i:b2i, 3));

    
    end 
 end

end 

toc % took 98 seconds for 3 replicates


% save these results
CrunsR = Cruns;
MrunsR = Mruns;
HrunsR = Hruns;

CmeansR = Cmeans;
MmeansR = Mmeans;
HmeansR = Hmeans;


%% plot results (Fig. S6)
% plot results

% note for some reason the code to plot all the panels together only works
% every other time it is run

fset21 = linspace(0.09, 0.13, 20); 

%allcols = parula(length(wset));

pgon = polyshape([flow flow fup fup],[2 -1 -1 2]); % bistability region

MrunsP = MrunsST0;
CrunsP = CrunsST0;

figure(1)
x0=10;
y0=10;
width=800;
height=900;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(3, 2); % (rows, columns)
t.TileSpacing = 'compact';
t.TileIndexing = 'rowmajor'; % default is rowmajor
nexttile
plot(xset(b1i:b2i),CrunsP, 'LineWidth',2, 'Color', Ccol)
ylabel('Initial cover','FontSize',18)
title('Example initial conditions','FontSize',18)
ylim([0 1])
hold on 
plot(xset(b1i:b2i),MrunsP, 'LineWidth',2, 'Color', Mcol)
hold off
legend('Coral','Macroalgae', 'Location','Northwest','NumColumns',2)
lgd = legend;
lgd.FontSize = 14;
legend boxoff
nexttile
Mmeans = MmeansST;
wset = round([length(xset), length(xset)/2, length(xset)/4, length(xset)/8, length(xset)/16, length(xset)/32, length(xset)/64, length(xset)/96, length(xset)/128, 0]);
allcols = parula(length(wset));
for k = 1:length(wset)
    plot(fset21,Mmeans(:,1,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) 
    hold on;
end
hold off
hold on
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
title('Effect of initial conditions on equilibria','FontSize',18)
hold off
ylabel('Equilibrium M cover','FontSize',18)
text(0.1125, 0.65, 'Bistable', 'Color', 'black','FontSize', 14)
xline([flow fup]) % bistability region
lgd = legend('1','1/2','1/4','1/8','1/16','1/32','1/64','1/96','1/128','0','','', 'Location','northwest', 'NumColumns', 2);
title(lgd,{'Initial coral patch width';'(fraction total space)'})
nexttile
% sine
wset = linspace(0.1,1,10);
Mmeans = MmeansS;
MrunsP = MrunsS0;
CrunsP = CrunsS0;
plot(xset(b1i:b2i),MrunsP, 'LineWidth',2, 'Color', Mcol)
ylim([0 1])
ylabel('Initial cover','FontSize',18)
hold on 
plot(xset(b1i:b2i),CrunsP, 'LineWidth',2, 'Color', Ccol)
hold off
nexttile
for k = 1:length(wset)
    plot(fset21,Mmeans(:,1,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) 
    hold on;
end
hold off
hold on
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
ylabel('Equilibrium M cover','FontSize',18)
xline([flow fup]) % bistability region
lgd = legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1','','', 'Location','northwest', 'NumColumns', 2);
title(lgd,{'Initial coral patch';'frequency'})
nexttile
wset = [1, 1000, 10000]; % random replicate
Mmeans = MmeansR;
MrunsP = MrunsR0;
CrunsP = CrunsR0;
plot(xset(b1i:b2i),MrunsP, 'LineWidth',2, 'Color', Mcol)
ylim([0 1])
xlabel('Location','FontSize',18)
ylabel('Initial cover','FontSize',18)
hold on 
plot(xset(b1i:b2i),CrunsP, 'LineWidth',2, 'Color', Ccol)
hold off
nexttile
for k = 1:length(wset)
    plot(fset21,Mmeans(:,1,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k+5,:)) 
    hold on;
end
hold off
hold on
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
xlabel('Fishing pressure','FontSize',18)
ylabel('Equilibrium M cover','FontSize',18)
xline([flow fup]) % bistability region
lgd = legend('1','2','3','','', 'Location','northwest');
title(lgd,{'Random replicate'})


%% max/min envelope for Fig. 1

% get the maxes at each fishing pressure for each initial conditions
MmxST = NaN(20, 10);
MmnST = NaN(20, 10);
CmxST = NaN(20, 10);
CmnST = NaN(20, 10);
MmxS = NaN(20, 10);
MmnS = NaN(20, 10);
CmxS = NaN(20, 10);
CmnS = NaN(20, 10);

for i = 1:20
    for j = 1:10
     MmxST(i, j) = max(MrunsST(1,b1i:b2i,i,1,j));
     MmnST(i, j) = min(MrunsST(1,b1i:b2i,i,1,j));
     CmxST(i, j) = max(CrunsST(1,b1i:b2i,i,1,j));
     CmnST(i, j) = min(CrunsST(1,b1i:b2i,i,1,j));
     MmxS(i, j) = max(MrunsS(1,b1i:b2i,i,1,j));
     MmnS(i, j) = min(MrunsS(1,b1i:b2i,i,1,j));
     CmxS(i, j) = max(CrunsS(1,b1i:b2i,i,1,j));
     CmnS(i, j) = min(CrunsS(1,b1i:b2i,i,1,j));
    end

end

% get rid of the 1st and 10th entries for the ST runs (these were
% homogeneous initial conditions with no patterns)
MmxST = MmxST(:, 2:9);
MmnST = MmnST(:, 2:9);
CmxST = CmxST(:, 2:9);
CmnST = CmnST(:, 2:9);


% now the overall maxes and mins across all initial conditions
Mmxall = NaN(20, 1);
Mmnall = NaN(20, 1);
Cmxall = NaN(20, 1);
Cmnall = NaN(20, 1);

for i = 1:20
    
Mmxall(i,1) = max([MmxST(i, :),MmxS(i, :)]);
Mmnall(i,1) = min([MmnST(i, :),MmnS(i, :)]);
Cmxall(i,1) = max([CmxST(i, :),CmxS(i, :)]);
Cmnall(i,1) = min([CmnST(i, :),CmnS(i, :)]);

end



%% save results

 save('code output/FigS6.mat','CmeansST', 'MmeansST', 'HmeansST', 'CmeansS', ...
     'MmeansS', 'HmeansS','CmeansR', 'MmeansR', 'HmeansR', 'Mmxall', ...
     'Mmnall', 'Cmxall', 'Cmnall')


%% load results

load('code output/FigS6.mat','CmeansST', 'MmeansST', 'HmeansST', 'CmeansS', ...
     'MmeansS', 'HmeansS','CmeansR', 'MmeansR', 'HmeansR', 'Mmxall', ...
     'Mmnall', 'Cmxall', 'Cmnall')


