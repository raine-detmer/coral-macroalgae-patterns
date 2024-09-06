% README: code for making Figure S7


%% set up
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];
Hcol = [0.9294 0.6941 0.1255];

flow = 0.1111; % lower tipping point (calculated in Fig2.m)
fup = 0.1229; % upper tipping point

%% parameter set up

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

len = 400;
xset = linspace(-len/2,len/2,800);

% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85;
C0low = 0.05;
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

ftest = 0.1106;

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


%% step H, homogeneous C and M, high vs low M eq
% this is the same as the simulations with varying step width in ICall.m

icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

icchoiceH = 2;

% start with default taxis and diffusion values
taxisC = -0.5;
diffs = [0.05, 0.05, 0.2, 0];

% values of fishing pressure
fset21 = linspace(0.09, 0.13, 20); 

ftest = fset21(end);

% low and high cover inside and outside patches
lowset = [0.05];
highset = [0.85];

% step widths
wset = round([length(xset), length(xset)/2, length(xset)/4, length(xset)/8, length(xset)/16, length(xset)/32, length(xset)/64, length(xset)/96, length(xset)/128, 0]);

% first just get the initial conditions for the example plot
C0widths = wset(5);
initC = stepfun(C0widths, xset); 

% starting near the high M equilibrium
C0high = 0.1;%highset(1);
C0low = 0.1;%lowset(1);
M0high = 0.4;%highset(1);
M0low = 0.4;%lowset(1);


H0widths = wset(5);
initH = stepfun(H0widths, xset); 

[solij] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, icchoiceH, initH); 

    % record full results
    CrunsHst0 = solij(1,b1i:b2i,2);
    MrunsHst0 = solij(1,b1i:b2i,1)+ solij(1,b1i:b2i,4);
    HrunsHst0 = solij(1,b1i:b2i,3);

    %% check these
%   tcheck = 1;
% 
%  figure(1)
% plot(xset, solij(tcheck,:,2), 'LineWidth',2, 'Color', Ccol)
% ylim([0 1.5])
% xlabel('Location', 'FontSize',18)
% ylabel('Prop. cover', 'FontSize',18)
% hold on 
% plot(xset, solij(tcheck,:,1) + solij(tcheck,:,4), 'LineWidth',2, 'Color', Mcol)
% plot(xset, solij(tcheck,:,3), 'LineWidth',2, 'Color', Hcol)
% hold off
% legend('Coral','Macroalgae', 'Location','Northwest','NumColumns',2)
% lgd = legend;
% lgd.FontSize = 14;

%% full runs, high M eq

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
 %for k = length(wset)

  H0widths = wset(k);
  initH = stepfun(H0widths, xset); 
   
    C0widths = wset(k);
    initC = stepfun(C0widths, xset); 

 for j = 1:length(highset) % for each step height

    
     C0high = 0.1;%highset(1);
     C0low = 0.1;%lowset(1);
     M0high = 0.4;%highset(1);
     M0low = 0.4;%lowset(1);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, icchoiceH, initH); 

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

toc % took 292 seconds


% save these results
CrunsHst = Cruns;
MrunsHst = Mruns;
HrunsHst = Hruns;
CmeansHst = Cmeans;
MmeansHst = Mmeans;
HmeansHst = Hmeans;


%% low M eq

% first just get the initial conditions for the example plot
C0widths = wset(5);
initC = stepfun(C0widths, xset); 

% starting near the high C equilibrium
C0high = 0.8;
C0low = 0.8;
M0high = 0.1;
M0low = 0.1;

H0widths = wset(5);
initH = stepfun(H0widths, xset); 

[solij] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, icchoiceH, initH); 

    % record full results
    
    CrunsHst02 = solij(1,b1i:b2i,2);
    MrunsHst02 = solij(1,b1i:b2i,1)+ solij(1,b1i:b2i,4);
    HrunsHst02 = solij(1,b1i:b2i,3);

%% full runs, low M eq

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
 %for k = length(wset)

  H0widths = wset(k);
  initH = stepfun(H0widths, xset); 
   
    C0widths = wset(k);
    initC = stepfun(C0widths, xset); 

 for j = 1:length(highset) % for each step height

      C0high = 0.8;
      C0low = 0.8;
      M0high = 0.1;
      M0low = 0.1;

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, icchoiceH, initH); 

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

toc % took 289 seconds


% NEW SECTION
% save these results
CrunsHst2 = Cruns;
MrunsHst2 = Mruns;
HrunsHst2 = Hruns;

CmeansHst2 = Cmeans;
MmeansHst2 = Mmeans;
HmeansHst2 = Hmeans;


%% step H and C and M

% first just get the initial conditions for the example plot
C0widths = wset(5);
initC = stepfun(C0widths, xset); 
% starting near the high M equilibrium
C0high = highset(1);
C0low = lowset(1);
M0high = highset(1);
M0low = lowset(1);

H0widths = wset(5);
initH = stepfun(H0widths, xset); 

[solij] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, icchoiceH, initH); 

    % record full results
    CrunsHst03 = solij(1,b1i:b2i,2);
    MrunsHst03 = solij(1,b1i:b2i,1)+ solij(1,b1i:b2i,4);
    HrunsHst03 = solij(1,b1i:b2i,3);

%% full runs, step H and C and M

lowset = [0.05];
highset = [0.85];

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
 %for k = length(wset)

  H0widths = wset(k);
  initH = stepfun(H0widths, xset); 
   
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
    [solij] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, icchoiceH, initH); 

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

toc % took 406 seconds


% save these results
CrunsHst3 = Cruns;
MrunsHst3 = Mruns;
HrunsHst3 = Hruns;
CmeansHst3 = Cmeans;
MmeansHst3 = Mmeans;
HmeansHst3 = Hmeans;

%% plotting
% plot results

fset21 = linspace(0.09, 0.13, 20); 

allcols = parula(length(wset));

pgon = polyshape([flow flow fup fup],[2 -1 -1 2]);

MrunsP = MrunsHst0;
CrunsP = CrunsHst0;
HrunsP = HrunsHst0;

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
ylabel('Initial abundance','FontSize',18)
title('Example initial conditions','FontSize',18)
ylim([0 1.5])
hold on 
plot(xset(b1i:b2i),MrunsP, 'LineWidth',2, 'Color', Mcol)
plot(xset(b1i:b2i),HrunsP, 'LineWidth',2, 'Color', Hcol)
hold off
legend('Coral','Macroalgae', 'Herbivores', 'Location','Northwest','NumColumns',2)
lgd = legend;
lgd.FontSize = 14;
legend boxoff
nexttile
Mmeans = MmeansHst;
wset = round([length(xset), length(xset)/2, length(xset)/4, length(xset)/8, length(xset)/16, length(xset)/32, length(xset)/64, length(xset)/96, length(xset)/128, 0]);
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
title(lgd,{'Initial herbivore ''patch'' width';'(fraction total space)'})
nexttile
% high C eq
Mmeans = MmeansHst2;
MrunsP = MrunsHst02;
CrunsP = CrunsHst02;
HrunsP = HrunsHst02;
plot(xset(b1i:b2i),MrunsP, 'LineWidth',2, 'Color', Mcol)
ylim([0 1.5])
ylabel('Initial abundance','FontSize',18)
hold on 
plot(xset(b1i:b2i),CrunsP, 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP, 'LineWidth',2, 'Color', Hcol)
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
lgd = legend('1','1/2','1/4','1/8','1/16','1/32','1/64','1/96','1/128','0','','', 'Location','northwest');
title(lgd,{'Initial herbivore ''patch'' width';'(fraction total space)'})
nexttile
% all are stepped
Mmeans = MmeansHst3;
MrunsP = MrunsHst03;
CrunsP = CrunsHst03;
HrunsP = HrunsHst03;
plot(xset(b1i:b2i),MrunsP, 'LineWidth',2, 'Color', Mcol)
ylim([0 1.5])
xlabel('Location','FontSize',18)
ylabel('Initial abundance','FontSize',18)
hold on 
plot(xset(b1i:b2i),CrunsP, 'LineWidth',2, 'Color', Ccol)
plot(xset(b1i:b2i),HrunsP, 'LineWidth',2, 'Color', Hcol)
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
xlabel('Fishing pressure','FontSize',18)
ylabel('Equilibrium M cover','FontSize',18)
xline([flow fup]) % bistability region
lgd = legend('1','1/2','1/4','1/8','1/16','1/32','1/64','1/96','1/128','0','','', 'Location','northwest', 'NumColumns', 2);
title(lgd,{'Initial patch width';'(fraction total space)'})

%exportgraphics(figure(1), 'FigS7.pdf', 'ContentType', 'vector');

% saving oversized pdfs:
% https://www.mathworks.com/matlabcentral/answers/1708745-missing-part-of-my-figure-when-saving-as-pdf

%% save results

save('code output/FigS7.mat','CmeansHst', 'MmeansHst', 'HmeansHst', 'CmeansHst2', ...
    'MmeansHst2', 'HmeansHst2', 'CmeansHst3','MmeansHst3', 'HmeansHst3')



%% load results

load('code output/FigS7.mat','CmeansHst', 'MmeansHst', 'HmeansHst', 'CmeansHst2', ...
    'MmeansHst2', 'HmeansHst2', 'CmeansHst3','MmeansHst3', 'HmeansHst3')




