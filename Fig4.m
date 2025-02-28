% README: code for making Figure 4

% takes several minutes to run, or can load stored output

load('code output/Fig4.mat','Mmeans1', 'Hprops1')


%% PDE setup
% parameter setup

% PDE parameters
diffs = [0.05,0.05,0.25, 0]; % diffusion rates 
taxisM = 0; % taxis rate toward macroalgae
taxisC = 0; % taxis rate toward coral
taxisT = 0; % taxis rate toward turf

taxisM1 = 0; % first herbivore population
taxisM2 = 0; % second herbivore population
taxisC1 = 0; % first herbivore population
taxisC2 = 0; % second herbivore population
taxisT1 = 0; % first herbivore population
taxisT2 = 0; % second herbivore population

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges

% space parameters
len = 400;
xset = linspace(-len/2,len/2,800);

% time parameters
t_end = 3*50000;
tset = linspace(0,t_end,2*2500); 


% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85;
C0low = 0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 3 (random)
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4 (step)
C0widths = round(length(xset)/64);  % step widths
initC = stepfun(C0widths, xset); 

% for icchoice = 5 (sinusoidal)
 ampC0 = (C0high-C0low)/2;
 ampM0 = (M0high-M0low)/2;
 period0 = 0.4;

% peak characteristics
pkthresh = 0.05; % min prominence that a peak has to have to count
dthresh = 0.25*len; % threshold distance from edge before a peak gets considered
b1 = xset(1) + dthresh; % lower boundary for peak consideration
b2 = xset(end)-dthresh; % upper boundary for peak consideration

% get the indeces of xset corresponding to these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));


%% simulation set up


% region of bistability
flow = 0.1674; % lower tipping point (calculated in Fig2.m)
fup = 0.1878; % upper tipping point

favg = 0.99*flow; % choose a fishing pressure just below the lower nonspatial tipping point

% now calculate equilibrium macroalgal cover at this fishing pressure

% state variables
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
phiH = 0.05; % recruitment rate


% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

    fi = favg; % fishing pressure

    % equations
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    
   % store the value of equilibrium macroalgal cover
   Mref = soli.Mi+soli.Mv;


%% check herbivore equilibria

% ((rH-favg) + sqrt((rH-favg)^2 + 4*dH*phiH))/(2*dH)
% 
% ((rH-favg) + sqrt((rH-favg)^2 + 4*2*dH*0.5*phiH))/(2*2*dH) + ((rH-favg) + sqrt((rH-favg)^2 + 4*2*dH*0.5*phiH))/(2*2*dH)
% 
% % choose fishing pressures to keep HA + HB constant
% % hhh = ((rH-favg) + sqrt((rH-favg)^2 + 4*2*dH*0.5*phiH))/(2*2*dH);
% % rH-(4*dH*hhh^2-phiH)/(2*hhh)
% % favg
% 
% Htot = ((rH-favg) + sqrt((rH-favg)^2 + 4*dH*phiH))/(2*dH);
% 
% % fishing pressures on population A
% f1set = linspace(0.3*favg, 3*favg, 10);
% HAeq = ((rH-f1set) + sqrt((rH-f1set).^2 + 4*2*dH*0.5*phiH))/(2*2*dH); % need .^ to operate on each element
% 
% HBeq = Htot - HAeq;
% 
% f2set = rH-(4*dH*HBeq.^2-phiH)./(2*HBeq)

%% PDE simulations

% total herbivore abundance in all cases
Htot = ((rH-favg) + sqrt((rH-favg)^2 + 4*dH*phiH))/(2*dH);
% set of fishing pressures on population A
f1set = linspace(0.3*favg, 3*favg, 10);
% corresponding biomass of population A
HAeq = ((rH-f1set) + sqrt((rH-f1set).^2 + 4*2*dH*0.5*phiH))/(2*2*dH); % need .^ to operate on each element
% biomass of population B
HBeq = Htot - HAeq;
% fishing pressure on population B that produces biomass equal to HBeq
f2set = rH-(4*dH*HBeq.^2-phiH)./(2*HBeq);

fset21 = favg; % fishing pressure on single population

% 4 different options for the second population
taxisC2set = [1.25 1.25 0 0]; % taxis towards coral for population B
diffH2set = [0.25 1 0.25 1]; % diffusion rates of population B

parset = taxisC2set; % iterate over elements in taxisC2set

taxisC1 = -1.25; % first population is strongly attracted to coral

% holding arrays
% full model output
Cruns = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));
Mruns = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));
H1runs = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));
H2runs = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));

% avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset),length(f1set));
Mmeans = NaN(length(fset21), length(parset),length(f1set));
H1means = NaN(length(fset21), length(parset),length(f1set));
H2means = NaN(length(fset21), length(parset),length(f1set));

% proportion of herbivores that are the pattern-driving population
Hprops = NaN(length(fset21), length(parset),length(f1set));

tic
for k = 1:length(parset) % for each element in parset
   
    taxisC2 = taxisC2set(k); % set the taxis rate of the second herbivore population
    diffs = [0.05,0.05,0.25, 0,diffH2set(k)]; % set the diffusion rate of the second herbivore population

 for j = 1:length(f1set) % for each element of f1set

    f1test = f1set(j); % fishing pressure on first population
    f2test = f2set(j); % fishing pressure on second population

    for i = 1:length(fset21) % for each average fishing pressure
    %for i = 1

     % run PDE
    [solij] = Briggs2HrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, f1test, f2test,diffs,taxisM1,taxisC1, taxisT1,taxisM2,taxisC2, taxisT2, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record full results
    Cruns(:, :,i, k, j) = solij(:,:,2);
    Mruns(:, :, i, k, j) = solij(:,:,1)+ solij(:,:,4);
    H1runs(:, :, i, k, j) = solij(:,:,3);
    H2runs(:, :, i, k, j) = solij(:,:,5);

    % record spatial averages at final time point
    Cmeans(i, k, j) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k, j) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    H1means(i, k, j) = mean(solij(end, b1i:b2i, 3));
    H2means(i, k, j) = mean(solij(end, b1i:b2i, 5));

    % record proportion of total herbivore abundance in population 1
    Hprops(i, k, j) = H1means(i, k, j)/(H1means(i, k, j)+H2means(i, k, j));


    end 
 end

end 

toc % 132 seconds


% save results

Cruns1 = Cruns;
Mruns1 = Mruns;
H1runs1 = H1runs;
H2runs1 = H2runs;

Cmeans1 = Cmeans;
Mmeans1 = Mmeans;
H1means1 = H1means;
H2means1 = H2means;

Hprops1 = Hprops;

%% plot results


% Mref = 0.01074925336433911319721804609953; % calculated above and pasted here


% colors
Ctx0 = [0.750    0.750    0.750];
Ctx1 = [0 0 0];

Mmeans = Mmeans1;
Hprops = Hprops1;

flow = 0.1674; % lower tipping point (calculated in Fig2.m)
fup = 0.1878; % upper tipping point
favg = 0.99*flow; % fishing pressure on single herbivore population

% re-define parameter sets here so this code chunk will run if just using
% the loaded results from Fig4.mat (will also have to run the first two
% code chunks)
% total herbivore abundance in single population
Htot = ((rH-favg) + sqrt((rH-favg)^2 + 4*dH*phiH))/(2*dH);
% set of fishing pressures on first population
f1set = linspace(0.3*favg, 3*favg, 10);
% corresponding biomass of first population
HAeq = ((rH-f1set) + sqrt((rH-f1set).^2 + 4*2*dH*0.5*phiH))/(2*2*dH); % need .^ to operate on each element
% biomass of second population
HBeq = Htot - HAeq;
% corresponding fishing pressure on second population
f2set = rH-(4*dH*HBeq.^2-phiH)./(2*HBeq);

% proportion of fishing pressure on first population
fratios = f1set./(f1set + f2set);

% H1 proportions
figure(1)
x0=10;
y0=10;
width=600;
height=850;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 1);
t.TileIndexing = 'rowmajor'; % default is rowmajor
t.TileSpacing = 'tight';
nexttile
plot(fratios,squeeze(Hprops(1,1,:)),'-', 'LineWidth',2.25, 'Color', Ctx1) % 'MarkerSize',7
ylim([0 1.1])
hold on
plot(fratios,squeeze(Hprops(1,2,:)),'--', 'LineWidth',2.75, 'Color', Ctx1) % 4th element for transparency
plot(fratios,squeeze(Hprops(1,3,:)),'-', 'LineWidth',2.25, 'Color', Ctx0) % 4th element for transparency
plot(fratios,squeeze(Hprops(1,4,:)),'--', 'LineWidth',3.75, 'Color', Ctx0) % 4th element for transparency
hold off
ylabel({'Proportion herbivores'; 'in population A'},'FontSize',22)
text(0.01, 1.05, 'a)', 'Color', 'black','FontSize', 18)
%legend elements
 hold on
lg2{1} = plot(nan, '-', 'LineWidth', 2.5,'Color', Ctx1);%'MarkerSize',7
lg2{2} = plot(nan, '--', 'LineWidth', 2.5,'Color', Ctx1);
lg2{3} = plot(nan, '-', 'LineWidth', 2.5,'Color', Ctx0);
lg2{4} = plot(nan, '--', 'LineWidth', 2.5,'Color', Ctx0);
hold off
legend([lg2{:}],{'\tau_{c} = +1.25, D_{H} = 0.25', '\tau_{c} = +1.25, D_{H} = 1', '\tau_{c} = 0, D_{H} = 0.25', '\tau_{c} = 0, D_{H} = 1'}, 'Location', 'northeast')
lgd = legend;
%title(lgd,{'Population two''s';'traits'})
title(lgd,{'Population B';'traits'})
lgd.FontSize = 14;
nexttile

% mean macroalgal cover
%figure(9)
plot(fratios,squeeze(Mmeans(1,1,:)),'-', 'LineWidth',2.25, 'Color', Ctx1) 
yline(double(Mref), 'Linewidth',1.5)
ylim([0 0.42])
hold on
plot(fratios,squeeze(Mmeans(1,2,:)),'--', 'LineWidth',2.75, 'Color', Ctx1) 
plot(fratios,squeeze(Mmeans(1,3,:)),'-', 'LineWidth',2.25, 'Color', Ctx0) 
plot(fratios,squeeze(Mmeans(1,4,:)),'--', 'LineWidth',3.75, 'Color', Ctx0) 
text(0.02, 0.035,{'Nonspatial'; 'equilibrium'},'FontSize',14)
hold off
xlabel(t,'Proportion fishing effort on population A','FontSize',22)
ylabel('Mean macroalgal cover','FontSize',22)
text(0.01, 0.4, 'b)', 'Color', 'black','FontSize', 18)


%% check distributions

% plotj = 5; % equal fishing pressure
% plotk = 2; % second group = avoids coral and has high diffusion
% 
% %plotj = 1; % lowest fishing pressure on focal group
% %plotk = 2; % second group = avoids coral and has high diffusion
% 
% figure(2)
% plot(xset(b1i:b2i),Cruns5(end, b1i:b2i, 1, plotk,plotj), 'LineWidth',2, 'Color', Ccol)
% ylim([-0.01 1.5])
% hold on
% plot(xset(b1i:b2i),Mruns5(end, b1i:b2i, 1, plotk,plotj), 'LineWidth',2, 'Color', Mcol)
% plot(xset(b1i:b2i),H1runs5(end, b1i:b2i, 1, plotk,plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
% plot(xset(b1i:b2i),H2runs5(end, b1i:b2i, 1, plotk,plotj), 'LineWidth',2, 'Color', [0.7176    0.2745    1.0000])
% hold off
% legend('Coral cover', 'Macroalgal cover', 'H1 biomass','H2 biomass', 'location', 'northeast', 'FontSize',14);

%% save results

save('code output/Fig4.mat','Mmeans1', 'Hprops1')



