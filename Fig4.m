% README: code for making Figure 4



%% PDE setup
% parameter setup

% PDE parameters
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = 0;%0; % taxis rate toward coral
taxisT = 0;

taxisM1 = 0;
taxisM2 = 0;
taxisC1 = 0;
taxisC2 = 0;
taxisT1 = 0;
taxisT2 = 0;

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges

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
C0widths = round(length(xset)/64);  % step widths
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
fset21 = linspace(0.07, 0.125, 20); % for higher bistability region

rH = 0.2; % herbivore growth rate
fhigh21 = 2*fset21; % if 2f > 0.2, fhigh = 0.2, otherwise fhigh = 2f
fhigh21(fhigh21>rH) = rH;



%% simulation set up

% first choose a fishing pressure just below the lower nonspatial tipping point

finitL = linspace(0.11, 0.115, 100);
finitU = linspace(0.12, 0.125, 100);

parset = [0.01, 0.1, 0.4, 0.4, 0.02, 0.01, 0.5, 0.2, 2, 2, 0.4, 0.2, 0.1];

tic
[lowtp, uptp] =  tpfun(parset, finitL, finitU);
toc % about 33 seconds

% save fishing pressure to use
favg = lowtp(1)-0.005*lowtp(1);
%favg = 0.1105+0.0000053053;

% now calculate equilibrium macroalgal cover at this fishing pressure

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


% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

    fi = favg;

    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
   % Cstars(1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
   % Mstars(1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);

   Mref = soli.Mi+soli.Mv;
   Mref = Mref(2);

%% PDE simulations
% get set of fishing pressures to use
%fhigh21 = 2*favg; % if 2f > 0.2, fhigh = 0.2, otherwise fhigh = 2f
%fhigh21(fhigh21>rH) = rH;
%fhigh21 = 0.2; % max high fishing pressure is 0.2
%2*favg-0.2 % min is 0.0210
f1set = round(linspace(0.021, 0.2, 9),4); % 9 entries and round a little to ensure ratio of 1 is included
f2set = 2*round(favg,4)-f1set;

fset21 = favg;

% 4 different options for the second population
taxisC2set = [1 1 0 0];
diffH2set = [0.2 1 0.2 1];

parset = taxisC2set;

taxisC1 = -1; % first population is strongly attracted to coral

% also record the proportion of herbivores that are the pattern-driving
% population

% holding arrays
Cruns = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));
Mruns = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));
H1runs = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));
H2runs = NaN(length(tset), length(xset),length(fset21), length(parset),length(f1set));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset),length(f1set));
Mmeans = NaN(length(fset21), length(parset),length(f1set));
H1means = NaN(length(fset21), length(parset),length(f1set));
H2means = NaN(length(fset21), length(parset),length(f1set));

Hprops = NaN(length(fset21), length(parset),length(f1set));

tic
for k = 1:length(parset) % for each step width
   
    taxisC2 = taxisC2set(k);
    diffs = [0.05,0.05,0.2, 0,diffH2set(k)];

 for j = 1:length(f1set) % for different fishing pressure ratios

    f1test = f1set(j);
    f2test = f2set(j);

    for i = 1:length(fset21) % for each fishing pressure
    %for i = 1

     % run PDE
    [solij] = Briggs2HrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, f1test, f2test,diffs,taxisM1,taxisC1, taxisT1,taxisM2,taxisC2, taxisT2, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(:, :,i, k, j) = solij(:,:,2);
    Mruns(:, :, i, k, j) = solij(:,:,1)+ solij(:,:,4);
    H1runs(:, :, i, k, j) = solij(:,:,3);
    H2runs(:, :, i, k, j) = solij(:,:,5);

    % record spatial averages at final time point
    % just take the averages from the middle to avoid edge effects
    Cmeans(i, k, j) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k, j) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    H1means(i, k, j) = mean(solij(end, b1i:b2i, 3));
    H2means(i, k, j) = mean(solij(end, b1i:b2i, 5));

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

parset = linspace(0,1,9); % just for getting colors

allcols2 = gray(length(parset(5:end))+1);
allcols2 = flip(allcols2(1:length(parset(5:end)),:));

% colors
Ctx0 = [0.750    0.750    0.750];
Ctx1 = [0 0 0];

Mmeans = Mmeans1;
Hprops = Hprops1;

favg = 0.1105+0.0000053053; % calculated above and pasted here
f1set = round(linspace(0.021, 0.2, 9),4); % 9 entries and round a little to ensure ratio of 1 is included
f2set = 2*round(favg,4)-f1set;

%fratios = f1set./f2set;
% update: plot f proportion instead
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
hold on
plot(fratios,squeeze(Hprops(1,2,:)),'--', 'LineWidth',2.75, 'Color', Ctx1) % 4th element for transparency
plot(fratios,squeeze(Hprops(1,3,:)),'-', 'LineWidth',2.25, 'Color', Ctx0) % 4th element for transparency
plot(fratios,squeeze(Hprops(1,4,:)),'--', 'LineWidth',3.75, 'Color', Ctx0) % 4th element for transparency
hold off
ylabel({'Proportion herbivores'; 'in population A'},'FontSize',22)
text(0.01, 0.97, 'a)', 'Color', 'black','FontSize', 18)
%legend elements
 hold on
lg2{1} = plot(nan, '-', 'LineWidth', 2.5,'Color', Ctx1);%'MarkerSize',7
lg2{2} = plot(nan, '--', 'LineWidth', 2.5,'Color', Ctx1);
lg2{3} = plot(nan, '-', 'LineWidth', 2.5,'Color', Ctx0);
lg2{4} = plot(nan, '--', 'LineWidth', 2.5,'Color', Ctx0);
hold off
legend([lg2{:}],{'\tau_{c} = +1, D_{H} = 0.2', '\tau_{c} = +1, D_{H} = 1', '\tau_{c} = 0, D_{H} = 0.2', '\tau_{c} = 0, D_{H} = 1'}, 'Location', 'northeast')
lgd = legend;
%title(lgd,{'Population two''s';'traits'})
title(lgd,{'Population B';'traits'})
lgd.FontSize = 14;
nexttile

% mean macroalgal cover
%figure(9)
plot(fratios,squeeze(Mmeans(1,1,:)),'-', 'LineWidth',2.25, 'Color', Ctx1) 
yline(double(Mref), 'Linewidth',1.5)
hold on
plot(fratios,squeeze(Mmeans(1,2,:)),'--', 'LineWidth',2.75, 'Color', Ctx1) 
plot(fratios,squeeze(Mmeans(1,3,:)),'-', 'LineWidth',2.25, 'Color', Ctx0) 
plot(fratios,squeeze(Mmeans(1,4,:)),'--', 'LineWidth',3.75, 'Color', Ctx0) 
text(0.02, 0.035,{'Nonspatial'; 'equilibrium'},'FontSize',14)
hold off
xlabel(t,'Proportion fishing effort on population A','FontSize',22)
ylabel('Mean macroalgal cover','FontSize',22)
text(0.01, 0.38, 'b)', 'Color', 'black','FontSize', 18)


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


%% load results
% load('code output/Fig4.mat','Mmeans1', 'Hprops1')

