% README: code for making Figures S10, S13, and S14

% sensitivity of all the models to whether or not there is external
% recruitment

%% setup
% plotting colors
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];

% external recruitment sets
extCs = [0.01, 0.01, 0, 0];
extMs = [0.01, 0, 0.01, 0];

%% Briggs model ODE

% define the symbols
syms Mi C H Mv

% define parameters
gTC = 0.1; %0.1
gamma = 0.4; %0.4
gTI = 0.4;%0.4
rM = 0.5; %0.5
gTV = 0.2;%0.2
dv = 2; % 2, grazing rate on vulnerable M
omega = 2; % 2
di = 0.4; % 0.4, grazing rate on invulnerable M
dC = 0.02;% default 0.05 

%phiM = 0.01; % default 0.0001
%phiC = 0.01; % default 0.001

% herbivore parameters
rH = 0.2;%0.1; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
f = 0.08; % herbivore fishing pressure

% set of fishing values
fset = linspace(0.05, 0.145, 100);


% holding arrays for each combination of external recruitment pars
Cups = NaN(length(fset), 4); % 4 recruitment par combinations
Cmids = NaN(length(fset), 4);
Clows = NaN(length(fset), 4);

Mups = NaN(length(fset), 4);
Mmids = NaN(length(fset), 4);
Mlows = NaN(length(fset), 4);

bend = NaN(1,4);
bstart = NaN(1,4);

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')


tic
for j = 1:length(extCs)
%for j = 3:4

    phiC = extCs(j);
    phiM = extMs(j);

% holding vector of eq values
Cstars = NaN(length(fset), 8);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 8);

%Mistars = NaN(length(fset), 4);
%Mvstars = NaN(length(fset), 4);

for i = 1:length(fset)%for each element of gset
    % get the eqns
    fi = fset(i);

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

if j==1
% need to rearrange the eq to get smooth lines when plotting
bend(j) = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% use vertcat to concatenate vertical vectors
% look at the Cstars to figure out how to piece these together
% for f on x axis:
Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = Cstars(:, 3);
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 3), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 3));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 3), Mstars(bend(j)+1:end, 1));
% note ups and lows are from the coral's perspective still

elseif j==2
bend(j) = find(isnan(Cstars(:, 5))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 5))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 5)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):bend(j), 3), Cstars(bend(j)+1:end, 5));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 6)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 3), Mstars(bend(j)+1:end, 5));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 4), Mstars(bend(j)+1:end, 3));


elseif j==3

bend(j) = find(isnan(Cstars(:, 4))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 4))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 3), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 3));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 5));
Mlows(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 3), Mstars(bend(j)+1:end, 1));


else

bend(j) = find(isnan(Cstars(:, 7))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 7))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bend(j), 6), Cstars(bend(j)+1:end, 7)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 7), Cstars(bstart(j):bend(j), 5), Cstars(bend(j)+1:end, 7));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 7), Cstars(bstart(j):end, 4));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 7)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 7), Mstars(bstart(j):bend(j), 5), Mstars(bend(j)+1:end, 7));
Mlows(:,j) = vertcat(Mstars(1:bstart(j)-1, 7), Mstars(bstart(j):bend(j), 6), Mstars(bend(j)+1:end, 5));


end

end

toc % about 70 seconds

%% check results
%bstart
% j = 4;
% 
% plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
% ylim([0 1])
% hold on
% plot(fset, Mups(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% plot(fset, Mmids(:,j),'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% hold off

%% Briggs PDE

% PDE parameters
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5;%0; % taxis rate toward coral
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

C0high = 0.85;
C0low = 0.05;%0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 3
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4
%C0widths = round(length(xset)/16);  % step widths
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
summ10 = 1; % 1 = record peak summaries, 0 = record all peaks

% get the indeces of these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));

%  values of fishing pressure
fset21 = linspace(0.07, 0.142, 20);


% holding arrays
Cruns = NaN(length(tset), length(xset),length(extCs),length(fset21));
Mruns = NaN(length(tset), length(xset),length(extCs),length(fset21));
Hruns = NaN(length(tset), length(xset),length(extCs),length(fset21));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(extCs),length(fset21));
Mmeans = NaN(length(extCs),length(fset21));
Hmeans = NaN(length(extCs),length(fset21));

%pksumm = NaN(6,3,length(extCs),length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

tic
for j = 1:length(extCs) % for each recruitment scenario
    % get the recruitment parameters
    phiC = extCs(j);
    phiM = extMs(j);

    for i = 1:length(fset21) % for each fishing pressure
    %for i = 1:5

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(:, :,j, i) = solij(:,:,2);
    Mruns(:, :, j, i) = solij(:,:,1)+ solij(:,:,4);
    Hruns(:, :, j, i) = solij(:,:,3);

    % record spatial averages at final time point
    Cmeans(j, i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(j, i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(j, i) = mean(solij(end, b1i:b2i, 3));

    
   
    end 

end 

toc % 144 seconds

%% check patterns

% figure(2)
% plot(xset, Cruns(end,:,4, 8), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
% hold on 
% plot(xset, Mruns(end,:,4, 8), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
% hold off

%% save all the Briggs results
CupsB = Cups;
CmidsB = Cmids;
ClowsB = Clows;
MupsB = Mups;
MmidsB = Mmids;
MlowsB = Mlows;

fsetB = fset;
bstartB = bstart;
bendB = bend;

fset21B = fset21;

CmeansB = Cmeans;
MmeansB = Mmeans;
HmeansB = Hmeans;
   
CrunsB = Cruns;
MrunsB = Mruns;
HrunsB = Hruns;


%% plot results

Cups = CupsB;
Cmids = CmidsB;
Clows = ClowsB;
Mups = MupsB;
Mmids = MmidsB;
Mlows = MlowsB;

fset = fsetB;
bstart = bstartB;
bend = bendB;

fset21 = fset21B;

Cmeans = CmeansB;
Mmeans = MmeansB;
   

figure(1)
x0=10;
y0=10;
width=1500;
height=700;%375;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 4);
t.TileSpacing = 'compact';
nexttile
% start with coral
Pcol = Ccol;
Pmeans = Cmeans;
Pups = Cups;
Pmids = Cmids;
Plows = Clows;
% C and M recruitment
j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'coral cover'},'FontSize',22)
title({'a) External recruitment of'; 'coral and macroalgae'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
legend('','','','ODE, stable', 'ODE, unstable','', 'PDE, spatial means', 'Location', 'southwest')
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
nexttile
j = 2;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'b) External recruitment of'; 'coral only'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 3;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'c) External recruitment of'; 'macroalgae only'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 4;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title('d) No external recruitment', 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
% now the macroalgae
Pcol = Mcol;
Pmeans = Mmeans;
Pups = Mups;
Pmids = Mmids;
Plows = Mlows;

j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'macroalgal cover'},'FontSize',22)
%title('C and M', 'FontSize',22)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
%text(0.113, 0.65, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 2;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 3;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
ylim([0 1])
xlim([fset(1), fset(end)])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 4;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
ylim([0 1])
xlim([fset(1), fset(end)])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off


%% Mumby model ODE

% define the symbols
syms M C H 

% define parameters
r = 1; % growth rate of corals over turf
a = 0.1; % overgrowth rate of corals by macroalgae
gamma = 0.8; % overgrowth rate of turf by macroalgae
d = 0.44; % coral death rate
gz = 1; % grazing rate

% herbivore growth rate and den dep mortality: carrying capacity = rH/dH
rH = 0.2;%0.1; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
%f = 0; % herbivore fishing pressure

% set of fishing values
fset = linspace(0.14, 0.19, 100);
%fset = linspace(0.14, 0.19, 80);


% holding arrays for each combination of external recruitment pars
Cups = NaN(length(fset), 4); % 4 recruitment par combinations
Cmids = NaN(length(fset), 4);
Clows = NaN(length(fset), 4);

Mups = NaN(length(fset), 4);
Mmids = NaN(length(fset), 4);
Mlows = NaN(length(fset), 4);

bend = NaN(1,4);
bstart = NaN(1,4);

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')


tic
for j = 1:length(extCs)
%for j = 1:2

    beta = extCs(j);
    alpha = extMs(j);

% holding vector of eq values
Cstars = NaN(length(fset), 8);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 8);

%Mistars = NaN(length(fset), 4);
%Mvstars = NaN(length(fset), 4);

for i = 1:length(fset)%for each element of gset
    % get the eqns
    fi = fset(i);

    eq1i = a*M*C-gz*H*M/(1-C)+gamma*M*(1-M-C)+alpha*(1-M-C) == 0;
    eq2i = r*(1-M-C)*C-d*C-a*M*C+beta*(1-M-C) ==0;
    eq3i = rH*H-dH*H*H-fi*H ==0;
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i],[M,C, H], [0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C);
    Mstars(i, 1:length(soli.M)) = sort(soli.M);
end

% process results

if j==1
% need to rearrange the eq to get smooth lines when plotting
bend(j) = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% use vertcat to concatenate vertical vectors
% look at the Cstars to figure out how to piece these together
% for f on x axis:
Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = Cstars(:, 3);
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 3), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 3));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 3), Mstars(bend(j)+1:end, 1));
% note ups and lows are from the coral's perspective still

elseif j==2
bend(j) = find(isnan(Cstars(:, 5))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 5))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 5)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):bend(j), 3), Cstars(bend(j)+1:end, 5));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 6)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 3), Mstars(bend(j)+1:end, 5));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 4), Mstars(bend(j)+1:end, 3));


elseif j==3

bend(j) = find(isnan(Cstars(:, 4))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 4))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 3), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 3));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 5));
Mlows(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 3), Mstars(bend(j)+1:end, 1));


else

bend(j) = find(isnan(Cstars(:, 7))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 7))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bend(j), 6), Cstars(bend(j)+1:end, 7)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 7), Cstars(bstart(j):bend(j), 5), Cstars(bend(j)+1:end, 7));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 7), Cstars(bstart(j):end, 4));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 7)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 7), Mstars(bstart(j):bend(j), 5), Mstars(bend(j)+1:end, 7));
Mlows(:,j) = vertcat(Mstars(1:bstart(j)-1, 7), Mstars(bstart(j):bend(j), 6), Mstars(bend(j)+1:end, 5));


end

end

toc % about 45 seconds

%% check results
%bstart

% j = 4;
% 
% plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
% ylim([0 1])
% hold on
% plot(fset, Mups(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% plot(fset, Mmids(:,j),'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% %plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% %plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% hold off

%% check C

% j = 4;
% 
% plot(fset, Clows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
% ylim([0 1])
% hold on
% plot(fset, Cups(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% plot(fset, Cmids(:,j),'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% %plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% %plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% hold off

%% Mumby PDE

fset21 = linspace(0.145, 0.185, 20);

% holding arrays
Cruns = NaN(length(tset), length(xset),length(extCs),length(fset21));
Mruns = NaN(length(tset), length(xset),length(extCs),length(fset21));
Hruns = NaN(length(tset), length(xset),length(extCs),length(fset21));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(extCs),length(fset21));
Mmeans = NaN(length(extCs),length(fset21));
Hmeans = NaN(length(extCs),length(fset21));

%pksumm = NaN(6,3,length(extCs),length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

tic
for j = 1:length(extCs) % for each recruitment scenario
    % get the recruitment parameters
    beta = extCs(j);
    alpha = extMs(j);

    for i = 1:length(fset21) % for each fishing pressure
    %for i = 1:5

    ftest = fset21(i);

     % run PDE
     [solij] = MumbyHPDE(r,a,gamma,gz, rH, dH, ftest, d, alpha,beta,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(:, :,j, i) = solij(:,:,2);
    Mruns(:, :, j, i) = solij(:,:,1);
    Hruns(:, :, j, i) = solij(:,:,3);

    % record spatial averages at final time point
    Cmeans(j, i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(j, i) = mean(solij(end, b1i:b2i, 1));
    Hmeans(j, i) = mean(solij(end, b1i:b2i, 3));

    
   
    end 

end 

toc % about 160 seconds

%% check results
% fplot = 5;
% jplot = 2;
% 
% figure(3)
% plot(xset, Cruns(end,:,jplot,fplot), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
% ylim([0 1])
% xlabel('Location','FontSize',22)
% ylabel('Prop. cover','FontSize',22)
% hold on 
% plot(xset, Mruns(end,:,jplot,fplot), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
% hold off
% %legend('Coral','Macroalgae', 'Location','Northwest','NumColumns',2)
% %lgd = legend;
% %lgd.FontSize = 14;

%% save all the Mumby results
CupsM = Cups;
CmidsM = Cmids;
ClowsM = Clows;
MupsM = Mups;
MmidsM = Mmids;
MlowsM = Mlows;

fsetM = fset;
bstartM = bstart;
bendM = bend;

fset21M = fset21;

CmeansM = Cmeans;
MmeansM = Mmeans;
HmeansM = Hmeans;
   
CrunsM = Cruns;
MrunsM = Mruns;
HrunsM = Hruns;

%% plot results

Cups = CupsM;
Cmids = CmidsM;
Clows = ClowsM;
Mups = MupsM;
Mmids = MmidsM;
Mlows = MlowsM;

fset = fsetM;
bstart = bstartM;
bend = bendM;

fset21 = fset21M;

Cmeans = CmeansM;
Mmeans = MmeansM;

figure(2)
x0=10;
y0=10;
width=1500;
height=700;%375;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 4);
t.TileSpacing = 'compact';
nexttile
% start with coral
Pcol = Ccol;
Pmeans = Cmeans;
Pups = Cups;
Pmids = Cmids;
Plows = Clows;
% C and M recruitment
j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'coral cover'},'FontSize',22)
title({'a) External recruitment of'; 'coral and macroalgae'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
%text(0.113, 0.65, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
legend('','','','ODE, stable', 'ODE, unstable','', 'PDE, spatial means', 'Location', 'northwest')
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
nexttile
j = 2;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'b) External recruitment of'; 'coral only'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 3;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'c) External recruitment of'; 'macroalgae only'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 4;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'d) No external recruitment'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
% now the macroalgae
Pcol = Mcol;
Pmeans = Mmeans;
Pups = Mups;
Pmids = Mmids;
Plows = Mlows;

j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'macroalgae cover'},'FontSize',22)
%title('C and M', 'FontSize',22)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
%text(0.113, 0.65, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 2;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 3;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
ylim([0 1])
xlim([fset(1), fset(end)])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 4;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
ylim([0 1])
xlim([fset(1), fset(end)])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off


%% van de Leemput model ODE

% define the symbols
syms M C H 

% define parameters
r = 1; % growth rate of corals over turf
a = 0.1; % overgrowth rate of corals by macroalgae
gamma = 0.8; % overgrowth rate of turf by macroalgae
d = 0.44; % coral death rate
gz = 1; % grazing rate
h = 2; % half-saturation constant for herbivore grazing


% herbivore growth rate and den dep mortality: carrying capacity = rH/dH
rH = 0.2;%0.1; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
%f = 0; % herbivore fishing pressure

% set of fishing values
fset = linspace(0.095, 0.17, 100);


% holding arrays for each combination of external recruitment pars
Cups = NaN(length(fset), 4); % 4 recruitment par combinations
Cmids = NaN(length(fset), 4);
Clows = NaN(length(fset), 4);

Mups = NaN(length(fset), 4);
Mmids = NaN(length(fset), 4);
Mlows = NaN(length(fset), 4);

bend = NaN(1,4);
bstart = NaN(1,4);


% turn off warning
warning('off','symbolic:numeric:NumericalInstability')


tic
for j = 1:length(extCs)
%for j = 3:4

    beta = extCs(j);
    alpha = extMs(j);

% holding vector of eq values
Cstars = NaN(length(fset), 8);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 8);

%Mistars = NaN(length(fset), 4);
%Mvstars = NaN(length(fset), 4);

for i = 1:length(fset)%for each element of gset
    % get the eqns
    fi = fset(i);

    eq1i = a*M*C-gz*H*M/(gz*h*M + 1)+gamma*M*(1-M-C)+alpha*(1-M-C) == 0;
    eq2i = r*(1-M-C)*C-d*C-a*M*C+beta*(1-M-C) ==0;
    eq3i = rH*H-dH*H*H-fi*H ==0;
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i],[M,C, H], [0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C);
    Mstars(i, 1:length(soli.M)) = sort(soli.M);
end

% process results

if j==1
% need to rearrange the eq to get smooth lines when plotting
bend(j) = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% use vertcat to concatenate vertical vectors
% look at the Cstars to figure out how to piece these together
% for f on x axis:
Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = Cstars(:, 3);
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 3), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 3));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 3), Mstars(bend(j)+1:end, 1));
% note ups and lows are from the coral's perspective still

elseif j==2
bend(j) = find(isnan(Cstars(:, 5))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 5))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 5)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):bend(j), 3), Cstars(bend(j)+1:end, 5));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 6)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 3), Mstars(bend(j)+1:end, 5));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 4), Mstars(bend(j)+1:end, 3));


elseif j==3

bend(j) = find(isnan(Cstars(:, 4))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 4))==0, 1, 'first' );% start of bistability region

Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 3), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 3));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 5));
Mlows(:,j) = vertcat(Mstars(1:bstart(j)-1, 4), Mstars(bstart(j):bend(j), 3), Mstars(bend(j)+1:end, 1));


else

%bstart(j) = find(Cstars(:, 5)==0, 1, 'last' )+1;% start of bistability region
bstart(j) = find(Cstars(:, 5)== Cstars(:, 6), 1, 'first' )+1;
bend(j) = find(isnan(Cstars(:, 7))==0, 1, 'last' );% end of bistability region


%bstart(j) = find(Mstars(:, 6)==1, 1, 'first' );% start of bistability region
%bend(j) = find(Cstars(:, 5)~= Cstars(:, 6), 1, 'last' );% end of bistability region


brk1 = find(isnan(Cstars(:, 6))==1, 1, 'last' );

Cups(:,j) = vertcat(Cstars(1:brk1, 5), Cstars(brk1+1:bend(j), 6),Cstars(bend(j)+1:end, 7)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
%Cups(:,j) = vertcat(Cstars(1:brk1, 5), Cstars(brk1+1:bend(j), 6),Cstars(bend(j)+1:end, 8)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 8), Cstars(bstart(j):bend(j), 5), Cstars(bend(j)+1:end, 7));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 8), Cstars(bstart(j):end, 4));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 7)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 8), Mstars(bstart(j):bend(j), 5), Mstars(bend(j)+1:end, 8));
Mlows(:,j) = vertcat(Mstars(1:bstart(j)-1, 8), Mstars(bstart(j):bend(j), 6), Mstars(bend(j)+1:end, 5));


end

end

toc % about 46 seconds

%% check results
%bstart

% j = 4;
% 
% plot(fset, Mlows(:,j),"-o",'Color', Mcol, 'LineWidth', 2.5) %Cups(:, ploti)
% ylim([0 1])
% hold on
% plot(fset, Mups(:,j),"-",'Color', Mcol,'LineWidth', 2.5)
% plot(fset, Mmids(:,j),"--",'Color', Mcol, 'LineWidth', 2.5) % unstable
% %plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% %plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% hold off

%% check C

% j = 4;
% 
% plot(fset, Clows(:,j),"-",'Color', Mcol, 'LineWidth', 2.5) %Cups(:, ploti)
% ylim([0 1])
% hold on
% plot(fset, Cups(:,j),"-",'Color', Mcol, 'LineWidth', 2.5)
% plot(fset, Cmids(:,j),"--o",'Color', Mcol, 'LineWidth', 2.5) % unstable
% %plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% %plot(fset21, Mmeans, '.','MarkerSize',30,'Color', Mcol)
% hold off


%% van de Leemput PDE

fset21 = linspace(0.096, 0.169, 20);

% holding arrays
Cruns = NaN(length(tset), length(xset),length(extCs),length(fset21));
Mruns = NaN(length(tset), length(xset),length(extCs),length(fset21));
Hruns = NaN(length(tset), length(xset),length(extCs),length(fset21));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(extCs),length(fset21));
Mmeans = NaN(length(extCs),length(fset21));
Hmeans = NaN(length(extCs),length(fset21));

%pksumm = NaN(6,3,length(extCs),length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

tic
for j = 1:length(extCs) % for each recruitment scenario
    % get the recruitment parameters
    beta = extCs(j);
    alpha = extMs(j);

    for i = 1:length(fset21) % for each fishing pressure
    %for i = 1:5

    ftest = fset21(i);

     % run PDE
     [solij] = altPDE(r,a,gamma,gz, rH, dH, ftest, h,d, alpha,beta,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    %Cruns(:, :, j, i) = solij(:,:,2);
    Cruns(:, :,j, i) = solij(:,:,2);
    Mruns(:, :, j, i) = solij(:,:,1);
    Hruns(:, :, j, i) = solij(:,:,3);

    % record spatial averages at final time point
    Cmeans(j, i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(j, i) = mean(solij(end, b1i:b2i, 1));
    Hmeans(j, i) = mean(solij(end, b1i:b2i, 3));

    
   
    end 

end 

toc % about 140 seconds

%% save van de Leemput model results
% to be able to plot the results for different models in the same figure

CupsL = Cups;
CmidsL = Cmids;
ClowsL = Clows;
MupsL = Mups;
MmidsL = Mmids;
MlowsL = Mlows;

fsetL = fset;
bstartL = bstart;
bendL = bend;

fset21L = fset21;

CmeansL = Cmeans;
MmeansL = Mmeans;
HmeansL = Hmeans;
   
CrunsL = Cruns;
MrunsL = Mruns;
HrunsL = Hruns;


%% plot results

Cups = CupsL;
Cmids = CmidsL;
Clows = ClowsL;
Mups = MupsL;
Mmids = MmidsL;
Mlows = MlowsL;

fset = fsetL;
bstart = bstartL;
bend = bendL;

fset21 = fset21L;

Cmeans = CmeansL;
Mmeans = MmeansL;

figure(3)
x0=10;
y0=10;
width=1500;
height=700;%375;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 4);
t.TileSpacing = 'compact';
nexttile
% start with coral
Pcol = Ccol;
Pmeans = Cmeans;
Pups = Cups;
Pmids = Cmids;
Plows = Clows;
% C and M recruitment
j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'coral cover'},'FontSize',22)
title({'a) External recruitment of'; 'coral and macroalgae'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
legend('','','','ODE, stable', 'ODE, unstable','', 'PDE, spatial means', 'Location', 'northwest')
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
nexttile
j = 2;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'b) External recruitment of'; 'coral only'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 3;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'c) External recruitment of'; 'macroalgae only'}, 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 4;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title('d) No external recruitment', 'FontSize',20)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
% now the macroalgae
Pcol = Mcol;
Pmeans = Mmeans;
Pups = Mups;
Pmids = Mmids;
Plows = Mlows;

j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'macroalgal cover'},'FontSize',22)
%title('C and M', 'FontSize',22)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
%text(0.113, 0.65, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 2;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
xlim([fset(1), fset(end)])
ylim([0 1])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 3;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
ylim([0 1])
xlim([fset(1), fset(end)])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off
nexttile
j = 4;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xline([fset(bstart(j)) fset(bend(j))]) % bistability region
ylim([0 1])
xlim([fset(1), fset(end)])
hold on
plot(fset, Pups(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Pmids(:,j),'Color', Pcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Plows(:,j),'Color', Pcol, "LineStyle","-", 'LineWidth', 2.5)
%plot the PDE means
plot(fset21, Pmeans(j,:), '.','MarkerSize',30,'Color', Pcol)
hold off

%% save output
% don't save the full runs since they take up a lot of space

save('code output/FigS10S13S14.mat','CupsB', 'CmidsB', 'ClowsB', 'MupsB', ...
    'MmidsB', 'MlowsB', 'fsetB', 'bstartB', 'bendB', 'fset21B', 'CmeansB', ...
    'MmeansB', 'CupsM', 'CmidsM', 'ClowsM', 'MupsM', ...
    'MmidsM', 'MlowsM', 'fsetM', 'bstartM', 'bendM', 'fset21M', 'CmeansM', ...
    'MmeansM','CupsL', 'CmidsL', 'ClowsL', 'MupsL', ...
    'MmidsL', 'MlowsL', 'fsetL', 'bstartL', 'bendL', 'fset21L', 'CmeansL', ...
    'MmeansL')

%% load output
% load('code output/sens2.mat','CupsB', 'CmidsB', 'ClowsB', 'MupsB', ...
%     'MmidsB', 'MlowsB', 'fsetB', 'bstartB', 'bendB', 'fset21B', 'CmeansB', ...
%     'MmeansB', 'CupsM', 'CmidsM', 'ClowsM', 'MupsM', ...
%     'MmidsM', 'MlowsM', 'fsetM', 'bstartM', 'bendM', 'fset21M', 'CmeansM', ...
%     'MmeansM','CupsL', 'CmidsL', 'ClowsL', 'MupsL', ...
%     'MmidsL', 'MlowsL', 'fsetL', 'bstartL', 'bendL', 'fset21L', 'CmeansL', ...
%     'MmeansL')



