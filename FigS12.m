% README: code for making Fig. S12


%% setup
% plotting colors
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];
Hcol = [0.9294 0.6941 0.1255]; % herbivores

% external recruitment sets
% extHs = [0, 0.001, 0.005, 0.01];

%extHs = [0, 0.001, 0.01, 0.1];
extHs = [0, 0.05, 0.1, 0.15];

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

phiM = 0.01; % default 0.0001
phiC = 0.01; % default 0.001

% herbivore parameters
rH = 0.2;%0.1; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
f = 0; % herbivore fishing pressure
phiH = 0.05; % external recruitment rate

% set of fishing values
%fset = linspace(0.05, 0.145, 100);
%fset = linspace(0.08, 0.175, 100);

fset1 = horzcat(linspace(0.05, 0.145, 100), linspace(0.145, 0.35, 20));

fset2 = horzcat(linspace(0.05, 0.16, 20), linspace(0.16, 0.35, 100));

% holding arrays for each combination of external recruitment pars
Cups = NaN(length(fset1), 4); % 4 recruitment par combinations
Cmids = NaN(length(fset1), 4);
Clows = NaN(length(fset1), 4);

Mups = NaN(length(fset1), 4);
Mmids = NaN(length(fset1), 4);
Mlows = NaN(length(fset1), 4);

bend = NaN(1,4);
bstart = NaN(1,4);

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')


tic
for j = 1:length(extHs)
%for j = 1
phiH = extHs(j);

if j < 2
    fset = fset1;
else
    fset = fset2;
end 


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
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
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

bend2 = find(isnan(Cstars(:, 2))==0, 1, 'last' ); % second equibrium switch

% use vertcat to concatenate vertical vectors
% look at the Cstars to figure out how to piece these together
% for f on x axis:
Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 2), Cstars(bstart(j):end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = Cstars(:, 3);
%Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 2));
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):bend2, 2), Cstars(bend2+1:end, 1));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 3), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 3));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 3), Mstars(bend(j)+1:end, 1));
% note ups and lows are from the coral's perspective still

else
bend(j) = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart(j) = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% use vertcat to concatenate vertical vectors
% look at the Cstars to figure out how to piece these together
% for f on x axis:
Cups(:,j) = vertcat(Cstars(1:bstart(j)-1, 1), Cstars(bstart(j):bend(j), 3), Cstars(bend(j)+1:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):bend(j), 2), Cstars(bend(j)+1:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Clows(:,j) = vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):end, 1));

Mups(:,j) = vertcat(Mstars(1:bend(j), 1), Mstars(bend(j)+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids(:,j) = vertcat(Mstars(1:bstart(j)-1, 3), Mstars(bstart(j):bend(j), 2), Mstars(bend(j)+1:end, 3));
Mlows(:,j) = vertcat(Mstars(1:bend(j), 3), Mstars(bend(j)+1:end, 1));

end

end

toc % about 136 seconds 


%% check results
%bstart
j = 3;

% plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
% ylim([0 1])
% hold on
% plot(fset, Mups(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% plot(fset, Mmids(:,j),'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
% hold off

fset = fset2;
plot(fset, Mlows(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
ylim([0 1])
%xlim([0.05 0.25])
hold on
plot(fset, Mups(:,j),'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids(:,j),'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
%plot(fset, vertcat(Cstars(1:bstart(j)-1, 4), Cstars(bstart(j):bend(j), 2), Cstars(bend(j)+1:end, 4)),'Color', Ccol, "LineStyle","--", 'LineWidth', 2.5)
hold off

%% Briggs PDE

% PDE parameters
diffs = [0.05,0.05,0.25, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.75;%0; % taxis rate toward coral
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
%fset1 = horzcat(linspace(0.05, 0.145, 100), linspace(0.146, 0.26, 20));
%fset2 = horzcat(linspace(0.05, 0.174, 20), linspace(0.175, 0.26, 100));


fset21 = linspace(0.07, 0.2, 20);
fsetL = linspace(0.07, 0.2, 20);
fsetU = linspace(0.15, 0.3, 20);

% holding arrays
Cruns = NaN(length(tset), length(xset),length(extHs),length(fset21));
Mruns = NaN(length(tset), length(xset),length(extHs),length(fset21));
Hruns = NaN(length(tset), length(xset),length(extHs),length(fset21));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(extHs),length(fset21));
Mmeans = NaN(length(extHs),length(fset21));
Hmeans = NaN(length(extHs),length(fset21));

%pksumm = NaN(6,3,length(extCs),length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

tic
for j = 1:length(extHs) % for each recruitment scenario
    % get the recruitment parameters
    phiH = extHs(j);

    if j < 3
        fset21 = fsetL;
    else
        fset21 = fsetU;
    end

    for i = 1:length(fset21) % for each fishing pressure
    %for i = 1:5

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

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

toc % 129 seconds

beep on 
beep

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

fset1 = horzcat(linspace(0.05, 0.145, 100), linspace(0.145, 0.35, 20));
fset2 = horzcat(linspace(0.05, 0.16, 20), linspace(0.16, 0.35, 100));


fset21 = linspace(0.07, 0.2, 20);
fsetL = linspace(0.07, 0.2, 20);
fsetU = linspace(0.15, 0.3, 20);


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
% no extH
fset = fset1;
fset21 = fsetL;
j = 1;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'coral cover'},'FontSize',22)
title({'a) \phi_H = 0'}, 'FontSize',20)
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
legend('','','','ODE, stable', 'ODE, unstable','', 'PDE, spatial means', 'Location', 'northeast')
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
nexttile
j = 2;
fset = fset2;
fset21 = fsetL;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'b) \phi_H = 0.05'}, 'FontSize',20)
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
fset = fset2;
fset21 = fsetU;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title({'c) \phi_H = 0.1'}, 'FontSize',20)
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
fset = fset2;
fset21 = fsetU;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
title('d) \phi_H = 0.15', 'FontSize',20)
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
fset = fset1;
fset21 = fsetL;
pgon = polyshape([fset(bstart(j)) fset(bstart(j)) fset(bend(j)) fset(bend(j))],[2 -1 -1 2]);
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel(t,'Fishing pressure (f)','FontSize',22)
ylabel({'Equilibrium'; 'macroalgal cover'},'FontSize',22)
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
fset = fset2;
fset21 = fsetL;
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
fset = fset2;
fset21 = fsetU;
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
fset = fset2;
fset21 = fsetU;
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

%% check full patterns

%j = 1:length(extHs)
%i = 1:length(fset21) 
j = 2; % 1, 2, 3, 4
i = 14;%11, 12, 13, 14

figure(2)
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1.75])
xlabel('Location','FontSize',22) % t for shared label
ylabel('Abundance','FontSize',22)
%title('A) \tau_C = -0.96, D_C = 1.16', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Mcol)
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Hcol)
hold off
legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


%% save output
% don't save the full runs since they take up a lot of space

save('code output/FigS12.mat','CupsB', 'CmidsB', 'ClowsB', 'MupsB', ...
    'MmidsB', 'MlowsB', 'fsetB', 'bstartB', 'bendB', 'fset21B', 'CmeansB', ...
    'MmeansB')

%% load output
% load('code output/FigS12.mat','CupsB', 'CmidsB', 'ClowsB', 'MupsB', ...
%     'MmidsB', 'MlowsB', 'fsetB', 'bstartB', 'bendB', 'fset21B', 'CmeansB', ...
%     'MmeansB')



