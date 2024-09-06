% README: code for making Figures S12, S15, and S16


%% general PDE set up

% PDE parameters
diffs = [0.05,0.05,0.2,0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
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
summ10 = 1;

% get the indeces of these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));

%% Briggs model set up

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

% region of bistability
flow0 = 0.1111; % lower tipping point (calculated in Fig2.m)
fup0 = 0.1229; % upper tipping point

%% Briggs model vary taxis
%  values of fishing pressure
fset21 = linspace(0.07, 0.13, 20); 

% values of taxis
txset = linspace(-1, 1,9);

parset = txset; % parameter set 

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset));
Mruns = NaN(1, length(xset),length(fset21), length(parset));
Hruns = NaN(1, length(xset),length(fset21), length(parset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset));
Mmeans = NaN(length(fset21), length(parset));
Hmeans = NaN(length(fset21), length(parset));

tic
for k = 1:length(parset) % for each step width
   
    taxisC = parset(k);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record results
    Cruns(1, :,i, k) = solij(end,:,2);
    Mruns(1, :, i, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k) = solij(end,:,3);

    % record spatial averages at final time point
    % just take the averages from the middle to avoid edge effects
    Cmeans(i, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k) = mean(solij(end, b1i:b2i, 3));


    end 

end 

toc % most current: took 280 seconds

% save results
Cruns1 = Cruns;
Mruns1 = Mruns;
Hruns1 = Hruns;

Cmeans1 = Cmeans;
Mmeans1 = Mmeans;
Hmeans1 = Hmeans;

%% Briggs model vary herbivore diffusion

taxisC = -0.5;

diffHset = linspace(0.1, 0.5, 5);

parset = diffHset; % parameter set 

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset));
Mruns = NaN(1, length(xset),length(fset21), length(parset));
Hruns = NaN(1, length(xset),length(fset21), length(parset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset));
Mmeans = NaN(length(fset21), length(parset));
Hmeans = NaN(length(fset21), length(parset));

tic
for k = 1:length(parset) % for each step width
   
    diffs = [0.05,0.05,parset(k), 0];

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :,i, k) = solij(end,:,2);
    Mruns(1, :, i, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k) = solij(end,:,3);


    % record spatial averages at final time point
    Cmeans(i, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k) = mean(solij(end, b1i:b2i, 3));


    end 

end 

toc % took about 135 seconds



% save results
Cruns2 = Cruns;
Mruns2 = Mruns;
Hruns2 = Hruns;

Cmeans2 = Cmeans;
Mmeans2 = Mmeans;
Hmeans2 = Hmeans;

%% plot results (Fig S12)

fset21 = linspace(0.07, 0.13, 20); % for higher bistability region
txset = linspace(-1, 1,9);
diffHset = linspace(0.1, 0.5, 5);

fset = linspace(0.05, 0.15, 100);

parset = txset; 

allcols = flip(sky(length(parset(1:4)) + 1));
allcols = allcols(1:length(parset(1:4)),:);

allcols2 = gray(length(parset(5:end))+1);
allcols2 = flip(allcols2(1:length(parset(5:end)),:));

pgon = polyshape([flow0 flow0 fup0 fup0],[2 -1 -1 2]);


Mmeans = Mmeans1;

figure(1)
x0=10;
y0=10;
width=950;
height=400;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(1, 2);
t.TileSpacing = 'compact';
nexttile
for k = 1:length(parset)
    
    if k <=4
    plot(fset21,Mmeans(:,k),'o-', 'LineWidth',2, 'MarkerSize',7, 'Color', allcols(k,:)) % 4th element for transparency
    else
    plot(fset21,Mmeans(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols2(k-4,:)) % 4th element for transparency 
    end

    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
text(0.112, 0.65, 'Bistable', 'Color', 'black','FontSize', 16)
hold off
xlabel(t,'Fishing pressure','FontSize',22)
ylabel(t,{'Equilibrium';'macroalgal cover'},'FontSize',22)
title('a) Effect of herbivore taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow0 fup0]) % bistability region
lgd = legend('-1','-0.75','-0.5','-0.25','0', '0.25', '0.5', '0.75', '1','','','', 'Location','northwest', 'NumColumns',2);
title(lgd,{'Taxis towards (-) or';'away from (+) coral'})

lgd.FontSize = 14;


parset = diffHset;
allcols = pink(length(parset)+1);
allcols = flip(allcols(1:length(parset),:));

Mmeans = Mmeans2;

nexttile
for k = 1:length(parset)
    plot(fset21,Mmeans(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
title('b) Effect of herbivore diffusion', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow0 fup0]) % bistability region
lgd = legend('0.1','0.2','0.3','0.4','0.5','','','', 'Location','northwest');
title(lgd,{'Herbivore'; 'diffusion rate'})
lgd.FontSize = 14;



%% Mumby model set up

% define the symbols
syms M C H 

% parameters
r = 1; % growth rate of corals over turf
a = 0.1; % overgrowth rate of corals by macroalgae
gamma = 0.8; % overgrowth rate of turf by macroalgae
d = 0.44; % coral death rate
gz = 1; % grazing rate

% herbivore growth rate and den dep mortality: carrying capacity = rH/dH
rH = 0.2;%0.1; % herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
%f = 0; % herbivore fishing pressure

alpha = 0.01; % 0.01
beta = 0.01;

% set of fishing values
fset = linspace(0.165, 0.18, 100);

% holding vector of eq values
Cstars = NaN(length(fset), 4);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 4);

Mistars = NaN(length(fset), 4);
Mvstars = NaN(length(fset), 4);

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

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% store the tipping points
flow = fset(bstart);
fup = fset(bend);

% equilibria
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3);
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still


%% Mumby taxis simulations

% reset defaults
diffs = [0.05,0.05,0.2,0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5;%0; % taxis rate toward coral
taxisT = 0;

%  values of fishing pressure
fset21 = linspace(0.13, 0.185, 20); % for higher bistability region

% values of coral taxis
txset = linspace(-1, 1,9);

parset = txset; % parameter set 



% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset));
Mruns = NaN(1, length(xset),length(fset21), length(parset));
Hruns = NaN(1, length(xset),length(fset21), length(parset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset));
Mmeans = NaN(length(fset21), length(parset));
Hmeans = NaN(length(fset21), length(parset));

tic
for k = 1:length(parset) % for each step width
   
    taxisC = parset(k);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = MumbyHPDE(r,a,gamma,gz, rH, dH, ftest, d, alpha,beta,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record results
    Cruns(1, :,i, k) = solij(end,:,2);
    Mruns(1, :, i, k) = solij(end,:,1);
    Hruns(1, :, i, k) = solij(end,:,3);

    % record spatial averages at final time point
    % just take the averages from the middle to avoid edge effects
    Cmeans(i, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k) = mean(solij(end, b1i:b2i, 1));
    Hmeans(i, k) = mean(solij(end, b1i:b2i, 3));


    end 

end 

toc % 292 seconds


% save results
Cruns3 = Cruns;
Mruns3 = Mruns;
Hruns3 = Hruns;

Cmeans3 = Cmeans;
Mmeans3 = Mmeans;
Hmeans3 = Hmeans;


%% Mumby: herbivore diffusion


% reset defaults
diffs = [0.05,0.05,0.2,0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5;%0; % taxis rate toward coral
taxisT = 0;

diffHset = linspace(0.1, 0.5, 5);

parset = diffHset; % parameter set 

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset));
Mruns = NaN(1, length(xset),length(fset21), length(parset));
Hruns = NaN(1, length(xset),length(fset21), length(parset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset));
Mmeans = NaN(length(fset21), length(parset));
Hmeans = NaN(length(fset21), length(parset));

tic
for k = 1:length(parset) % for each step width
   
    diffs = [0.05,0.05,parset(k), 0];

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = MumbyHPDE(r,a,gamma,gz, rH, dH, ftest, d, alpha,beta,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record results
    Cruns(1, :,i, k) = solij(end,:,2);
    Mruns(1, :, i, k) = solij(end,:,1);
    Hruns(1, :, i, k) = solij(end,:,3);


    % record spatial averages at final time point
    Cmeans(i, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k) = mean(solij(end, b1i:b2i, 1));
    Hmeans(i, k) = mean(solij(end, b1i:b2i, 3));


    end 

end 

toc % took about 135 seconds

% save results
Cruns4 = Cruns;
Mruns4 = Mruns;
Hruns4 = Hruns;

Cmeans4 = Cmeans;
Mmeans4 = Mmeans;
Hmeans4 = Hmeans;


%% plot taxis and diffusion together
%  values of fishing pressure
fset21 = linspace(0.13, 0.185, 20); % for higher bistability region

% values of coral taxis
txset = linspace(-1, 1,9);
% values of herbivore diffusion
diffHset = linspace(0.1, 0.5, 5);

parset = txset; 

allcols = flip(sky(length(parset(1:4)) + 1));
allcols = allcols(1:length(parset(1:4)),:);

allcols2 = gray(length(parset(5:end))+1);
allcols2 = flip(allcols2(1:length(parset(5:end)),:));

pgon = polyshape([flow flow fup fup],[2 -1 -1 2]);

MmeansP = Mmeans3;

figure(2)
x0=10;
y0=10;
width=950;
height=400;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(1, 2);
t.TileSpacing = 'compact';
nexttile
for k = 1:length(parset)
    
    if k <=4
    plot(fset21,MmeansP(:,k),'o-', 'LineWidth',2, 'MarkerSize',7, 'Color', allcols(k,:)) % 4th element for transparency
    else
    plot(fset21,MmeansP(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols2(k-4,:)) % 4th element for transparency 
    end

    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 1])
xlim([min(fset21) max(fset21)])
text(0.167, 0.9, 'Bistable', 'Color', 'black','FontSize', 16)
hold off
xlabel(t,'Fishing pressure','FontSize',22)
ylabel(t,{'Equilibrium';'macroalgal cover'},'FontSize',22)
title('a) Effect of herbivore taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
lgd = legend('-1','-0.75','-0.5','-0.25','0', '0.25', '0.5', '0.75', '1','','','', 'Location','northwest', 'NumColumns',2);
title(lgd,{'Taxis towards (-) or';'away from (+) coral'})

lgd.FontSize = 14;

parset = diffHset;
allcols = pink(length(parset)+1);
allcols = flip(allcols(1:length(parset),:));

MmeansP = Mmeans4;

nexttile
for k = 1:length(parset)
    plot(fset21,MmeansP(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 1])
xlim([min(fset21) max(fset21)])
hold off
title('b) Effect of herbivore diffusion', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
lgd = legend('0.1','0.2','0.3','0.4','0.5','','','', 'Location','northwest');
title(lgd,{'Herbivore'; 'diffusion rate'})
lgd.FontSize = 14;


%% van de Leemput model set up

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

alpha = 0.01;% 0.025
beta = 0.01;

fset = linspace(0.125, 0.15, 100); % bistable region is around 0.145

% holding vector of eq values
Cstars = NaN(length(fset), 4);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 4);

Mistars = NaN(length(fset), 4);
Mvstars = NaN(length(fset), 4);

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

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% store the tipping points
flow2 = fset(bstart);
fup2 = fset(bend);

% equilibria
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3);
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still

%% van de Leemput: coral taxis 

%  values of fishing pressure
fset21 = linspace(0.05, 0.16, 20); % for higher bistability region

% values of coral taxis
txset = linspace(-1, 1,9);

% reset defaults
diffs = [0.05,0.05,0.2,0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5;%0; % taxis rate toward coral
taxisT = 0;

parset = txset; % parameter set 

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset));
Mruns = NaN(1, length(xset),length(fset21), length(parset));
Hruns = NaN(1, length(xset),length(fset21), length(parset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset));
Mmeans = NaN(length(fset21), length(parset));
Hmeans = NaN(length(fset21), length(parset));

tic
for k = 1:length(parset) % for each step width
   
    taxisC = parset(k);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = altPDE(r,a,gamma,gz, rH, dH, ftest, h,d, alpha,beta,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record results
    Cruns(1, :,i, k) = solij(end,:,2);
    Mruns(1, :, i, k) = solij(end,:,1);
    Hruns(1, :, i, k) = solij(end,:,3);

    % record spatial averages at final time point
    % just take the averages from the middle to avoid edge effects
    Cmeans(i, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k) = mean(solij(end, b1i:b2i, 1));
    Hmeans(i, k) = mean(solij(end, b1i:b2i, 3));


    end 

end 

toc % 205 seconds


% save results
Cruns5 = Cruns;
Mruns5 = Mruns;
Hruns5 = Hruns;

Cmeans5 = Cmeans;
Mmeans5 = Mmeans;
Hmeans5 = Hmeans;


%% van de Leemput: herbivore diffusion

%  values of fishing pressure
fset21 = linspace(0.05, 0.16, 20); % for higher bistability region


% reset defaults
diffs = [0.05,0.05,0.2,0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5; % taxis rate toward coral
taxisT = 0;

diffHset = linspace(0.1, 0.5, 5);

parset = diffHset; % parameter set 

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset));
Mruns = NaN(1, length(xset),length(fset21), length(parset));
Hruns = NaN(1, length(xset),length(fset21), length(parset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset));
Mmeans = NaN(length(fset21), length(parset));
Hmeans = NaN(length(fset21), length(parset));

tic
for k = 1:length(parset) % for each step width
   
    diffs = [0.05,0.05,parset(k), 0];

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = altPDE(r,a,gamma,gz, rH, dH, ftest, h,d, alpha,beta,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record results
    Cruns(1, :,i, k) = solij(end,:,2);
    Mruns(1, :, i, k) = solij(end,:,1);
    Hruns(1, :, i, k) = solij(end,:,3);


    % record spatial averages at final time point
    Cmeans(i, k) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k) = mean(solij(end, b1i:b2i, 1));
    Hmeans(i, k) = mean(solij(end, b1i:b2i, 3));


    end 

end 

toc % about 115 seconds

% save results
Cruns6 = Cruns;
Mruns6 = Mruns;
Hruns6 = Hruns;

Cmeans6 = Cmeans;
Mmeans6 = Mmeans;
Hmeans6 = Hmeans;



%% plot taxis and diffusion together
%  values of fishing pressure
fset21 = linspace(0.05, 0.16, 20); % for higher bistability region

% values of coral taxis
txset = linspace(-1, 1,9);
% values of herbivore diffusion
diffHset = linspace(0.1, 0.5, 5);


parset = txset; 

allcols = flip(sky(length(parset(1:4)) + 1));
allcols = allcols(1:length(parset(1:4)),:);

allcols2 = gray(length(parset(5:end))+1);
allcols2 = flip(allcols2(1:length(parset(5:end)),:));

pgon = polyshape([flow2 flow2 fup2 fup2],[2 -1 -1 2]);

MmeansP = Mmeans5;

figure(3)
x0=10;
y0=10;
width=950;
height=400;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(1, 2);
t.TileSpacing = 'compact';
nexttile
for k = 1:length(parset)
    
    if k <=4
    plot(fset21,MmeansP(:,k),'o-', 'LineWidth',2, 'MarkerSize',7, 'Color', allcols(k,:)) % 4th element for transparency
    else
    plot(fset21,MmeansP(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols2(k-4,:)) % 4th element for transparency 
    end

    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 1])
xlim([min(fset21) max(fset21)])
text(0.129, 0.9, 'Bistable', 'Color', 'black','FontSize', 16)
hold off
xlabel(t,'Fishing pressure','FontSize',22)
ylabel(t,{'Equilibrium';'macroalgal cover'},'FontSize',22)
title('a) Effect of herbivore taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow2 fup2]) % bistability region
lgd = legend('-1','-0.75','-0.5','-0.25','0', '0.25', '0.5', '0.75', '1','','','', 'Location','northwest', 'NumColumns',2);
title(lgd,{'Taxis towards (-) or';'away from (+) coral'})
lgd.FontSize = 14;


parset = diffHset;
allcols = pink(length(parset)+1);
allcols = flip(allcols(1:length(parset),:));

MmeansP = Mmeans6;

nexttile
for k = 1:length(parset)
    plot(fset21,MmeansP(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 1])
xlim([min(fset21) max(fset21)])
hold off
title('b) Effect of herbivore diffusion', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow2 fup2]) % bistability region
lgd = legend('0.1','0.2','0.3','0.4','0.5','','','', 'Location','northwest');
title(lgd,{'Herbivore'; 'diffusion rate'})
lgd.FontSize = 14;





%% save results
save('code output/FigS12S15S16.mat','Cmeans1','Mmeans1', 'Hmeans1', 'Cmeans2', ...
    'Mmeans2', 'Hmeans2', 'Cmeans3','Mmeans3', 'Hmeans3', 'Cmeans4', ...
    'Mmeans4', 'Hmeans4', 'Cmeans5','Mmeans5', 'Hmeans5', 'Cmeans6','Mmeans6', ...
    'Hmeans6')


%% load results
% load('code output/FigS12S15S16.mat','Cmeans1','Mmeans1', 'Hmeans1', 'Cmeans2', ...
%     'Mmeans2', 'Hmeans2', 'Cmeans3','Mmeans3', 'Hmeans3', 'Cmeans4', ...
%     'Mmeans4', 'Hmeans4', 'Cmeans5','Mmeans5', 'Hmeans5', 'Cmeans6','Mmeans6', ...
%     'Hmeans6')

