% README: simulations with stochastic external coral and macroalgal
% recruitment to how adding stochasticity to these parameters affects
% pattern formation

%% setup
% plotting colors
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];
Hcol = [0.9294 0.6941 0.1255]; % herbivores


%% ODE pars

% define the symbols
%syms Mi C H Mv

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
%f = 0.08; % herbivore fishing pressure


%% PDE parameter set up

% PDE parameters
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5; % taxis rate toward coral
taxisT = 0;

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges


% space
len = 400;
xset = linspace(-len/2,len/2,800);

% time
t_end = 2500;
tset = linspace(0,t_end,2500); 

% stochastic parameters
tstoch = tset(1:20:end); % approx. timepoints for stochastic recruitment
%tstoch = tset;

% boundaries of regions with different stochastic recruitment values
xstoch = [-100, 0, 100]; % xset goes from -200 to 200

% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85;
C0low = 0.05;%0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 3
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4
C0widths = round(length(xset)/16);  % patch widths
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

ftest = 0.11;%0.11;


% set of taxis values
txset = [0, -0.5];

% random replicates
rset = [5, 500, 5000];

%% step IC, stoch rec
% with and without taxis

% t_end = 10;
% tset = linspace(0,t_end,10); 
% tstoch = tset(9);


icchoice = 4; 

% holding arrays
Cruns = NaN(1, length(xset),length(txset), length(rset));
Mruns = NaN(1, length(xset), length(txset), length(rset));
Hruns = NaN(1, length(xset), length(txset), length(rset));


% outer loop: set random, stochmat
% inner loop: taxis

tic
for k = 1:length(rset) % for each replicate

    rng(rset(k)) % set seed

% calculate the stochastic recruitment matrices: random number btw 0 and 2
 %stochmatC = rand(length(tstoch), length(xstoch)+1);
stochmatC = 0 + (2-0)*rand(length(tstoch), length(xstoch)+1);


%stochmatM = rand(length(tstoch), length(xstoch)+1);
stochmatM = 0 + (2-0)*rand(length(tstoch), length(xstoch)+1);


 for mm = 1:length(txset) % for each taxis

    taxisC = txset(mm);

     % run PDE
     %rng(rset(k)) % random k
    [solij] = BriggsHrPDEStoch(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, stochmatC, stochmatM, xstoch, tstoch); 

    % record full results
    Cruns(1, :, mm, k) = solij(end,:,2);
    Mruns(1, :, mm, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, mm, k) = solij(end,:,3);
    
 end

end 

toc % 587 seconds (about 10 min)

% save results
CrunsSs = Cruns;
MrunsSs = Mruns;
HrunsSs = Hruns;

%% test plot
% plotk = 1;
% plotmm = 1;
% 
% figure(1)
% plot(xset,CrunsSs(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
% ylim([-0.01 1.75])
% xlabel('Location','FontSize',22) % t for shared label
% ylabel('Abundance','FontSize',22)
% %title('A) \tau_C = -0.96, D_C = 1.16', 'FontSize',16)
% ax = gca;
% ax.TitleHorizontalAlignment = 'left';
% hold on
% plot(xset,MrunsSs(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
% plot(xset,HrunsSs(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Hcol)
% hold off
% legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


%% check timeseries
% check to make sure stochasticity is working

% run full simulation without taxis
% taxisC = txset(1);
% 
%      % run PDE
%      %rng(rset(k)) % random k
%     [solij] = BriggsHrPDEStoch(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, stochmatC, stochmatM, xstoch, tstoch); 
% 
%     % record full results
%     CrunsTest = solij(:,:,2);
%     MrunsTest = solij(:,:,1)+ solij(:,:,4);
%     HrunsTest = solij(:,:,3);
    

% figure(2)
% plot(tset,CrunsTest(:, 100))
% hold on
% plot(tset,CrunsTest(:, 300))
% plot(tset,CrunsTest(:, 750))
% hold off
% 
% figure(3)
% plot(tset,MrunsTest(:, 100))
% hold on
% plot(tset,MrunsTest(:, 300))
% plot(tset,MrunsTest(:, 750))
% hold off

%% ran IC, stoch rec
% with and without taxis

icchoice = 3;
 

% holding arrays
Cruns = NaN(1, length(xset),length(txset), length(rset));
Mruns = NaN(1, length(xset), length(txset), length(rset));
Hruns = NaN(1, length(xset), length(txset), length(rset));


% outer loop: set random, stochmat
% inner loop: taxis

tic
for k = 1:length(rset) % for each replicate

    rng(rset(k)) % set seed

% calculate the stochastic recruitment matrices
% stochmatC = rand(length(tstoch), length(xstoch)+1);
stochmatC = 0 + (2-0)*rand(length(tstoch), length(xstoch)+1);


% stochmatM = rand(length(tstoch), length(xstoch)+1);
stochmatM = 0 + (2-0)*rand(length(tstoch), length(xstoch)+1);


 for mm = 1:length(txset) % for each taxis

    taxisC = txset(mm);

     % run PDE
     rng(rset(k)) % set seed again for random initial conditions
    [solij] = BriggsHrPDEStoch(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, stochmatC, stochmatM, xstoch, tstoch); 

    % record full results
    Cruns(1, :, mm, k) = solij(end,:,2);
    Mruns(1, :, mm, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, mm, k) = solij(end,:,3);
    
 end

end 

toc % 599 seconds (about 10 min)

% save results
CrunsRs = Cruns;
MrunsRs = Mruns;
HrunsRs = Hruns;


%% test plot
% plotk = 3;
% plotmm = 2;
% 
% figure(1)
% plot(xset,CrunsRs(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
% ylim([-0.01 1.75])
% xlabel('Location','FontSize',22) % t for shared label
% ylabel('Abundance','FontSize',22)
% %title('A) \tau_C = -0.96, D_C = 1.16', 'FontSize',16)
% ax = gca;
% ax.TitleHorizontalAlignment = 'left';
% hold on
% plot(xset,MrunsRs(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
% plot(xset,HrunsRs(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Hcol)
% hold off
% legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


%% step IC, no stoch rec
% with and without taxis

icchoice = 4; 

% holding arrays
Cruns = NaN(1, length(xset),length(txset), length(rset));
Mruns = NaN(1, length(xset), length(txset), length(rset));
Hruns = NaN(1, length(xset), length(txset), length(rset));

% for initial conditions
CrunsS0 = NaN(1, length(xset),length(txset), length(rset));
MrunsS0 = NaN(1, length(xset), length(txset), length(rset));
HrunsS0 = NaN(1, length(xset), length(txset), length(rset));


% outer loop: set random, stochmat
% inner loop: taxis

tic
%for k = 1:length(rset) % for each replicate
for k = 1

    %rng(rset(k)) % set seed

 for mm = 1:length(txset) % for each taxis

    taxisC = txset(mm);

     % run PDE
     %rng(rset(k)) % random k
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :, mm, k) = solij(end,:,2);
    Mruns(1, :, mm, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, mm, k) = solij(end,:,3);

    if mm == 1 
        % record initial conditions
    CrunsS0(1, :, mm, k) = solij(1,:,2);
    MrunsS0(1, :, mm, k) = solij(1,:,1)+ solij(1,:,4);
    HrunsS0(1, :, mm, k) = solij(1,:,3);

    end 
    
 end

end 

toc % 4 seconds

% save results
CrunsSd = Cruns;
MrunsSd = Mruns;
HrunsSd = Hruns;

%% ran IC, no stoch rec
% with and without taxis

icchoice = 3;

% holding arrays
Cruns = NaN(1, length(xset),length(txset), length(rset));
Mruns = NaN(1, length(xset), length(txset), length(rset));
Hruns = NaN(1, length(xset), length(txset), length(rset));

% for initial conditions
CrunsR0 = NaN(1, length(xset),length(txset), length(rset));
MrunsR0 = NaN(1, length(xset), length(txset), length(rset));
HrunsR0 = NaN(1, length(xset), length(txset), length(rset));


% outer loop: set random, stochmat
% inner loop: taxis

tic
for k = 1:length(rset) % for each replicate
%for k = 1

    %rng(rset(k)) % set seed

 for mm = 1:length(txset) % for each taxis

    taxisC = txset(mm);

     % run PDE
     rng(rset(k)) % random k
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :, mm, k) = solij(end,:,2);
    Mruns(1, :, mm, k) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, mm, k) = solij(end,:,3);

    if mm == 1 
        % record initial conditions
    CrunsR0(1, :, mm, k) = solij(1,:,2);
    MrunsR0(1, :, mm, k) = solij(1,:,1)+ solij(1,:,4);
    HrunsR0(1, :, mm, k) = solij(1,:,3);

    end 
    
 end

end 

toc % 8 seconds

% save results
CrunsRd = Cruns;
MrunsRd = Mruns;
HrunsRd = Hruns;

%% test plot
% plotk = 1;
% plotmm = 2;
% 
% figure(1)
% plot(xset,CrunsRd(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
% ylim([-0.01 1.75])
% xlabel('Location','FontSize',22) % t for shared label
% ylabel('Abundance','FontSize',22)
% %title('A) \tau_C = -0.96, D_C = 1.16', 'FontSize',16)
% ax = gca;
% ax.TitleHorizontalAlignment = 'left';
% hold on
% plot(xset,MrunsRd(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
% plot(xset,HrunsRd(end, :, plotmm, plotk), 'LineWidth',2, 'Color', Hcol)
% hold off
% legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


%% example of stochastic recruitment values


% holding matrixes for the timeseries in each region
phiCmat = NaN(4, length(tset));


% get the matrix for 
rng(rset(1)) % set seed

% calculate the stochastic recruitment matrices
% random number btw [a,b] = a + (b-a)*rand()
stochmatC = 0 + (2-0)*rand(length(tstoch), length(xstoch)+1);
%stochmatC = randi([0,2], [length(tstoch), length(xstoch)+1]);


% calculate the phiC values
for i = 1:length(tset)

    t = tset(i);

    % get the timepoint in tstoch closest to the current timepoint

    tstoch_near = tstoch(abs(t-tstoch)==min(abs(t-tstoch)));
    tstoch_near = tstoch_near(1);
    
    % tstoch(find(abs(t-tstoch)==min(abs(t-tstoch))));

    %tmin = max(0, tstoch_near - 5*(tset(2)-tset(1)));
   % tmax = min(tstoch_near + 5*(tset(2)-tset(1)), max(tset(end)));

   tmin = max(0, tstoch_near - 5);
   tmax = tstoch_near + 5;

   %if (t >= tstoch_near*0.99) && (t <= tstoch_near*1.01)% if t is close enough to a stochastic timepoint
   if (t >= tmin) && (t <= tmax)% if t is close enough to a stochastic timepoint
   % get the stochastic recruitment values corresponding to this element of
   % tstoch
   stochmat_near = stochmatC(tstoch==tstoch_near, :);
        
            phiCmat(i,1) = phiC*stochmat_near(1);
        
            phiCmat(i,2) = phiC*stochmat_near(2);
         
            phiCmat(i,3) = phiC*stochmat_near(3);
        
            phiCmat(i,4) = phiC*stochmat_near(4);

  else
    % if t is not close enough
            phiCmat(i,1) = phiC;
        
            phiCmat(i,2) = phiC;
         
            phiCmat(i,3) = phiC;
        
            phiCmat(i,4) = phiC;
   end


end

%% plot phiC example

% just one region
% figure(1)
% x0=10;
% y0=10;
% width=400;
% height=200;
% set(gcf,'position',[x0,y0,width,height])
% plot(tset, phiCmat(:,1))
% xlabel('Time','FontSize',18) % t for shared label
% ylabel('\phi_C','FontSize',18)

% all four regions

figure(1)
x0=10;
y0=10;
%width=400;
%height=900;
width=700;
height=600;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 2);
t.TileSpacing = 'tight';
t.TileIndexing = 'columnmajor'; % default is 'rowmajor'
nexttile
plot(tset, phiCmat(:,1))
title("Region 1", 'FontSize',18)
xlabel(t, 'Time','FontSize',22) % t for shared label
ylabel(t, 'External coral recruitment \phi_C','FontSize',22)
nexttile
plot(tset, phiCmat(:,2))
title("Region 2", 'FontSize',18)
nexttile
plot(tset, phiCmat(:,3))
title("Region 3", 'FontSize',18)
nexttile
plot(tset, phiCmat(:,4))
title("Region 4", 'FontSize',18)



%% plot initial conditions

plotk = 1;
plotmm = 1;

figure(2)
x0=10;
y0=10;
width=500;
height=100;
set(gcf,'position',[x0,y0,width,height])
plot(xset,CrunsR0(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
%xlabel('Location','FontSize',22) % t for shared label
%ylabel('Abundance','FontSize',22)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset,MrunsR0(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
plot(xset,HrunsR0(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Hcol)
hold off
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);

figure(3)
x0=10;
y0=10;
width=500;
height=100;
set(gcf,'position',[x0,y0,width,height])
plot(xset,CrunsS0(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
%xlabel('Location','FontSize',22) % t for shared label
%ylabel('Abundance','FontSize',22)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset,MrunsS0(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
plot(xset,HrunsS0(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Hcol)
hold off
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


%% plot simulation output


len = 400;
xset = linspace(-len/2,len/2,800);

xstoch = [-100, 0, 100]; 

Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];

CrunsT = CrunsSs;
MrunsT = MrunsSs;
HrunsT = HrunsSs;

figure(4)
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(3, 2);
t.TileSpacing = 'tight';
t.TileIndexing = 'columnmajor'; % default is 'rowmajor'
plotmm = 2; % with taxis
plotk  = 1; % first replicate
nexttile
plot(xset,CrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
title('Step initial conditions','FontSize',18)
text(-185, 0.9, 'Region 1', 'Color', 'black','FontSize', 16)
text(-85, 0.9, 'Region 2', 'Color', 'black','FontSize', 16)
text(15, 0.9, 'Region 3', 'Color', 'black','FontSize', 16)
text(115, 0.9, 'Region 4', 'Color', 'black','FontSize', 16)
%subtitle('no taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
xline(xstoch) % spatial regions
plot(xset,MrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
%text(-90, 1.5, 'Year 10', 'FontSize',16)
hold off

plotk  = 2; % second replicate
nexttile
plot(xset,CrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
%title('Attraction to coral','FontSize',18)
%subtitle('no taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
xline(xstoch) % spatial regions
plot(xset,MrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
%text(-90, 1.5, 'Year 10', 'FontSize',16)
hold off

plotk  = 3; % third replicate
nexttile
plot(xset,CrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
%title('Attraction to coral','FontSize',18)
%subtitle('no taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
xline(xstoch) % spatial regions
plot(xset,MrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
%text(-90, 1.5, 'Year 10', 'FontSize',16)
hold off


CrunsT = CrunsRs;
MrunsT = MrunsRs;
HrunsT = HrunsRs;

plotmm = 2; % with taxis
plotk  = 1; % first replicate
nexttile
plot(xset,CrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
title('Random initial conditions','FontSize',18)
%subtitle('no taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset,MrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
%text(-90, 1.5, 'Year 10', 'FontSize',16)
xline(xstoch) % spatial regions
hold off
legend('Coral cover', 'Macroalgal cover', 'location', 'northeast', 'FontSize',14);


plotk  = 2; % second replicate
nexttile
plot(xset,CrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
%title('Attraction to coral','FontSize',18)
%subtitle('no taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
xline(xstoch) % spatial regions
plot(xset,MrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
%text(-90, 1.5, 'Year 10', 'FontSize',16)
hold off


plotk  = 3; % third replicate
nexttile
plot(xset,CrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
%title('Attraction to coral','FontSize',18)
%subtitle('no taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
xline(xstoch) % spatial regions
plot(xset,MrunsT(1, :, plotmm, plotk), 'LineWidth',2, 'Color', Mcol)
%text(-90, 1.5, 'Year 10', 'FontSize',16)
hold off


%% save results

save('code output/RevStochExtRec.mat','CrunsSs', 'MrunsSs', 'HrunsSs', 'CrunsRs', ...
     'MrunsRs', 'HrunsRs')

%% load results

% load('code output/RevStochExtRec.mat','CrunsSs', 'MrunsSs', 'HrunsSs', 'CrunsRs', ...
%      'MrunsRs', 'HrunsRs')
