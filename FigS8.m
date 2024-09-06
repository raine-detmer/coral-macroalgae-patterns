% README: code for making Figure S8

%% set up
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


% PDE set up
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



%% vary amplitude of C and M

C0widths = round(length(xset)/64);
initC = stepfun(C0widths, xset); 

lowset = linspace(0,0.25, 10);
highset = 0.5-lowset;

% just look at a fishing pressure in the lower part of the balloon
fset21 = 0.10;

% holding arrays
% record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(highset));
Mmeans = NaN(length(fset21), length(highset));
Hmeans = NaN(length(fset21), length(highset));

tic
 for j = 1:length(highset) % for each step height

     C0high = highset(j);
     C0low = lowset(j);
     M0high = highset(j);
     M0low = lowset(j);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    % record spatial averages at final time point
    Cmeans(i, j) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, j) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, j) = mean(solij(end, b1i:b2i, 3));


    end 
 end

toc % took 13 seconds

% save results to plot
Mmeans1 = Mmeans;



%% vary initial M amplitude, low or high C 
%REMEMBER M + C needs to be less than or equal to 1 though... so need C to
%go between 0 and 0.5

fixedset = [0.05, 0.45];

% record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(highset), length(fixedset));
Mmeans = NaN(length(fset21), length(highset), length(fixedset));
Hmeans = NaN(length(fset21), length(highset), length(fixedset));

tic
for k = 1:length(fixedset) % for each step width
   
     C0high = fixedset(k);
     C0low = fixedset(k);

 for j = 1:length(highset) % for each step height

     M0high = highset(j);
     M0low = lowset(j);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record spatial averages at final time point
    Cmeans(i, j, k) = mean(solij(end, :, 2));
    Mmeans(i, j, k) = mean(solij(end, :, 1)+ solij(end,:,4));
    Hmeans(i, j, k) = mean(solij(end, :, 3));


    end 
 end

end 

toc % took 23 seconds

% save results to plot
Mmeans2 = Mmeans;



%% vary initial C amplitude, high or low M
%REMEMBER M + C needs to be less than or equal to 1, so need C to
%go between 0 and 0.5

fixedset = [0.05, 0.45];

% record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(highset), length(fixedset));
Mmeans = NaN(length(fset21), length(highset), length(fixedset));
Hmeans = NaN(length(fset21), length(highset), length(fixedset));

tic
for k = 1:length(fixedset) % for each step width
   
     M0high = fixedset(k);
     M0low = fixedset(k);
     
 for j = 1:length(highset) % for each step height

     C0high = highset(j);
     C0low = lowset(j);
     %M0high = highset(j);
     %M0low = lowset(j);

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record spatial averages at final time point
    Cmeans(i, j, k) = mean(solij(end, :, 2));
    Mmeans(i, j, k) = mean(solij(end, :, 1)+ solij(end,:,4));
    Hmeans(i, j, k) = mean(solij(end, :, 3));


    end 
 end

end 

toc % took 26 seconds

% save results to plot
Mmeans3 = Mmeans;


%% plot results

figure(1)
x0=10;
y0=10;
width=900;
height=300;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(1, 3);
t.TileSpacing = 'compact';
nexttile
plot(flip(highset-lowset), flip(Mmeans1(i, :)), '.-', 'MarkerSize',18,'LineWidth',2,'Color',[0.4667 0.6745 0.1882])
ylim([0 0.5])
xlim([0 0.5])
xlabel(t,'Initial patch height','FontSize',18) % t for shared label
ylabel(t,{'Equilibrium mean'; 'macroalgal cover'},'FontSize',18)
title({'a) Vary initial coral and'; 'macroalgal patch heights'},'FontSize',14)

nexttile
plot(flip(highset-lowset), flip(Mmeans3(i, :,1)), '.-', 'MarkerSize',18,'LineWidth',2,'Color',[0.4667 0.6745 0.1882])
ylim([0 0.5])
hold on 
plot(flip(highset-lowset), flip(Mmeans3(i, :,2)), '.--', 'MarkerSize',18,'LineWidth',2,'Color',[0.4667 0.6745 0.1882])
hold off
title({'b) Vary initial coral'; 'patch height'},'FontSize',14)
legend('M_0 = 0.05', 'M_0 = 0.45', 'location','northwest','FontSize',14)

nexttile
plot(flip(highset-lowset), flip(Mmeans2(i, :,1)), '.-', 'MarkerSize',18,'LineWidth',2,'Color',[0.4667 0.6745 0.1882])
ylim([0 0.5])
hold on 
plot(flip(highset-lowset), flip(Mmeans2(i, :,2)), '.--', 'MarkerSize',18,'LineWidth',2,'Color',[0.4667 0.6745 0.1882])
hold off
title({'c) Vary initial macroalgal'; 'patch height'},'FontSize',14)
legend('C_0 = 0.05', 'C_0 = 0.45', 'location','northwest','FontSize',14)


