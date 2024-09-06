%README: code for making Figure S11

% sensitivity to coral and macroalgal diffusion


%% plot colors
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];

%% set up

% bistability limits
flow = 0.1111; % lower tipping point (calculated in Fig2.m)
fup = 0.1229; % upper tipping point

% nonspatial parameters
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

%  values of fishing pressure
fset21 = linspace(0.07, 0.13, 20); % for higher bistability region

diffCM = linspace(0, 0.1,5); % default is 0.05 for both, note default Hdiff is 0.2

%% vary C diffusion

parset = diffCM; % parameter set 

txCset = [-0.5,0];% herbivore taxis

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset), length(txCset));
Mruns = NaN(1, length(xset),length(fset21), length(parset), length(txCset));
Hruns = NaN(1, length(xset),length(fset21), length(parset), length(txCset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset), length(txCset));
Mmeans = NaN(length(fset21), length(parset), length(txCset));
Hmeans = NaN(length(fset21), length(parset), length(txCset));

tic

for z = 1:length(txCset)

    taxisC = txCset(z);

for k = 1:length(parset) % for each step width
   
    diffs = [0.05,parset(k),0.2, 0];

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :,i, k, z) = solij(end,:,2);
    Mruns(1, :, i, k, z) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k, z) = solij(end,:,3);


    % record spatial averages at final time point
    Cmeans(i, k,z) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k,z) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k,z) = mean(solij(end, b1i:b2i, 3));


    end 
end

end 

toc % took about 414 seconds (7 or so minutes)

% save results
CrunsCd =Cruns;
MrunsCd = Mruns;
HrunsCd = Hruns;

CmeansCd = Cmeans;
MmeansCd = Mmeans;
HmeansCd = Hmeans;




%% vary M diffusion

parset = diffCM; % parameter set 

txCset = [-0.5,0];% herbivore taxis

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset), length(txCset));
Mruns = NaN(1, length(xset),length(fset21), length(parset), length(txCset));
Hruns = NaN(1, length(xset),length(fset21), length(parset), length(txCset));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(length(fset21), length(parset), length(txCset));
Mmeans = NaN(length(fset21), length(parset), length(txCset));
Hmeans = NaN(length(fset21), length(parset), length(txCset));

tic

for z = 1:length(txCset)

    taxisC = txCset(z);

for k = 1:length(parset) % for each step width
   
    diffs = [parset(k),0.05,0.2, 0];

    for i = 1:length(fset21) % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record full results
    Cruns(1, :,i, k, z) = solij(end,:,2);
    Mruns(1, :, i, k, z) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k, z) = solij(end,:,3);


    % record spatial averages at final time point
    Cmeans(i, k,z) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k,z) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k,z) = mean(solij(end, b1i:b2i, 3));


    end 
end

end 

toc % took about 311 seconds (5 or so minutes)

% save results
CrunsMd =Cruns;
MrunsMd = Mruns;
HrunsMd = Hruns;

CmeansMd = Cmeans;
MmeansMd = Mmeans;
HmeansMd = Hmeans;


%% plot results

fset21 = linspace(0.07, 0.13, 20); % for higher bistability region
diffCM = linspace(0, 0.1,5); % default is 0.05 for both, note default Hdiff is 0.2


parset = diffCM;

allcols = pink(length(parset)+1);
allcols = flip(allcols(1:length(parset),:));

pgon = polyshape([flow flow fup fup],[2 -1 -1 2]);


figure(1)
x0=10;
y0=10;
width=950;
height=700;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 2);
t.TileSpacing = 'compact';
nexttile
Mmeans = MmeansCd(:, :,1);
for k = 1:length(parset)
    plot(fset21,Mmeans(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
xlabel(t,'Fishing pressure','FontSize',22)
ylabel(t,'Equilibrium macroalgal cover','FontSize',22)
title('a) Coral diffusion, with herbivore taxis', 'FontSize',16)
text(0.112, 0.65, 'Bistable', 'Color', 'black','FontSize', 16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
lgd = legend('0','0.025','0.05','0.075','0.1','','','', 'Location','northwest','NumColumns',2);
title(lgd,{'Coral'; 'diffusion rate'})
lgd.FontSize = 14;
nexttile
Mmeans = MmeansCd(:, :,2);
for k = 1:length(parset)
    plot(fset21,Mmeans(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
title('b) Coral diffusion, no herbivore taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region

nexttile
Mmeans = MmeansMd(:, :,1);
for k = 1:length(parset)
    plot(fset21,Mmeans(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
xlabel(t,'Fishing pressure','FontSize',22)
ylabel(t,'Equilibrium macroalgal cover','FontSize',22)
title('c) Macroalgal diffusion, with herbivore taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region
lgd = legend('0','0.025','0.05','0.075','0.1','','','', 'Location','northwest');
title(lgd,{'Macroalgal'; 'diffusion rate'})
lgd.FontSize = 14;
nexttile
Mmeans = MmeansMd(:, :,2);
for k = 1:length(parset)
    plot(fset21,Mmeans(:,k),'.-', 'LineWidth',2, 'MarkerSize',18, 'Color', allcols(k,:)) % 4th element for transparency
    hold on;
end
hold off
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.75])
hold off
title('d) Macroalgal diffusion, no herbivore taxis', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xline([flow fup]) % bistability region





%% save results

save('code output/FigS11.mat','CmeansCd','MmeansCd', 'HmeansCd', 'CmeansMd', ...
    'MmeansMd', 'HmeansMd')

%% load results

% load('code output/FigS11.mat','CmeansCd','MmeansCd', 'HmeansCd', 'CmeansMd', ...
%     'MmeansMd', 'HmeansMd', 'CmeansCMd', 'MmeansCMd', 'HmeansCMd')
