% README: code for making Figure 3



%% model set up

% parameters

% define the symbols
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


% first get the tipping point more precisely
% set of fishing values
fset2 = linspace(0.11, 0.115, 50);

% holding vector of eq values
Cstars2 = NaN(length(fset2), 4);

for i = 1:length(fset2)%for each element of gset
    % get the eqns
    fi = fset2(i);

    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
    Cstars2(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
end

% process results

% get the tipping point
bend2 = find(isnan(Cstars2(:, 3))==0, 1, 'last' );% end of bistability region
bstart2 = find(isnan(Cstars2(:, 3))==0, 1, 'first' );% start of bistability region

%% PDE set up

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


%% vary taxis just past tipping point and record peak metrics

% parameter set up
ftest = fset2(bstart2)-0.005*fset2(bstart2);

txset2 = linspace(0, 1, 20);

parset = txset2; % parameter set 

% initial conditions
parset2 = [round(length(xset)/2), round(length(xset)/64)];

% holding arrays
Cruns = NaN(1, length(xset),1, length(parset), length(parset2));
Mruns = NaN(1, length(xset),1, length(parset), length(parset2));
Hruns = NaN(1, length(xset),1, length(parset), length(parset2));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(1, length(parset), length(parset2));
Mmeans = NaN(1, length(parset), length(parset2));
Hmeans = NaN(1, length(parset), length(parset2));

summ10 = 1; % 1 = record peak summaries, 0 = record all peaks
pksumm = NaN(6,3,1,length(parset), length(parset2)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M


tic
for z = 1:length(parset2)

 C0widths = parset2(z);  % step widths
initC = stepfun(C0widths, xset); 


for k = 1:length(parset) % for each step width
   
    taxisC = -1*parset(k);

   for i = 1:1

     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    
    % record output
    Cruns(1, :,i, k, z) = solij(end,:,2);
    Mruns(1, :, i, k, z) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k, z) = solij(end,:,3);

    % record spatial averages at final time point
    % update to only record means in middle region
    Cmeans(i, k, z) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k, z) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k, z) = mean(solij(end, b1i:b2i, 3));

    % record peak metrics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
      [npks0, npks, pklambdas,pkwidths,pkproms,pkheights] = peakfun(Cvalsijk,Mvalsijk,summ10,xset, pkthresh, b1, b2);

            % record these
            % note: record peak density not number of peaks
            pksumm(6,1,i,k, z) = npks0/(b2-b1);
            pksumm(5,1,i, k, z) = npks/(b2-b1);
            pksumm(1,1:length(pklambdas),i,k, z) = pklambdas;
            pksumm(2,1:length(pkwidths),i,k, z) = pkwidths;
            pksumm(3,1:length(pkproms),i,k, z) = pkproms;
            pksumm(4,1:length(pkheights),i,k, z) = pkheights;


    end 

end 

end

toc % most current: 220 seconds

%% save results
Cruns1 = Cruns;
Mruns1 = Mruns;
Hruns1 = Hruns;

Cmeans1 = Cmeans;
Mmeans1 = Mmeans;
Hmeans1 = Hmeans;

pksumm1 = pksumm;


%% vary diffusion just past tipping point (for recording peak metrics)

taxisC = -0.5;

ftest = fset2(bstart2)-0.005*fset2(bstart2);

diffHset2 = linspace(0.1, 1.4, 25);

parset = diffHset2; % parameter set 


% holding arrays
Cruns = NaN(1, length(xset),1, length(parset), length(parset2));
Mruns = NaN(1, length(xset),1, length(parset), length(parset2));
Hruns = NaN(1, length(xset),1, length(parset), length(parset2));

% also record avg abundance at final timepoint for each parameter combination
Cmeans = NaN(1, length(parset), length(parset2));
Mmeans = NaN(1, length(parset), length(parset2));
Hmeans = NaN(1, length(parset), length(parset2));

summ10 = 1; % 1 = record peak summaries, 0 = record all peaks
pksumm = NaN(6,3,1,length(parset), length(parset2)); % record characteristics of middle three peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M


tic
for z = 1:length(parset2)

 C0widths = parset2(z);  % step widths
initC = stepfun(C0widths, xset); 


for k = 1:length(parset) % for each step width
   
   diffs = [0.05,0.05,parset(k), 0];

   for i = 1:1

      % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record output
    Cruns(1, :,i, k, z) = solij(end,:,2);
    Mruns(1, :, i, k, z) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k, z) = solij(end,:,3);


    % record spatial averages at final time point
    % update to only record means in middle region
    Cmeans(i, k, z) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k, z) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k, z) = mean(solij(end, b1i:b2i, 3));

    % record peak metrics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
      [npks0, npks, pklambdas,pkwidths,pkproms,pkheights] = peakfun(Cvalsijk,Mvalsijk,summ10,xset, pkthresh, b1, b2);

            % record these
            % update: record peak density not number of peaks
            pksumm(6,1,i,k, z) = npks0/(b2-b1);
            pksumm(5,1,i, k, z) = npks/(b2-b1);
            pksumm(1,1:length(pklambdas),i,k, z) = pklambdas;
            pksumm(2,1:length(pkwidths),i,k, z) = pkwidths;
            pksumm(3,1:length(pkproms),i,k, z) = pkproms;
            pksumm(4,1:length(pkheights),i,k, z) = pkheights;


    end 
% end

end 

end

toc % took 184 seconds

%% save results
Cruns2 = Cruns;
Mruns2 = Mruns;
Hruns2 = Hruns;

Cmeans2 = Cmeans;
Mmeans2 = Mmeans;
Hmeans2 = Hmeans;

pksumm2 = pksumm;


%% calculate reference C cover
% for adding a horizontal line at nonspatial equilibrium

syms Mi C H Mv

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

    % get the eqns
    fi = fset2(bstart2)-0.005*fset2(bstart2);

    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
   % Cstars(1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
   % Mstars(1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);

   Cref = soli.C;
   Cref = Cref(2);
%% plot patch characteristics together
% plot the taxis and diffusion columns on the same plot

txset2 = linspace(0, 1, 20);
diffHset2 = linspace(0.1, 1.4, 25);

%uisetcolor

C1 = [0.0118    0.6588    0.6588];
C2 = [0.1412    0.0824    0.9294];

Cref = 0.7972676185331950869975426238645; 

pksumm = pksumm1;
Cmeans = Cmeans1;

% clear gcf
% clear gca

% panel for means, number of peaks, peak wavelengths, and peak heights
figure(5)
x0=10;
y0=10;
width=600;
height=900;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 2);
t.TileIndexing = 'columnmajor'; % default is rowmajor
t.TileSpacing = 'compact';
nexttile
% start with patch density
plot(txset2, squeeze(pksumm(5,1,1,:,1)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(txset2) max(txset2)])
ylim([0 0.14])
ylabel({'Coral patch';'density'},'FontSize',22)
hold on 
plot(txset2, squeeze(pksumm(5,1,1,:,2)),'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
hold off

nexttile
% now patch widths
plot(txset2, squeeze(pksumm(2,1,1,:,1)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(txset2) max(txset2)])
ylim([0 25])
ylabel({'Coral patch';'width'},'FontSize',22)
hold on 
plot(txset2, squeeze(pksumm(2,2,1,:,1)),'Color',C1,"LineStyle","none",'Marker','.', 'LineWidth', 2.5)
plot(txset2, squeeze(pksumm(2,3,1,:,1)),'Color',C1,"LineStyle",":", 'LineWidth', 2.5)
plot(txset2, squeeze(pksumm(2,1,1,:,2)),'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
plot(txset2, squeeze(pksumm(2,2,1,:,2)),'Color',C2,"LineStyle","none",'Marker','.','LineWidth', 2.5)
plot(txset2, squeeze(pksumm(2,3,1,:,2)),'Color',C2,"LineStyle",":", 'LineWidth', 2.5)
hold off
nexttile
% now patch height (max coral cover)
plot(txset2, squeeze(pksumm(4,1,1,:,1)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(txset2) max(txset2)])
ylim([0 1])
ylabel({'Max. coral';'cover in patch'},'FontSize',22)
hold on 
plot(txset2, squeeze(pksumm(4,2,1,:,1)),'Color',C1,"LineStyle","none",'Marker','.', 'LineWidth', 2.5)
plot(txset2, squeeze(pksumm(4,3,1,:,1)),'Color',C1,"LineStyle",":", 'LineWidth', 2.5)
plot(txset2, squeeze(pksumm(4,1,1,:,2)),'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
plot(txset2, squeeze(pksumm(4,2,1,:,2)),'Color',C2,"LineStyle","none",'Marker','.','LineWidth', 2.5)
plot(txset2, squeeze(pksumm(4,3,1,:,2)),'Color',C2,"LineStyle",":", 'LineWidth', 2.5)
hold off
nexttile
% now mean coral cover
plot(txset2, Cmeans(1,:,1), 'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(txset2) max(txset2)])
yline(double(Cref), 'Linewidth',1.5)
text(0.3, 0.85,'Nonspatial equilibrium','FontSize',14)
ylim([0 1])
xlabel('Taxis towards coral','FontSize',22)
ylabel({'Mean';'Coral cover'},'FontSize',22)
hold on 
plot(txset2, Cmeans(1,:,2), 'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
hold off

% diffusion
nexttile
pksumm = pksumm2;
Cmeans = Cmeans2;
% start with patch density
plot(diffHset2, squeeze(pksumm(5,1,1,:,1)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(diffHset2) max(diffHset2)])
ylim([0 0.14])
hold on 
plot(diffHset2, squeeze(pksumm(5,1,1,:,2)),'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
hold off
lgd = legend('1/2','1/64', 'Location','northeast');
title(lgd,{'Initial patch width';'(fraction total space)'})
lgd.FontSize = 9;
nexttile
% now patch widths
plot(diffHset2, squeeze(pksumm(2,1,1,:,1)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(diffHset2) max(diffHset2)])
ylim([0 25])
hold on 
plot(diffHset2, squeeze(pksumm(2,2,1,:,1)),'Color',C1,"LineStyle","none",'Marker','.', 'LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(2,3,1,:,1)),'Color',C1,"LineStyle",":", 'LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(2,1,1,:,2)),'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(2,2,1,:,2)),'Color',C2,"LineStyle","none",'Marker','.','LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(2,3,1,:,2)),'Color',C2,"LineStyle",":", 'LineWidth', 2.5)
hold off
nexttile
% now patch height (max coral cover)
plot(diffHset2, squeeze(pksumm(4,1,1,:,1)),'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(diffHset2) max(diffHset2)])
ylim([0 1])
hold on 
plot(diffHset2, squeeze(pksumm(4,2,1,:,1)),'Color',C1,"LineStyle","none",'Marker','.', 'LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(4,3,1,:,1)),'Color',C1,"LineStyle",":", 'LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(4,1,1,:,2)),'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(4,2,1,:,2)),'Color',C2,"LineStyle","none",'Marker','.','LineWidth', 2.5)
plot(diffHset2, squeeze(pksumm(4,3,1,:,2)),'Color',C2,"LineStyle",":", 'LineWidth', 2.5)
hold off
nexttile
% now mean coral cover
plot(diffHset2, Cmeans(1,:,1), 'Color',C1,"LineStyle","-", 'LineWidth', 2.5)
xlim([min(diffHset2) max(diffHset2)])
ylim([0 1])
yline(double(Cref), 'Linewidth',1.5)
xlabel('Herbivore diffusion rate','FontSize',22)
hold on 
plot(diffHset2, Cmeans(1,:,2), 'Color',C2,"LineStyle","-", 'LineWidth', 2.5)
hold off


%% save results

save('code output/Fig3.mat','Cmeans1','pksumm1', 'Cmeans2','pksumm2')

%% load results

% load('code output/Fig3.mat','Cmeans1','pksumm1', 'Cmeans2','pksumm2')
