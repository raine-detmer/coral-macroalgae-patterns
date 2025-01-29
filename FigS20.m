% README: code for making Figure S20

% range of patterns as a function of external herbivore recruitment

% region of bistability changes with phiH so need to calculate the region
% of bistability of each, and then make the y-axis the distance below the
% tipping point the patterns extend expressed as a percentage of the range
% of bistability or something like that

%% get the boundaries of the region of bistability

% external recruitment set
phiHset = linspace(0, 0.2, 10);


% and record 
% ftest1 = fset2(bstart2)-0.001*fset2(bstart2); % for initial test of patterns


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


% lower bound of fishing pressures
% 1: 0-0.044, 2: 0.0667 to 0.111, 3: 0.111 to 0.2
fsetL = [0.1, 0.16, 0.22];

% upper bound of fishing pressures
fsetU = [0.2, 0.28, 0.4];


%holding arrays for upper and lower boundaries of bistability
ftps = NaN(2, length(phiHset)); % 1 = lower, 2 = upper


% turn off warning
warning('off','symbolic:numeric:NumericalInstability')


tic
for j = 1:length(phiHset)

phiH = phiHset(j);

if j < 4
    fset = linspace(fsetL(1), fsetU(1), 200);
elseif j>3 && j<7
    fset = linspace(fsetL(2), fsetU(2), 200);
else
    fset = linspace(fsetL(3), fsetU(3), 200);
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

bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

ftps(1, j) = fset(bstart); % 1 = lower, 2 = upper
ftps(2, j) = fset(bend); % 1 = lower, 2 = upper

end

toc % about 730 seconds 


%% PDE parameter set up

% PDE parameters
diffs = [0.05,0.05,0.25, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.75; % taxis rate toward coral
taxisT = 0;

diric = 0; % 0 = Neumann boundaries for constant habitat. 1 = Dirichlet boundaries for loss at the edges


% space
len = 400;
xset = linspace(-len/2,len/2,800);

% time
% default for equilibration
%t_end = 3*50000;
%tset = linspace(0,t_end,2*2500); 

% since here we just care about whether there are patterns and not whether
% they have fully equilibrated, reduce the simulation length a bit
t_end = 5000;
tset = linspace(0,t_end,2500); 

% initial conditions
icchoice = 4; % 1 = low coral, 2 = high coral, 3 = random, 4 = step function, 5 = sin function

C0high = 0.85;
C0low = 0.05;%0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 3
rnsize = 1; % magnitude of random variation (0-1)

% for icchoice = 4
C0widths = round(length(xset)/64);  % patch widths
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
%bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));



errortol = 0.0005; % error tolerance for binary search algorithm

pkN = 2; % number of peaks (in M or C) needed to count as patterns

% set of initial conditions
%parset2 = [round(length(xset)/2), round(length(xset)/64)];

% set of taxis values
txset = [-0.5, -0.75, -1];



%% get the range of fishing pressures with patterns


% reset defaults
diffs = [0.05,0.05,0.25, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisC = -0.75;%0; % taxis rate toward coral

parset = phiHset; % external herbivore recruitment
parset2 = txset; % taxis values

% holding vector for limits
flims = NaN(1, length(parset), length(parset2)); % 1 = lower limit, middle = phiH, third = taxis

tic
for z = 1:length(parset2)
% for z = 1

    taxisC = parset2(z);

for k = 1:length(parset) % for each external recruitment rate
% for k = 1
   
    phiH = parset(k);

    ftest1 = ftps(1, k)-0.005*ftps(1, k);
  
   % first test if there are patterns just past the tipping point
    ftest = ftest1;
     % run PDE
    [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

      % record peak metrics
     Cvalsijk = solij(end, :, 2);
     Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
     npks1 = mxpks;

if npks1 >= pkN % if there was at least one patch, calculate region of fishing pressures over which patches occur       
             % for lower boundary
             fstart = 0;
             fend = ftest1;

             while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % get the peak characteristics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks = mxpks;

    if(npks < pkN) % if there weren't patterns
        fstart = fmid; % fmid was too low, so make the midpoint the new lower bound
    else % if there were patterns
        fend = fmid; % fmid was too high, so make the midpoint the new upper bound
    end

             end

             flims(1,k,z) = fmid; % store this
            % pkdiff(1,k,z) = npks0-npks; % to check if there are peaks where C < M


   
end

end

   
end



   

toc % 280 seconds


flims1 = flims;


%% save results
save('code output/FigS20.mat','ftps','flims1')


%% plot results

load('code output/FigS20.mat','ftps','flims1')

% width of bistability region
bistab = ftps(2, :)-ftps(1, :); 

% extent of patterns below lower tipping point
fprop1 = (ftps(1, :)-flims1(:,:,1));
fprop2 = (ftps(1, :)-flims1(:,:,2));
fprop3 = (ftps(1, :)-flims1(:,:,3));

Mcol = [0.4667 0.6745 0.1882];
phiHset = linspace(0, 0.2, 10);

figure(1)
%plot(phiHset, fprop1, 'Col', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(phiHset, horzcat(fprop1(1:6), repelem(0, 4)), 'Col', Mcol, "LineStyle","-", 'LineWidth', 2.5)
xlabel('Rate of external herbivore recruitment (\phi_H)','FontSize',20) % t for shared label
ylabel({'Extent of Busse balloon'},'FontSize',20)
xlim([min(phiHset), max(phiHset)])
ylim([0, 1.15*max(horzcat(fprop1, fprop2, fprop3))])
hold on
%plot(phiHset(7:10), repelem(0, 4), 'Col', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(phiHset, fprop2, 'Col', Mcol, "LineStyle","--", 'LineWidth', 2.5)
plot(phiHset, fprop3, 'Col', Mcol, "LineStyle",":", 'LineWidth', 2.5)
hold off
lgd = legend('-0.5', '-0.75', '-1', 'location', 'northeast', 'FontSize',14);
title(lgd,{'Taxis towards coral'})


% % extent of patterns below lower tipping point as a fraction of the width
% % of the region of bistability
% fprop1 = (ftps(1, :)-flims1(:,:,1))./bistab;
% fprop2 = (ftps(1, :)-flims1(:,:,2))./bistab;
% fprop3 = (ftps(1, :)-flims1(:,:,3))./bistab;

% %plot(phiHset, fprop1, 'Col', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% plot(phiHset, horzcat(fprop1(1:6), repelem(0, 4)), 'Col', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% xlabel('Rate of external herbivore recruitment (\phi_H)','FontSize',20) % t for shared label
% %ylabel({'Extent of spatial patterns beyond'; 'tipping point (prop. of bistability range)'},'FontSize',20)
% ylabel({'Busse balloon extent'; '(relative to region of bistability)'},'FontSize',20)
% %ylabel({'Range of fishing pressures with patterns'; '(prop. to range of bistability)'},'FontSize',22)
% xlim([min(phiHset), max(phiHset)])
% ylim([0, max(horzcat(fprop1, fprop2, fprop3))])
% hold on
% %plot(phiHset(7:10), repelem(0, 4), 'Col', Mcol, "LineStyle","-", 'LineWidth', 2.5)
% plot(phiHset, fprop2, 'Col', Mcol, "LineStyle","--", 'LineWidth', 2.5)
% plot(phiHset, fprop3, 'Col', Mcol, "LineStyle",":", 'LineWidth', 2.5)
% hold off
% lgd = legend('-0.5', '-0.75', '-1', 'location', 'northeast', 'FontSize',14);
% title(lgd,{'Taxis towards coral'})

%% check simulations

% ftest = 0.0972-errortol;
% 
% taxisC = -0.5;%parset2(1);
% 
% phiH = 0; %parset(1);
% 
%  [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 
% 
%  plot(xset(b1i:b2i), solij(end, b1i:b2i, 2))
%  ylim([0 1])
