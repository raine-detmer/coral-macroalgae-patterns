% README: code for making Figure S14

%% setup
% plotting colors
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];
Hcol = [0.9294 0.6941 0.1255]; % herbivores


%% get the lower boundary of bistability


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
f = 0; % herbivore fishing pressure
phiH = 0.05; % external recruitment

% calculate get the tipping point more precisely
% set of fishing values
fset2 = linspace(0.165, 0.169, 50);

% holding vector of eq values
Cstars2 = NaN(length(fset2), 4);

for i = 1:length(fset2)%for each element of gset
    % get the fishing pressure
    fi = fset2(i);

    % solve the equations
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
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

%% get the upper boundary of bistability

% set of fishing values
fset3 = linspace(0.185, 0.189, 50);

% holding vector of eq values
Cstars3 = NaN(length(fset3), 4);

for i = 1:length(fset3)%for each element of gset
    % get the fishing pressure
    fi = fset3(i);

    % solve the equations
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C cover
    Cstars3(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
end

% process results

% get the tipping point
bend3 = find(isnan(Cstars3(:, 3))==0, 1, 'last' );% end of bistability region
bstart3 = find(isnan(Cstars3(:, 3))==0, 1, 'first' );% start of bistability region

%% store these fishing pressures
flow = fset2(bstart2); % lower boundary of bistability
fup = fset3(bend3); % upper boundary of bistability

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
t_end = 50000;
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


% taxis and diffusion sets
%txset = linspace(0, 1, 9);
%diffHset = linspace(0.05, 1, 9); % don't go lower than 0.05 bc that's the diff values for C and M


errortol = 0.0005; % error tolerance for binary search algorithm

pkN = 2; % number of peaks (in M or C) needed to count as patterns

ftest1 = fset2(bstart2)-0.01*fset2(bstart2); % for initial test of patterns

% set of initial conditions
parset2 = [round(length(xset)/2), round(length(xset)/64)];

% set of fishing pressures in the region of bistability
%birange = flip(ftest1:0.0002:fset3(bend3));


%% tests

%txset2 = [-0.6428, -3.8589]*1.75;%*1.05;
txset2 = [-1.45, -7];%*1.05;
%diffHset2 = [1.1556, 6.6833]; 
diffHset2 = [1.2, 6.75]; 

C0widths = round(length(xset)/64);  % step widths
initC = stepfun(C0widths, xset); 

% holding arrays
Cruns = NaN(length(tset), length(xset),length(txset2));
Mruns = NaN(length(tset), length(xset),length(txset2));
Hruns = NaN(length(tset), length(xset),length(txset2));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(txset2));
Mmeans = NaN(length(txset2));
Hmeans = NaN(length(txset2));

%pksumm = NaN(6,3,length(extCs),length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

tic
for j = 1:length(txset2) % for each recruitment scenario
    % get the movement parameters
    diffs = [0.05,0.05,diffHset2(j), 0]; % diffusion rates 
    taxisC = txset2(j);%0; % taxis rate toward coral


    ftest = ftest1;

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

toc

 %% plot results

 % left panel of Fig. S14

 j = 1;

figure(1)
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

figure(2)
j = 2;
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1.75])
xlabel('Location','FontSize',22) % t for shared label
ylabel('Abundance','FontSize',22)
%title('B) \tau_C = -5.79, D_C = 8.89', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Mcol)
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Hcol)
hold off
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


%% repeat with second set of initial conditions

% takes a while to fully equilibrate so increase t

t_end = 5*50000;
tset = linspace(0,t_end,2*2500); 


C0widths = round(length(xset)/2);  % step widths
initC = stepfun(C0widths, xset); 

% holding arrays
Cruns = NaN(length(tset), length(xset),length(txset2));
Mruns = NaN(length(tset), length(xset),length(txset2));
Hruns = NaN(length(tset), length(xset),length(txset2));

% also record avg abundance for each parameter combination
Cmeans = NaN(length(txset2));
Mmeans = NaN(length(txset2));
Hmeans = NaN(length(txset2));

%pksumm = NaN(6,3,length(extCs),length(fset21)); % record characteristics of middle two peaks
% 1 = wavelength (dist btw peaks), 2 = widths, 3 = prominance, 4 = absolute
% height, 5 = number of peaks (where C>M), 6 = number of peaks even if C<M

tic
for j = 1:length(txset2) % for each recruitment scenario
    % get the movement parameters
    diffs = [0.05,0.05,diffHset2(j), 0]; % diffusion rates 
    taxisC = txset2(j);%0; % taxis rate toward coral


    ftest = ftest1;

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

toc

 %% plot results

 % left panel of Fig S14

 j = 1;

figure(3)
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
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);


figure(4)
j = 2;
plot(xset(b1i:b2i),Cruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Ccol)
ylim([-0.01 1.75])
xlabel('Location','FontSize',22) % t for shared label
ylabel('Abundance','FontSize',22)
%title('B) \tau_C = -5.79, D_C = 8.89', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),Mruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Mcol)
plot(xset(b1i:b2i),Hruns(end, b1i:b2i, j, i), 'LineWidth',2, 'Color', Hcol)
hold off
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northeast', 'FontSize',14);



%% taxis and diffusion operating diagrams

% negative taxis (attraction to coral)

% for each value of diffusion, find the lowest value of taxis for which
% there are still patterns and do this for large and small initial patch
% widths 

% reset t
t_end = 50000;
tset = linspace(0,t_end,2500); 


parset2 = [round(length(xset)/2), round(length(xset)/64)];
%parset2 = round(length(xset)/2);

% binary search algorithm 

% reset defaults
diffs = [0.05,0.05,0.25, 0]; % diffusion rates 
taxisC = -0.75;%0; % taxis rate toward coral


ftest = fset2(bstart2)-0.01*fset2(bstart2);

%diffHset3 = linspace(0.05, 1.5, 10);

diffHset3 = linspace(0.05, 10, 10);
% make 0.2 (default) is included
%diffHset3 = sort([diffHset3, 0.2, 1]);


parset = diffHset3; % parameter set 

% holding arrays
mntx = NaN(1, length(parset), length(parset2)); % min value of taxis for which there are patterns

%errortol = 0.001; % use the same error tolerance as for the fishing
%pressure range diagrams

tic
for z = 1:length(parset2)

 C0widths = parset2(z);  % step widths
initC = stepfun(C0widths, xset); 


for k = 1:length(parset) % for each step width

    diffs = [0.05,0.05,parset(k), 0];

        if k == 1 || isnan(mntx(1,k-1,z))% if this is the first diff level 
        txstart = 0;
        %txend = 5*parset(k);
        txend = -5*parset(k);

        else % know that as diff increases, the lower boundary should get higher so can make the initial lower bound higher
       % txstart = max(mntx(1,k-1,z) - 10*errortol,0); % mntx values are positive
       % txend = 5*parset(k);
        txstart = min(mntx(1,k-1,z) + 10*errortol,0); % mntx values are neg
        txend = -5*parset(k);
        end

while abs(txend-txstart) >= errortol

   % taxisC = -1*(txend + txstart)/2; % calculate the value of taxis in the middle
    taxisC = (txend + txstart)/2; % calculate the value of taxis in the middle
    %taxisC
   
    % run the pde with this level of taxis
    [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record peak metrics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
      [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks = mxpks;

      %npks

    if(npks >= pkN) % if there are peaks
       % txend = -1*taxisC; % taxis was too high, so make the midpoint the new upper bound
        txend = taxisC; % taxis was too high, so make the midpoint the new upper bound
    else % if there weren't any peaks
         %txstart = -1*taxisC; % taxis was too low, so make the midpoint the new lower bound
         txstart = taxisC; % taxis was too low, so make the midpoint the new lower bound
    end
    
end

             mntx(1,k,z) = taxisC;

end 
end

toc % 1285 seconds (21 min)


% save results
mntx1 = mntx;


%% test peaks
% taxisC = -0.1934;
% [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 
% Cvalsijk = solij(end, :, 2);
%       Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
%       [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
%       npks = mxpks;




%% plot operating diagram

flow = 0.1111; % lower tipping point (calculated in txdiff12)
fup = 0.1229; % upper tipping point

%fref = fset2(bstart2)-0.005*fset2(bstart2);
%fref = (flow+0.000022449)-0.005*(flow+0.000022449);

C1 = [0.0118    0.6588    0.6588];
C2 = [0.1412    0.0824    0.9294];

%txset = linspace(0, 1, 9);
diffHset3 = linspace(0.05, 10, 10);

% test points
% txset2 = [-0.6428, -3.8589]*1.75;
% diffHset2 = [1.1556, 7.7889]; 

txset2 = [-1.45, -7];%*1.05;
diffHset2 = [1.2, 6.75]; 



mntx = mntx1;

figure(5)
x0=10;
y0=10;
width=910;
height=500;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(2, 5); % (rows, columns)
t.TileSpacing = 'compact';
ax1 = nexttile(1,[2,3]); % nexttile(x, [r,c]) means put the upper left corner of the
% axes in tile x and then make the plot span r rows and c columns
fillup = 100*repelem(max(-mntx(1, :, 1)), length(diffHset3));
filldown = repelem(0, length(diffHset3));
% fill in the region between the maxes and mins
btwx = [diffHset3, fliplr(diffHset3)];
btwy2 = [-mntx(1, :, 1), fliplr(fillup)];
plot(ax1, diffHset3, -mntx(1, :, 1),'Color', C1)
ylim([min(diffHset3) max(diffHset3)])
xlim([min(diffHset3) max(diffHset3)])
xlabel('Herbivore diffusion rate','FontSize',19)
ylabel('Taxis towards coral','FontSize',19)
text(diffHset2(1), txset2(1)*-1, 'A', 'Color', [0 0 0],'FontSize', 17)
text(diffHset2(2), txset2(2)*-1, 'B', 'Color', [0 0 0],'FontSize', 17)
hold on
fill(btwx, btwy2, C1, 'FaceAlpha',0.1, 'EdgeColor', C1);
text(5, 1, 'No patterns', 'Color', [0 0 0],'FontSize', 18)
text(2, 8, 'Patterns', 'Color', [0 0 0],'FontSize', 18)
hold off
% next initial conditions
fillup = 100*repelem(max(-mntx(1, :, 2)), length(diffHset3));
btwy2 = [-mntx(1, :, 2), fliplr(fillup)];
hold on
plot(diffHset3, -mntx(1, :, 2),'Color', C2)
fill(btwx, btwy2, C2, 'FaceAlpha',0.1, 'EdgeColor', C2);
hold off
% add the lines again and annotations so they're on top
hold on
plot(diffHset3, -mntx(1, :, 1),'Color', C1,'LineWidth', 2.5)
plot(diffHset3, -mntx(1, :, 2),'Color', C2,'LineWidth', 2.5)
lnCol = [0.9098    0.0745    0.0745];
hold off
% legend elements
hold on
lg{1} = plot(nan, 'Color', C1, "LineStyle","-", 'LineWidth', 2.5);
lg{2} = plot(nan, 'Color', C2, "LineStyle","-", 'LineWidth', 2.5);
hold off
legend([lg{:}],{'1/2', '1/64'}, 'Location', 'northeast')
lgd = legend;
title(lgd,{'Initial patch width';'(fraction total space)'})
lgd.FontSize = 14;

%% save everything

save('code output/FigS14.mat','mntx1')

%% load everything

% load('code output/FigS14.mat','mntx1')



