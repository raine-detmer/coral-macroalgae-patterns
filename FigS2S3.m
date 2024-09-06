% README: code for making Fig. S2 and Fig. S3

% plot setup
Mcol = [0.4667 0.6745 0.1882];
Ccol = [0.3020 0.7451 0.9333];

flow = 0.1111; % lower tipping point (calculated in Fig3.m)
fup = 0.1229; % upper tipping point
pgon = polyshape([flow flow fup fup],[2 -1 -1 2]);


vec = [100; 0]; % nodes at which to place the colors
raw = [Mcol; Ccol];
N = 256;
CMmap = interp1(vec, raw, linspace(100,0, N),'pchip');


%% single step IC

C0high = 0.85;
C0low = 0.05;
M0high = 0.85;
M0low = 0.05;

%% parameter set up

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
rH = 0.2;
dH = 0.1; 

%% ODE equilibria

% define the symbols
syms Mi C H Mv


% set of fishing values
fset = linspace(0.05, 0.15, 100);

% holding vector of eq values
Cstars = NaN(length(fset), 4);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 4);

Mistars = NaN(length(fset), 4);
Mvstars = NaN(length(fset), 4);

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

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

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% use vertcat to concatenate vertical vectors
% look at the Cstars to figure out how to piece these together
% for f on x axis:
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3);
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));
% note ups and lows are from the coral's perspective still


% note this gave the following message:
% Warning: Solution does not reach accuracy goal. Verify solution manually to determine if it is
%satisfactory. 
% I spot checked a few near the tipping point in Mathematica and they seem to be accurate

%w =  warning('query','last') 
%w.identifier % figure out how the warning message is
%identified

% turn it off
%warning('off','symbolic:numeric:NumericalInstability')


%% PDE set up

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

%C0high = 0.85;
%C0low = 0.05;%0.05;
%M0high = 0.85;
%M0low = 0.05;

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
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));


% initial conditions set
initCset = NaN(5, length(xset));
%initCset(1, :) = repelem(-9999, length(cxset)); % no high coral anywhere
initCset(1, 1) = -9999; % no high coral anywhere
% initCset(2, 1:200) = xset(1:200); % 200 = 25% of xset
% initCset(3, 1:400) = xset(1:400); % 400 = 50% of xset
% initCset(4, 1:600) = xset(1:600); % 600 = 75% of xset
initCset(2, 1:40) = xset(1:40); % 40 = 5% of xset
initCset(3, 1:400) = xset(1:400); % 400 = 50% of xset
initCset(4, 1:760) = xset(1:760); % 760 = 95% of xset
initCset(5, :) = xset; % all high coral

%  values of fishing pressure
fset21 = linspace(0.07, 0.13, 20); % for higher bistability region


%% check initial conditions

% taxisC = 0.5;
% 
% k = 4;
% initC = initCset(k, find(isnan(initCset(k,:))==0));
% 
%  %ftest = fset21(14);
%  ftest = fset21(15);
%      % run PDE
%     [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 
% 
%     figure(1)
%     plot(xset, solij(1,:,2))
%     ylim([0 1])
% 
%     figure(2)
%     plot(xset, solij(1,:,1))
%     ylim([0 1])

    %% end
    % figure(1)
    % plot(xset, solij(end,:,2))
    % ylim([0 1])


%% PDE sims, different taxis values

%fset21 = [0.07];

%txset = [-1 -0.5 0 0.5 1];
txset = [-0.5 0 0.5];

%parset = [0 0.25 0.5 0.75 1]; % parameter set (% coral dominance)
parset = [0 0.05 0.5 0.95 1]; % parameter set (% coral dominance)

% holding arrays
Cruns = NaN(1, length(xset),length(fset21), length(parset), length(txset));
Mruns = NaN(1, length(xset),length(fset21), length(parset), length(txset));
Hruns = NaN(1, length(xset),length(fset21), length(parset), length(txset));

% avg abundance at each timepoint
Cmeans = NaN(length(fset21), length(parset), length(txset));
Mmeans = NaN(length(fset21), length(parset), length(txset));
Hmeans = NaN(length(fset21), length(parset), length(txset));

% initial conditions
C0runs = NaN(1, length(xset),length(fset21), length(parset), length(txset));
M0runs = NaN(1, length(xset),length(fset21), length(parset), length(txset));
H0runs = NaN(1, length(xset),length(fset21), length(parset), length(txset));



% also record the average cover initial time point (for the ODE sims)
% C0means = NaN(1, length(parset));
% Mv0means = NaN(1, length(parset));
% Mi0means = NaN(1, length(parset));
% H0means = NaN(1, length(parset));

tic
for k = 1:length(parset) % for each step width
   
    initC = initCset(k, find(isnan(initCset(k,:))==0));

 for j = 1:length(txset) % for each taxis value

     taxisC = txset(j);

    for i = 1:length(fset21) % for each fishing pressure
     %for i = 1 % for each fishing pressure

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    Cruns(1, :,i, k, j) = solij(end,:,2);
    Mruns(1, :, i, k, j) = solij(end,:,1)+ solij(end,:,4);
    Hruns(1, :, i, k, j) = solij(end,:,3);

    C0runs(1, :,i, k, j) = solij(1,:,2);
    M0runs(1, :, i, k, j) = solij(1,:,1)+ solij(1,:,4);
    H0runs(1, :, i, k, j) = solij(1,:,3);

    % record spatial averages at final time point
    % just take the averages from the middle to avoid edge effects
    Cmeans(i, k, j) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i, k, j) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i, k, j) = mean(solij(end, b1i:b2i, 3));

    % if i ==1 && j==1 % for the first fishing pressure
    %  % record the spatial averages at the first timepoint
    % C0means(k) = mean(solij(1, :, 2));
    % Mv0means(k) = mean(solij(1,:,4));
    % Mi0means(k) = mean(solij(1, :, 1));
    % H0means(k) = mean(solij(1, :, 3));
    % 
    % end


    end 
 end

end 

toc % most current: took 518 seconds (8.6 min)


% save results
Cruns1 = Cruns;
Mruns1 = Mruns;
Hruns1 = Hruns;

% record spatial averages at final time point
% just take the averages from the middle to avoid edge effects
 Cmeans1 = Cmeans;
 Mmeans1 = Mmeans;
 Hmeans1 = Hmeans;

 %% save results
 %save('code output/FigS2S3.mat','Cruns1', 'Mruns1', 'Hruns1', 'Cmeans1', ...
   % 'Mmeans1', 'Hmeans1', 'C0runs', 'M0runs', 'H0runs')


 


%% check distribution

k = 1;
ptx = 1;
pf = 14;
plot(xset,Cruns1(1, :, pf, k, ptx), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
xlabel('Location','FontSize',14) % t for shared label
ylabel('Abundance','FontSize',14)
title('0% coral dominance', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset,Mruns1(1, :, pf, k, ptx), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset,Hruns1(1, :, pf, k, ptx), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

%% check when patterns start

%t_end = 3*50000;
%tset = linspace(0,t_end,2*2500); 

t_end = 10000;
tset = linspace(0,t_end,2*2500); 

k = 1;  
initC = initCset(k, find(isnan(initCset(k,:))==0));

j = 1;
taxisC = txset(j);

 i = 14; 

        ftest = fset21(i);
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ...
        ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, ...
        C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    Ctest = solij(:,:,2);
    Mtest = solij(:,:,1)+ solij(:,:,4);
    Htest = solij(:,:,3);

% surface plot
figure(4)
surf(Ctest,'FaceAlpha',1, 'EdgeColor','none')
%colormap summer
colormap(CMmap)
view(0,90)

%% distribution plots

Hcol = [0.9294 0.6941 0.1255]; % herbivores

ttest = 2000; 

figure(5)
plot(xset,Ctest(ttest,:), 'LineWidth',2, 'Color', Ccol)
hold on
plot(xset,Mtest(ttest,:), 'LineWidth',2, 'Color', Mcol)
plot(xset,Htest(ttest,:), 'LineWidth',2, 'Color', Hcol)
hold off
ylim([-0.01 1])

% figure(6)
% plot(xset,Ctest(ttest-100,:), 'LineWidth',2, 'Color', Ccol)
% hold on
% plot(xset,Mtest(ttest-100,:), 'LineWidth',2, 'Color', Mcol)
% plot(xset,Htest(ttest-100,:), 'LineWidth',2, 'Color', Hcol)
% hold off
% ylim([-0.01 1])

% figure(6)
% plot(xset,Ctest(end,:), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
% ylim([-0.01 1])

%% plot results

% plot insets:
% https://www.mathworks.com/help/matlab/ref/axes.html
% but this messes up tiledlayout

figure(1)
k = 1;
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.7])
xlim([0.065 0.135])
title('a) 0% coral dominance', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
text(0.112, 0.6, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
scatter(fset21, Mmeans1(:, k, 1), 90, 'o','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
scatter(fset21, Mmeans1(:, k, 2), 250, '.', 'MarkerEdgeColor', Mcol)
scatter(fset21, Mmeans1(:, k, 3), 90, '^','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
hold off
% legend
hold on
lg{1} = plot(nan, 'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5);
lg{2} = plot(nan, 'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5);
lg{3} = scatter(nan, 100, 'o','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25);
lg{4} = plot(nan, '.','MarkerSize',12,'Color', Mcol);
lg{5} = scatter(nan, 100, '^','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25);
legend([lg{:}],{'ODE, stable', 'ODE, unstable','PDE means, \tau_{c} = -0.5','PDE means, \tau_{c} = 0', 'PDE means, \tau_{c} = 0.5'}, 'Location', 'southwest')
hold off
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
axes('Position',[.2 .7 .2 .2])
box on
plot(xset,C0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.2])
xlabel('Location','FontSize',14) % t for shared label
%ylabel('Abundance','FontSize',14)
hold on
plot(xset,M0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
text(-160, 1, 'Macroalgal cover', 'Color', Mcol, 'FontSize', 12)
text(-160, 0.15, 'Coral cover', 'Color', Ccol, 'FontSize', 12)
hold off
box off


%% remaining inset panels
% k = 2 through k = 5

figure(2)
k = 2; 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.7])
xlim([0.065 0.135])
title('b) 5% coral dominance', 'FontSize',16)
%title('c) 50% coral dominance', 'FontSize',16)
%title('d) 95% coral dominance', 'FontSize',16)
%title('e) 100% coral dominance', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
scatter(fset21, Mmeans1(:, k, 1), 90, 'o','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
scatter(fset21, Mmeans1(:, k, 2), 250, '.', 'MarkerEdgeColor', Mcol)
scatter(fset21, Mmeans1(:, k, 3), 90, '^','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
hold off
axes('Position',[.2 .7 .2 .2])
box on
plot(xset,C0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.2])
xlabel('Location','FontSize',14) 
hold on
plot(xset,M0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
box off


figure(3)
k = 3; 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.7])
xlim([0.065 0.135])
%title('b) 5% coral dominance', 'FontSize',16)
title('c) 50% coral dominance', 'FontSize',16)
%title('d) 95% coral dominance', 'FontSize',16)
%title('e) 100% coral dominance', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
scatter(fset21, Mmeans1(:, k, 1), 90, 'o','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
scatter(fset21, Mmeans1(:, k, 2), 250, '.', 'MarkerEdgeColor', Mcol)
scatter(fset21, Mmeans1(:, k, 3), 90, '^','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
hold off
axes('Position',[.2 .7 .2 .2])
box on
plot(xset,C0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.2])
xlabel('Location','FontSize',14) 
hold on
plot(xset,M0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
box off

figure(4)
k = 4; 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.7])
xlim([0.065 0.135])
%title('b) 5% coral dominance', 'FontSize',16)
%title('c) 50% coral dominance', 'FontSize',16)
title('d) 95% coral dominance', 'FontSize',16)
%title('e) 100% coral dominance', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
scatter(fset21, Mmeans1(:, k, 1), 90, 'o','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
scatter(fset21, Mmeans1(:, k, 2), 250, '.', 'MarkerEdgeColor', Mcol)
scatter(fset21, Mmeans1(:, k, 3), 90, '^','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
hold off
axes('Position',[.2 .7 .2 .2])
box on
plot(xset,C0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.2])
xlabel('Location','FontSize',14) 
hold on
plot(xset,M0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
box off

figure(5)
k = 5; 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
ylim([0 0.7])
xlim([0.065 0.135])
%title('b) 5% coral dominance', 'FontSize',16)
%title('c) 50% coral dominance', 'FontSize',16)
%title('d) 95% coral dominance', 'FontSize',16)
title('e) 100% coral dominance', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(fset, Mups,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5) % unstable
plot(fset, Mlows,'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5)
scatter(fset21, Mmeans1(:, k, 1), 90, 'o','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
scatter(fset21, Mmeans1(:, k, 2), 250, '.', 'MarkerEdgeColor', Mcol)
scatter(fset21, Mmeans1(:, k, 3), 90, '^','MarkerEdgeColor', Mcol, 'MarkerFaceColor', Mcol, 'MarkerFaceAlpha', 0.2, 'LineWidth', 1.25)
hold off
axes('Position',[.2 .7 .2 .2])
box on
plot(xset,C0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.2])
xlabel('Location','FontSize',14) 
hold on
plot(xset,M0runs(1, :, 1, k, 1), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
box off

%% effects of taxis on pattern dynamics

C0high = 0.85;
C0low = 0.05;
M0high = 0.85;
M0low = 0.05;

% for icchoice = 4
C0widths = round(length(xset)/64);  % step widths
initC = stepfun(C0widths, xset); 

txset = [-0.5 0 0.5];

% time
t_end = 2000;
tset = linspace(0,t_end,2000); 

% holding arrays
Cruns = NaN(length(tset), length(xset), length(txset));
Mruns = NaN(length(tset), length(xset),length(txset));
Hruns = NaN(length(tset), length(xset),length(txset));


tic

 for j = 1:length(txset) % for each taxis value

     taxisC = txset(j);

        %ftest = fset21(11); % just about at the tipping point
        ftest = fset21(13); % just before tipping point 
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    Cruns(:, :, j) = solij(:,:,2);
    Mruns(:, :, j) = solij(:,:,1)+ solij(:,:,4);
    Hruns(:, :, j) = solij(:,:,3);
   
 end

toc % most current: took 1560 seconds (26 min)


% save results
CrunsT = Cruns;
MrunsT = Mruns;
HrunsT = Hruns;

%% plot results

%tpset = [1, 2, 4, 8]; % timepoints to plot; year 0   30.0060   90.0180  210.0420

tpset = [10, 30, 80, 180]; 

figure(2)
x0=10;
y0=10;
width=700;
height=900;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(4, 3);
t.TileSpacing = 'tight';
t.TileIndexing = 'columnmajor'; % default is 'rowmajor'
plotj= 1;
tt  = tpset(1); 
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
xlabel(t, 'Location','FontSize',22) % t for shared label
ylabel(t, 'Abundance','FontSize',22)
title('Attraction to coral','FontSize',18)
%subtitle('Year 10', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
text(-90, 1.5, 'Year 10', 'FontSize',16)
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northwest', 'FontSize',14, 'NumColumns', 2);
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'southwest', 'FontSize',14);
%lgd = legend;

tt  = tpset(2); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 30', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
text(-90, 1.5, 'Year 30', 'FontSize',16)
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

tt  = tpset(3); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
% title('Year 80', 'FontSize',16)
% ax = gca;
% ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
text(-90, 1.5, 'Year 80', 'FontSize',16)
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

tt  = tpset(4); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
% title('Year 180', 'FontSize',16)
% ax = gca;
% ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
text(-90, 1.5, 'Year 180', 'FontSize',16)
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

plotj= 2;
tt  = tpset(1); 
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
title('Ignoring coral','FontSize',18)
%subtitle('Year 10', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off
legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northwest', 'FontSize',11, 'NumColumns', 1);
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'southwest', 'FontSize',14);
%lgd = legend;

tt  = tpset(2); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 30', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

tt  = tpset(3); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 80', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

tt  = tpset(4); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 180', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

plotj= 3;
tt  = tpset(1); 
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
title('Avoiding coral','FontSize',18)
%subtitle('Year 10', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off
%legend('Coral cover', 'Macroalgal cover', 'Herbivore biomass', 'location', 'northwest', 'FontSize',14, 'NumColumns', 2);

tt  = tpset(2); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 30', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

tt  = tpset(3); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 80', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off

tt  = tpset(4); % flat by 5 for attraction
nexttile
plot(xset(b1i:b2i),CrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
ylim([-0.01 1.75])
%title('Year 180', 'FontSize',16)
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';
hold on
plot(xset(b1i:b2i),MrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
hold on
plot(xset(b1i:b2i),HrunsT(tt, b1i:b2i, plotj), 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
hold off
