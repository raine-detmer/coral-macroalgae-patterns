% README: showing that rH doesn't have a large effect on the Busse balloon,
% other than shifting the range of fishing pressures over which
% bistability/the balloon occur

%% set up

% plot colors
Mcol = [0.4667 0.6745 0.1882]; % macroalgae
Ccol = [0.3020 0.7451 0.9333]; % coral
Hcol = [0.9294 0.6941 0.1255]; % herbivores

% region of bistability
% flow = 0.1111; % lower tipping point (calculated in Fig2.m)
% fup = 0.1229; % upper tipping point

%% ODE (non-spatial) bifurcation diagram

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
%rH = 0.2;% herbivore growth rate
dH = 0.1; % dens dep herbivore mortality
f = 0; % herbivore fishing pressure
phiH = 0.05;

% set of growth rates
rHset = [0.75*0.2, 0.2, 1.25*0.2];

% set of fishing pressure values
fset = linspace(0.09, 0.24, 200); % 200

% holding vectors for equilibrium values
Cstars = NaN(length(fset), 4); % coral
Mstars = NaN(length(fset), 4); % macroalgae (vuln + invuln)

% turn off warning
warning('off','symbolic:numeric:NumericalInstability')

% first rH
rH = rHset(1);

for i = 1:length(fset)%for each fishing pressure
    
    fi = fset(i);

    % get the eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C and M cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

% process results

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% region of bistability
flow = fset(bstart);
fup = fset(bend);

% use vertcat to concatenate vertical vectors
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3); % unstable eq
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));

% save for this rH
flow1 = flow;
fup1 = fup;
Cups1 = Cups;
Cmids1 = Cmids;
Clows1 = Clows;
Mups1 = Mups;
Mmids1 = Mmids;
Mlows1 = Mlows;


% second rH (default)
rH = rHset(2);

% holding vectors for equilibrium values
Cstars = NaN(length(fset), 4); % coral
Mstars = NaN(length(fset), 4); % macroalgae (vuln + invuln)

for i = 1:length(fset)%for each fishing pressure
    
    fi = fset(i);

    % get the eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C and M cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

% process results

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% region of bistability
flow = fset(bstart);
fup = fset(bend);

% use vertcat to concatenate vertical vectors
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3); % unstable eq
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));

% save for this rH
flow2 = flow;
fup2 = fup;
Cups2 = Cups;
Cmids2 = Cmids;
Clows2 = Clows;
Mups2 = Mups;
Mmids2 = Mmids;
Mlows2 = Mlows;


% third rH
rH = rHset(3);

% holding vectors for equilibrium values
Cstars = NaN(length(fset), 4); % coral
Mstars = NaN(length(fset), 4); % macroalgae (vuln + invuln)

for i = 1:length(fset)%for each fishing pressure
    
    fi = fset(i);

    % get the eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real
    % store the values of the eq C and M cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

% process results

% need to rearrange the eq to get smooth lines when plotting
bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% region of bistability
flow = fset(bstart);
fup = fset(bend);

% use vertcat to concatenate vertical vectors
Cups = vertcat(Cstars(1:bstart-1, 2), Cstars(bstart:end, 4)); % need to make sure the length stays the same so concatenate with NaNs from Cstars(3,)
Cmids = Cstars(:, 3); % unstable eq
Clows = vertcat(Cstars(1:bstart-1, 4), Cstars(bstart:end, 2));

Mups = vertcat(Mstars(1:bend, 1), Mstars(bend+1:end, 3)); 
Mmids = vertcat(Mstars(1:bstart-1, 3), Mstars(bstart:bend, 2), Mstars(bend+1:end, 3));
Mlows = vertcat(Mstars(1:bend, 3), Mstars(bend+1:end, 1));

% save for this rH
flow3 = flow;
fup3 = fup;
Cups3 = Cups;
Cmids3 = Cmids;
Clows3 = Clows;
Mups3 = Mups;
Mmids3 = Mmids;
Mlows3 = Mlows;

%% PDE spatial means

% PDE parameters
diffs = [0.05,0.05,0.25, 0]; % diffusion rates of MI, C, H, and Mv
taxisM = 0; 
taxisC = -0.75; % taxis rate toward coral
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
summ10 = 1; % 1 = record peak summaries, 0 = record all peaks

% get the indeces of these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));


% first rH
rH = rHset(1);

%  values of fishing pressure
fset21 = linspace(0.098, 0.14, 15); 

% record avg abundance for each parameter combination
Cmeans = NaN(length(fset21));
Mmeans = NaN(length(fset21));
Hmeans = NaN(length(fset21));

    tic
    for i = 1:length(fset21) % for each fishing pressure

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record spatial averages at final time point
    Cmeans(i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i) = mean(solij(end, b1i:b2i, 3));
   
    end 

    toc % 31 seconds

 % save results
 fset1 = fset21;
 Cmeans1 = Cmeans;
 Mmeans1 = Mmeans;
 Hmeans1 = Hmeans;
 

 % second rH
rH = rHset(2);

%  values of fishing pressure
fset21 = linspace(0.148, 0.19, 15); 

% record avg abundance for each parameter combination
Cmeans = NaN(length(fset21));
Mmeans = NaN(length(fset21));
Hmeans = NaN(length(fset21));

    tic
    for i = 1:length(fset21) % for each fishing pressure

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record spatial averages at final time point
    Cmeans(i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i) = mean(solij(end, b1i:b2i, 3));
   
    end 

    toc % 31 seconds

 % save results
 fset2 = fset21;
 Cmeans2 = Cmeans;
 Mmeans2 = Mmeans;
 Hmeans2 = Hmeans;
 

  % third rH
rH = rHset(3);

%  values of fishing pressure
fset21 = linspace(0.198, 0.24, 15); 

% record avg abundance for each parameter combination
Cmeans = NaN(length(fset21));
Mmeans = NaN(length(fset21));
Hmeans = NaN(length(fset21));

    tic
    for i = 1:length(fset21) % for each fishing pressure

    ftest = fset21(i);

     % run PDE
     [solij] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice, phiH); 

    % record spatial averages at final time point
    Cmeans(i) = mean(solij(end, b1i:b2i, 2));
    Mmeans(i) = mean(solij(end, b1i:b2i, 1)+ solij(end,b1i:b2i,4));
    Hmeans(i) = mean(solij(end, b1i:b2i, 3));
   
    end 

    toc % 31 seconds

 % save results
 fset3 = fset21;
 Cmeans3 = Cmeans;
 Mmeans3 = Mmeans;
 Hmeans3 = Hmeans;
 
 %% color picker

 uisetcolor

 %  0.1608    0.5020    0.0353

 %0.7490    0.8902    0.5882

 %% plot results

 Mcol1 = [0.7490    0.8902    0.5882];
 Mcol2 = [0.4667    0.6745    0.1882];
 Mcol3 = [0.1608    0.5020    0.0353];

figure(1)
x0=10;
y0=10;
width=600;
height=500;
set(gcf,'position',[x0,y0,width,height])
pgon = polyshape([flow1 flow1 fup1 fup1],[2 -1 -1 2]); % polygon around region of bistability
xline([flow1 fup1]) % bistability region
% start with bifurcation diagram for macroalgae
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
xlabel('Fishing pressure (f)','FontSize',22) % t for shared label
ylabel('Equilibrium macroalgal cover','FontSize',22)
%title('Sensitivity to r_H', 'FontSize',16)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylim([0 0.85])
xlim([min(fset), max(fset)])
% first rH
Mlows = Mlows1;
Mmids = Mmids1;
Mups = Mups1;
Mmeans = Mmeans1;
hold on
plot(fset, Mlows,'Color', Mcol1, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
%text(0.113, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol1, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol1, "LineStyle","--", 'LineWidth', 2.5) % unstable
%plot the PDE means
plot(fset1, Mmeans, '.','MarkerSize',30,'Color', Mcol1)
hold off
% second rH
Mlows = Mlows2;
Mmids = Mmids2;
Mups = Mups2;
Mmeans = Mmeans2;
hold on
pgon = polyshape([flow2 flow2 fup2 fup2],[2 -1 -1 2]); % polygon around region of bistability
xline([flow2 fup2]) % bistability region
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
plot(fset, Mlows,'Color', Mcol2, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
%text(0.113, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol2, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol2, "LineStyle","--", 'LineWidth', 2.5) % unstable
%plot the PDE means
plot(fset2, Mmeans, '.','MarkerSize',30,'Color', Mcol2)
hold off
% third rH
Mlows = Mlows3;
Mmids = Mmids3;
Mups = Mups3;
Mmeans = Mmeans3;
hold on
pgon = polyshape([flow3 flow3 fup3 fup3],[2 -1 -1 2]); % polygon around region of bistability
xline([flow3 fup3]) % bistability region
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
plot(fset, Mlows,'Color', Mcol3, "LineStyle","-", 'LineWidth', 2.5) %Cups(:, ploti)
%text(0.113, 0.75, 'Bistable', 'Color', 'black','FontSize', 16)
plot(fset, Mups,'Color', Mcol3, "LineStyle","-", 'LineWidth', 2.5)
plot(fset, Mmids,'Color', Mcol3, "LineStyle","--", 'LineWidth', 2.5) % unstable
%plot the PDE means
plot(fset3, Mmeans, '.','MarkerSize',30,'Color', Mcol3)
hold off
% add labels
text(0.095, 0.34, 'r_H = 0.15', 'Color', Mcol1,'FontSize', 16) %'Color', [0, 0, 0]
text(0.146, 0.34, 'r_H = 0.2', 'Color', Mcol2,'FontSize', 16)
text(0.195, 0.34, 'r_H = 0.25', 'Color', Mcol3,'FontSize', 16)
hold on
% make legend
lg{1} = plot(nan, 'Color', Mcol, "LineStyle","-", 'LineWidth', 2.5);
lg{2} = plot(nan, 'Color', Mcol, "LineStyle","--", 'LineWidth', 2.5);
lg{3} = plot(nan, '.','MarkerSize',35,'Color', Mcol);
legend([lg{:}],{'ODE, stable', 'ODE, unstable','PDE, means'}, 'Location', 'northwest')
hold off
legend('boxoff')
lgd = legend;
lgd.FontSize = 14;
