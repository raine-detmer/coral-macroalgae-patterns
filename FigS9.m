% README: code for making Figure S9

% sensitivity of the busse balloon to non-spatial parameters


%% PDE parameters
% default PDE parameters
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
C0widths = round(length(xset)/16);  % step widths
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

% get the indeces of these boundaries (will use these for intervals to take
% spatial averages)
b1i = find(abs(xset-b1)==min(abs(xset-b1)));
b2i = find(abs(xset-b2)==min(abs(xset-b2)));

summ10 = 1; % 1 = record peak summaries, 0 = record all peaks

%% simulations

% estimated upper and lower tipping points in Mathematica
% for the lower tipping point:
% lowest across all parameters is 0.109 (for increase in rH)
% highest is 0.114 (for rH)

% for the upper tipping point:
% lowest is 0.12 (for rH)
% highest is 0.125 (for rH)

%finitL = linspace(0.108, 0.115, 25); % difference between consecutive vals is 3e-4
%finitU = linspace(0.119, 0.126, 25);

% UPDATE: use range of values to search for the lower tipping point for each parameter set
% (based on rough estimate of tipping point values from Mathematica)
% first column = lower end of range, 2nd column = upper end of range
% default parameter values
finitLD = [0.167, 0.168]; 
% 10% increase
% phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di,rH, dH, phiH
finitLi = [0.168, 0.169; 0.171, 0.172; 0.162, 0.163; 0.156, 0.157; 0.165, 0.166; 0.166, 0.167;
   0.163, 0.164; 0.167, 0.168; 0.169, 0.17; 0.165, 0.166; 0.179, 0.18; 0.187, 0.188; 0.158, 0.159; 0.181, 0.186]; 
% 10% decrease
finitLd= [0.165, 0.166; 0.163, 0.164; 0.172, 0.173; 0.179, 0.18; 0.169, 0.17; 0.167, 0.168;
    0.171, 0.172; 0.167, 0.168; 0.165, 0.166; 0.169, 0.17; 0.154, 0.155; 0.147, 0.148; 0.176, 0.177; 0.141, 0.146]; 


% default parameters
phiC = 0.01; 
gTC = 0.1; 
gamma = 0.4; 
gTI = 0.4;
dC = 0.02;
phiM = 0.01; 
rM = 0.5; 
gTV = 0.2;
dv = 2; 
omega = 2; 
di = 0.4; 
rH = 0.2;
dH = 0.1; 
phiH = 0.05;

% default parameter set
%phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, phiH
pardef = [0.01, 0.1, 0.4, 0.4, 0.02, 0.01, 0.5, 0.2, 2, 2, 0.4, 0.2, 0.1, 0.05];
% phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di,rH, dH


% structure for each parameter: matrix with 3 rows and 2 columns, where
% first row = default spatial + nonspatial eq, 2nd = decrease spatial +
% nonspatial eq, and 3rd = increase spatial + nonspatial
% and then have this repeated for all the parameters
% and only want to run the defaults once

% get the parameter sets
parsens = NaN(length(pardef),2);

for i = 1:length(pardef)

    for j = 1:2

       if j == 1
           %parsens(i,j) = 0.99*pardef(i); % 10% decrease
           parsens(i,j) = 0.9*pardef(i); % 10% decrease
       else 
           %parsens(i,j) = 1.01*pardef(i); % 10% increase
           parsens(i,j) = 1.1*pardef(i); % 10% increase
       end

    end

end


% M and C, upper and lower tipping points
% holding vectors for everything
%Clowtp = NaN(3,2,length(pardef)); % 3 for default, decrease, increase; 2 for nonspatial and spatial mean
Mlowtp = NaN(3,2,length(pardef));

% fishing pressures used
flowtp = NaN(3, 1,length(pardef));

% get the defaults
pars = pardef;
[Dlowtp] =  tpfun(pars, finitLD); % output is [fishing pressure just below tipping point, nonspatial eq C cover at this fishing pressure, and nonspatial M cover at this fishing pressure]

ftest = Dlowtp(1);

% ftest = Dlowtp(1)-0.005;
[Dsol1] = BriggsHrPDEextH(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice, pars(14)); 

%% plot test

% figure(1)
% plot(xset, Dsol1(end,:,2), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
% hold on 
% plot(xset, Dsol1(end,:,1) + Dsol1(end,:,4), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
% hold off

% check individual
% 1 = phiC, 6 = phiM
% i = 1;
% pars = pardef;
% pars(i) = parsens(i,1);% 1st col = 10% dec, 2nd col = 10% inc
%             % get the range of fishing pressures
% 
% % decrease
% finitL = finitLd(i, :);
% [lowtp] =  tpfun(pars, finitL); %uptp = [fup, Ceq, Meq];

% increase
% pars = pardef;
% pars(i) = parsens(i,2);% 1st col = 10% dec, 2nd col = 10% inc
%             % get the range of fishing pressures
% finitL = finitLi(i, :);
% [lowtp] =  tpfun(pars, finitL); %uptp = [fup, Ceq, Meq];



%% run simulations


tic
%for i = 1:length(pardef)
    for i = 14


    % test without phiC (1) and phiM (6)
    %if i ~= 1 && j ~= 6

    for j = 1:3

        if j==1 % default value

            % Clowtp(j,1,i) = Dlowtp(2);
            % Clowtp(j,2,i) = mean(Dsol1(end, b1i:b2i, 2));

            Mlowtp(j,1,i) = Dlowtp(2);
            Mlowtp(j,2,i) = mean(Dsol1(end, b1i:b2i, 1)+ Dsol1(end,b1i:b2i,4));

            flowtp(j,1,i) = Dlowtp(1);

        else

            pars = pardef;
            pars(i) = parsens(i,j-1);% 1st col = 10% dec, 2nd col = 10% inc
            % get the range of fishing pressures

            if j==2 % decrease

                finitL = finitLd(i, :);
  
            else % increase
                finitL = finitLi(i, :);
            end
           
            [lowtp] =  tpfun(pars, finitL); %uptp = [fup, Ceq, Meq];

            flowtp(j,1,i) = lowtp(1);

            ftest = lowtp(1);
            [sol1] = BriggsHrPDEextH(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice, pars(14)); 

            % Clowtp(j,1,i) = lowtp(2);
            % Clowtp(j,2,i) = mean(sol1(end, b1i:b2i, 2));

            Mlowtp(j,1,i) = lowtp(2);
            Mlowtp(j,2,i) = mean(sol1(end, b1i:b2i, 1)+ sol1(end,b1i:b2i,4));

        end

    end

    %end
end

toc % took 365 seconds (about 6min)

%% save results
save('code output/FigS9.mat','Mlowtp', 'flowtp')


%% load results
%load('code output/FigS9.mat','Clowtp', 'Cuptp', 'Mlowtp', 'Muptp')



%% plot results

%ymn = -2.5;
%ymx = 2.5;
ymn = -10;%-1.2;
ymx = 10;%1.2;


Mcol = [0.4667 0.6745 0.1882];

barcol = Mcol;
edgecol = Mcol;

tpset = Mlowtp(:,:,1); % 1 = default

defht = abs(tpset(1,2)-tpset(1,1));

titlesz = 18; % font size for the title

labs = {'-10'; '10'};

figure(1)
x0=10;
y0=10;
width=600;
height=900;
set(gcf,'position',[x0,y0,width,height])
t=tiledlayout(5, 3);
t.TileSpacing = 'compact';
nexttile
tpset = Mlowtp(:,:,1); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
ylabel(t,'Percent change in balloon height relative to default','FontSize',22)
xlabel(t,'Percent change in parameter relative to default', 'FontSize',22)
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'a) \phi_C', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);%,'HorizontalAlignment','left'

nexttile
tpset = Mlowtp(:,:,2); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
%text(xl(1),yl(2),'b) g_T_C', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);
text(0.25,yl(2),'b) g_T_C','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,3); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'c) \gamma', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,4); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'d) g_T_I', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,5); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'e) d_C', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,6); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'f) \phi_M', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,7); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'g) r_M', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,8); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
%text(xl(1),yl(2),'h) g_T_V', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);
text(0.25,yl(2),'h) g_T_V', 'VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,9); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'i) d_V', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,10); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'j) \omega', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,11); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'k) d_I', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,12); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'l) r_H','HorizontalAlignment','right', 'VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,13); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'m) d_H', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

nexttile
tpset = Mlowtp(:,:,14); 
tpset2 = tpset(2:3,2)-tpset(2:3,1);
bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
ax = gca;
ax.XAxis.FontSize = 16; 
yline(0,'LineWidth',1.5)
ylim([ymn ymx])
text(xl(1),yl(2),'o) \phi_H', 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',titlesz);

%phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH

% plot individual panel
% Mcol = [0.4667 0.6745 0.1882];
% barcol = Mcol;
% edgecol = Mcol;
% tpset = Mlowtp(:,:,1); % 1 = default
% labs = {'-10'; '10'};
% defht = abs(tpset(1,2)-tpset(1,1));
% tpset = Mlowtp(:,:,12); 
% tpset2 = tpset(2:3,2)-tpset(2:3,1);
% bar(labs,(tpset2-defht)/defht*100,'FaceColor',barcol,'EdgeColor',edgecol)
% ax = gca;
% ax.XAxis.FontSize = 16; 
% yline(0,'LineWidth',1.5)
% ylim([-10 10])
% %text(xl(1),yl(2),'l) r_H','HorizontalAlignment','right', 'VerticalAlignment','top','FontSize',titlesz);


%% double check the patterns

% select parameters, then run the PDE at the value just past the tipping
% point with and without taxis and plot results


i = 12;
j = 2; % 1 = decrease, 2 = increase

pars = pardef;
pars(i) = parsens(i,j);% 1st col = 10% dec, 2nd col = 10% inc

ftest = flowtp(j+1,1,i);

pars(i)
ftest

% with taxis
[sol1t] = BriggsHrPDE(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice); 

% no taxis
[sol1nt] = BriggsHrPDE(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    0, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice); 


figure(3)
plot(xset, sol1t(end,:,2), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
hold on 
plot(xset, sol1t(end,:,1) + sol1t(end,:,4), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off

figure(4)
plot(xset, sol1nt(end,:,2), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
hold on 
plot(xset, sol1nt(end,:,1) + sol1nt(end,:,4), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
hold off
