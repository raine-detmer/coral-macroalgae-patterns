% README: code for making Figure S9

% sensitivity of the busse balloon to non-spatial parameters


%% PDE parameters
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

% from Mathematica
% for the lower tipping point:
% lowest across all parameters is 0.109 (for increase in rH)
% highest is 0.114 (for rH)

% for the upper tipping point:
% lowest is 0.12 (for rH)
% highest is 0.125 (for rH)

finitL = linspace(0.108, 0.115, 25); % difference between consecutive vals is 3e-4
finitU = linspace(0.119, 0.126, 25);

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

% default parameter set
%phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH
pardef = [0.01, 0.1, 0.4, 0.4, 0.02, 0.01, 0.5, 0.2, 2, 2, 0.4, 0.2, 0.1];
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
           parsens(i,j) = 0.99*pardef(i); % 10% decrease
       else 
           parsens(i,j) = 1.01*pardef(i); % 10% increase
       end

    end

end


% M and C, upper and lower tipping points
% holding vectors for everything
Clowtp = NaN(3,2,length(pardef)); % 3 for default, decrease, increase; 2 for nonspatial and spatial mean
Mlowtp = NaN(3,2,length(pardef));
Cuptp = NaN(3,2,length(pardef));
Muptp = NaN(3,2,length(pardef));

% get the defaults
pars = pardef;
[Dlowtp, Duptp] =  tpfun(pars, finitL, finitU); %uptp = [fup, Ceq, Meq];

ftest = Dlowtp(1);
[Dsol1] = BriggsHrPDE(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice); 

ftest = Duptp(1);
[Dsol2] = BriggsHrPDE(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice);


tic
for i = 1:length(pardef)

    for j = 1:3

        if j==1

            Clowtp(j,1,i) = Dlowtp(2);
            Clowtp(j,2,i) = mean(Dsol1(end, b1i:b2i, 2));
            Cuptp(j,1,i) = Duptp(2);
            Cuptp(j,2,i) = mean(Dsol2(end, b1i:b2i, 2));

            Mlowtp(j,1,i) = Dlowtp(3);
            Mlowtp(j,2,i) = mean(Dsol1(end, b1i:b2i, 1)+ Dsol1(end,b1i:b2i,4));
            Muptp(j,1,i) = Duptp(3);
            Muptp(j,2,i) = mean(Dsol2(end, b1i:b2i, 1)+ Dsol2(end,b1i:b2i,4));

        else

            pars = pardef;
            pars(i) = parsens(i,j-1);
            [lowtp, uptp] =  tpfun(pars, finitL, finitU); %uptp = [fup, Ceq, Meq];

            ftest = lowtp(1);
            [sol1] = BriggsHrPDE(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice); 

            ftest = uptp(1);
            [sol2] = BriggsHrPDE(pars(1), pars(2), pars(3), pars(4), pars(5), pars(6), pars(7), ...
    pars(8), pars(9), pars(10), pars(11), pars(12), pars(13), ftest,diffs,taxisM, ...
    taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ...
    ampC0, ampM0, period0, icchoice); 

            Clowtp(j,1,i) = lowtp(2);
            Clowtp(j,2,i) = mean(sol1(end, b1i:b2i, 2));
            Cuptp(j,1,i) = uptp(2);
            Cuptp(j,2,i) = mean(sol2(end, b1i:b2i, 2));

            Mlowtp(j,1,i) = lowtp(3);
            Mlowtp(j,2,i) = mean(sol1(end, b1i:b2i, 1)+ sol1(end,b1i:b2i,4));
            Muptp(j,1,i) = uptp(3);
            Muptp(j,2,i) = mean(sol2(end, b1i:b2i, 1)+ sol2(end,b1i:b2i,4));

        end

    end

end

toc % took 365 seconds (about 6min)

%% save results
save('code output/FigS9.mat','Clowtp', 'Cuptp', 'Mlowtp', 'Muptp')


%% load results
%load('code output/FigS9.mat','Clowtp', 'Cuptp', 'Mlowtp', 'Muptp')



%% plot results

%ymn = -2.5;
%ymx = 2.5;
ymn = -1.2;
ymx = 1.2;


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


%phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH
