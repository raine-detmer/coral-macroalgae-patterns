% README: code for making Figure 2


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
f = 0.08; % herbivore fishing pressure

% calculate get the tipping point more precisely
% set of fishing values
fset2 = linspace(0.11, 0.115, 50);

% holding vector of eq values
Cstars2 = NaN(length(fset2), 4);

for i = 1:length(fset2)%for each element of gset
    % get the fishing pressure
    fi = fset2(i);

    % solve the equations
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

%% get the upper boundary of bistability

% set of fishing values
fset3 = linspace(0.12, 0.125, 50);

% holding vector of eq values
Cstars3 = NaN(length(fset3), 4);

for i = 1:length(fset3)%for each element of gset
    % get the fishing pressure
    fi = fset3(i);

    % solve the equations
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = rH*H-dH*H*H-fi*H ==0; %H
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
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisM = 0; 
taxisC = -0.5; % taxis rate toward coral
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
txset = linspace(0, 1, 9);
diffHset = linspace(0.05, 1, 9); % don't go lower than 0.05 bc that's the diff values for C and M


errortol = 0.0005; % error tolerance for binary search algorithm

pkN = 2; % number of peaks (in M or C) needed to count as patterns

ftest1 = fset2(bstart2)-0.001*fset2(bstart2); % for initial test of patterns

% set of initial conditions
parset2 = [round(length(xset)/2), round(length(xset)/64)];

% set of fishing pressures in the region of bistability
birange = flip(ftest1:0.0002:fset3(bend3));

%% get the range of fishing pressures with patterns


% reset defaults
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisC = -0.5;%0; % taxis rate toward coral

parset = txset; % parameter set 

% holding vector for limits
flims = NaN(2, length(parset), length(parset2)); % 1 = lower, 2 = upper, middle = taxis, third = initial conditions

tic
for z = 1:length(parset2)

 C0widths = parset2(z);  % step widths
initC = stepfun(C0widths, xset); 


for k = 1:length(parset) % for each step width
   
    taxisC = -1*parset(k);
  
   % first test if there are patterns just past the tipping point
    ftest = ftest1;
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

      % record peak metrics
     Cvalsijk = solij(end, :, 2);
     Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
     npks1 = mxpks;

if npks1 >= pkN % if there was at least one patch, calculate region of fishing pressures over which patches occur       
             % for lower boundary
             if k == 1 || isnan(flims(1,k-1,z)) % if this is the first taxis level with patterns
             fstart = 0;
             fend = ftest1;

             else % know that as taxis increases, the lower boundary should get lower so can make the initial upper bound lower
             fstart = 0;
             fend = flims(1,k-1,z) + errortol;
             end

             while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

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


   % now do the upper bound
   % for upper bifurcation boundary
fstart = ftest1;
fend = 1.05*fup;


while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record peak metrics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks = mxpks;

    if(npks < pkN) % if there weren't patterns
        fend = fmid; % fmid was too high, so make the midpoint the new upper bound
    else % if there were patterns
        fstart = fmid; % fmid was too low, so make the midpoint the new lower bound
    end

end
      

   flims(2,k,z) = fmid; % store this
   
end


if npks1 < pkN % if there weren't peaks at the test point
    bbend = length(birange);
    % check the whole range of bistability to find turing before tipping
    for bb = 1:length(birange)
        ftest = birange(bb);
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % get the peak characteristics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks = mxpks;

      if npks >= pkN
          break
      end

    end % end of bb for loop

    bbend = bb;
    if bbend < length(birange) % if there were patterns (the for loop broke)
        flims(2,k,z) = birange(bbend); % store the upper boundary (turing before tipping point)

        % get the lower boundary with a binary search algorithm
        fstart = ftest1;
        fend = birange(bbend);
    while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

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

    end 

end % end of if npks1 < pkN

end 

end

toc % 1741 seconds


flims1 = flims;

% test loop breaking
% for kk = 1:10
% 
%     if kk > 5
%         break
%     end
% 
% end
% 
% kk


%% repeat for diffusion

% reset defaults
diffs = [0.05,0.05,0.2, 0]; % diffusion rates, changed from diff to diffs bc otherwise diff() function doesn't work 
taxisC = -0.5;%0; % taxis rate toward coral

parset = diffHset; % parameter set 

% holding vector for limits
flims = NaN(2, length(parset), length(parset2)); % 1 = lower, 2 = upper, middle = taxis, third = initial conditions

tic
for z = 1:length(parset2)

 C0widths = parset2(z);  % step widths
initC = stepfun(C0widths, xset); 


for k = 1:length(parset) % for each step width
   
     diffs = [0.05,0.05,parset(k), 0];
   
   % first test if there are patterns just past the tipping point
    ftest = ftest1;
     % run PDE
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

      % record peak metrics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks1 = mxpks;

if npks1 >= pkN % if there was at least one patch, calculate region of fishing pressures over which patches occur
                
             % for lower boundary
             if k == 1 || isnan(flims(1,k-1,z)) % if this is the first diff level with patterns
             fstart = 0;
             fend = ftest1;

             else % know that as diff increases, the lower boundary should get higher so can make the initial lower bound higher
             fstart = flims(1,k-1,z) - errortol;
             fend = ftest1;
             end
           

             while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

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
             

   % now do the upper bound
   % for upper bifurcation boundary
fstart = ftest1;
fend = 1.05*fup;


while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % record peak metrics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
     [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks = mxpks;

    if(npks < pkN) % if there weren't patterns
        fend = fmid; % fmid was too high, so make the midpoint the new upper bound
    else % if there were patterns
        fstart = fmid; % fmid was too low, so make the midpoint the new lower bound
    end

end
      

   flims(2,k,z) = fmid; % store this
   
end


if npks1 < pkN % if there weren't peaks at the test point
    bbend = length(birange);
    % check the whole range of bistability to find turing before tipping
    for bb = 1:length(birange)
        ftest = birange(bb);
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

    % get the peak characteristics
     Cvalsijk = solij(end, :, 2);
      Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
      [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
      npks = mxpks;

      if npks >= pkN
          break
      end

    end % end of bb for loop

    bbend = bb;
    if bbend < length(birange) % if there were patterns (the for loop broke)
        flims(2,k,z) = birange(bbend); % store the upper boundary (turing before tipping point)

        % get the lower boundary with a binary search algorithm
        fstart = ftest1;
        fend = birange(bbend);
    while abs(fend-fstart) >= errortol

    fmid = (fend + fstart)/2; % calculate the fishing pressure
    ftest = fmid;
    % run the pde with this fishing pressure
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

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
           

    end 

end % end of if npks1 < pkN


end 

end

toc % 1085 seconds

flims2 = flims;


%% taxis and diffusion operating diagrams

% negative taxis (attraction to coral)

% for each value of diffusion, find the lowest value of taxis for which
% there are still patterns and do this for large and small initial patch
% widths 

parset2 = [round(length(xset)/2), round(length(xset)/64)];
%parset2 = round(length(xset)/2);

% binary search algorithm 

% reset defaults
diffs = [0.05,0.05,0.2, 0]; % diffusion rates 
taxisC = -0.5;%0; % taxis rate toward coral


ftest = fset2(bstart2)-0.005*fset2(bstart2);

diffHset3 = linspace(0.05, 1.5, 10);
% make sure 1 (max value of diffHset) and 0.2 (default) are included
diffHset3 = sort([diffHset3, 0.2, 1]);


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
    [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 

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

toc % 937 seconds


% save results
mntx1 = mntx;


%% test peaks
% taxisC = -0.1934;
% [solij] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 
% Cvalsijk = solij(end, :, 2);
%       Mvalsijk = solij(end, :, 1)+ solij(end,:,4);
%       [mxpks] = peakfun2(Cvalsijk, Mvalsijk, xset, pkthresh, b1, b2);
%       npks = mxpks;




%% plot everything together

flow = 0.1111; % lower tipping point (calculated in txdiff12)
fup = 0.1229; % upper tipping point

%fref = fset2(bstart2)-0.005*fset2(bstart2);
fref = (flow+0.000022449)-0.005*(flow+0.000022449);

C1 = [0.0118    0.6588    0.6588];
C2 = [0.1412    0.0824    0.9294];

txset = linspace(0, 1, 9);
diffHset = linspace(0.05, 1, 9); % don't go lower than 0.05 bc that's the diff values for C and M
diffHset3 = linspace(0.05, 1.5, 10);
diffHset3 = sort([diffHset3, 0.2, 1]);


mntx = mntx1;

figure(1)
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
ylim([min(diffHset3) 1.2])
%ylim([0 1.2])
xlim([min(diffHset3) 1.2])
%xlabel('Herbivore diffusion rate (m^2 yr^{-1})','FontSize',19)
%ylabel('Taxis towards coral (m^2 C^{-1} yr^{-1})','FontSize',19)
xlabel('Herbivore diffusion rate','FontSize',19)
ylabel('Taxis towards coral','FontSize',19)
text(0.2, 1.45, 'a)', 'Color', [0 0 0],'FontSize', 16)
hold on
fill(btwx, btwy2, C1, 'FaceAlpha',0.1, 'EdgeColor', C1);
text(0.6, 0.2, 'No patterns', 'Color', [0 0 0],'FontSize', 18)
text(0.45, 0.85, 'Patterns', 'Color', [0 0 0],'FontSize', 18)
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
%lnCol = [0 0 0];
line([0.2 0.2], [txset(1) txset(end)], 'Color', lnCol, 'LineStyle', '-', 'LineWidth', 1)
line([diffHset(1) diffHset(end)], [0.5 0.5], 'Color', lnCol, 'LineStyle', '-', 'LineWidth', 1)
% add markers
plot(0.2, diffHset(1), 'o','Color',lnCol,'MarkerSize',7,'MarkerEdgeColor',lnCol, 'LineWidth', 1.5)
%plot(0.2, txset(1), 'o','Color',lnCol,'MarkerSize',7,'MarkerEdgeColor',lnCol, 'LineWidth', 1.5)
plot(0.2, txset(end), 'o','Color',lnCol,'MarkerSize',7.5, 'MarkerFaceColor',lnCol, 'MarkerEdgeColor',lnCol)
plot(diffHset(1), 0.5, 'diamond','Color',lnCol,'MarkerSize',7,'MarkerEdgeColor',lnCol, 'LineWidth', 1.5)
plot(diffHset(end), 0.5, 'diamond','Color',lnCol,'MarkerSize',7.5, 'MarkerFaceColor',lnCol, 'MarkerEdgeColor',lnCol)
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
% nexttile
ax2 = nexttile(4,[1,2]); 
flims = flims1;
% fill in the region between lower and upper limits
% there can't be any NaNs here
btwx = [txset(3:end), fliplr(txset(3:end))];
btwy = [flims(1,3:end,1), fliplr(flims(2,3:end,1))];
% polygon for region of bistability
pgon = polyshape([2 -1 -1 2], [flow flow fup fup]);
plot(ax2, txset, squeeze(flims(1,:,1)),'Color',C1, "LineStyle","-", 'LineWidth', 2.5)
%xlim([min(txset) max(txset)])
xlim([min(diffHset) max(diffHset)])
%xlim([0 max(txset)])
ylim([0.06 1.1*fup])
%text(0.02, 0.13, 'b)', 'Color', [0 0 0],'FontSize', 16)
text(0.07, 0.13, 'b)', 'Color', [0 0 0],'FontSize', 16)
hold off
%text(0.01, 0.119, 'Bistable', 'Color', 'black','FontSize', 14)
text(0.08, 0.119, 'Bistable', 'Color', 'black','FontSize', 14)
hold on
ylabel('Fishing pressure','FontSize',19)
xlabel('Taxis towards coral','FontSize',19)
%ylabel('Fishing pressure (yr^{-1})','FontSize',17)
%xlabel('Taxis towards coral (m^2 C^{-1} yr^{-1})','FontSize',17)
% add the lines marking the region of bistability
yline([flow fup]) % bistability region
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
plot(txset, squeeze(flims(2,:,1)),'Color',C1, "LineStyle","-", 'LineWidth', 2.5)
fill(btwx, btwy, C1, 'FaceAlpha',0.05, 'EdgeColor', C1);
hold off
%second set of initial conditions
btwx = [txset(2:end), fliplr(txset(2:end))];
btwy = [flims(1,2:end,2), fliplr(flims(2,2:end,2))];
hold on 
plot(txset, squeeze(flims(1,:,2)),'Color',C2, "LineStyle","-", 'LineWidth', 2.5)
plot(txset, squeeze(flims(2,:,2)),'Color',C2, "LineStyle","-", 'LineWidth', 2.5)
fill(btwx, btwy, C2, 'FaceAlpha',0.05, 'EdgeColor', C2);
% add reference markers
plot(diffHset(1), fref, 'o','Color',lnCol,'MarkerSize',7,'MarkerEdgeColor',lnCol, 'LineWidth', 1.5)
%plot(txset(1), fref, 'o','Color',lnCol,'MarkerSize',7,'MarkerEdgeColor',lnCol, 'LineWidth', 1.5)
plot(txset(end), fref, 'o','Color',lnCol,'MarkerSize',7.5, 'MarkerFaceColor',lnCol, 'MarkerEdgeColor',lnCol)
hold off
% next tile
ax3 = nexttile(9,[1,2]); 
flims = flims2;
btwx = [diffHset, fliplr(diffHset)];
btwy = [flims(1,:,1), fliplr(flims(2,:,1))];
% polygon for region of bistability
pgon = polyshape([2 -1 -1 2], [flow flow fup fup]);
plot(ax3, diffHset, squeeze(flims(1,:,1)),'Color',C1, "LineStyle","-", 'LineWidth', 2.5)
xlim([min(diffHset) max(diffHset)])
%xlim([0 max(diffHset)])
ylim([0.06 1.1*fup])
text(0.07, 0.13, 'c)', 'Color', [0 0 0],'FontSize', 16)
hold off
%text(0.01, 0.115, 'Bistable', 'Color', 'black','FontSize', 16)
hold on
%ylabel('Fishing pressure (yr^{-1})','FontSize',17)
%xlabel('Herbivore diffusion rate (m^2 yr^{-1})','FontSize',17)
ylabel('Fishing pressure','FontSize',19)
xlabel('Herbivore diffusion rate','FontSize',19)
% add the lines marking the region of bistability
yline([flow fup]) % bistability region
hold on 
plot(pgon,'FaceColor','black','FaceAlpha',0.025)
plot(diffHset, squeeze(flims(2,:,1)),'Color',C1, "LineStyle","-", 'LineWidth', 2.5)
fill(btwx, btwy, C1, 'FaceAlpha',0.05, 'EdgeColor', C1);
hold off
% second set of initial conditions
btwx = [diffHset, fliplr(diffHset)];
btwy = [flims(1,:,2), fliplr(flims(2,:,2))];
hold on 
plot(diffHset, squeeze(flims(1,:,2)),'Color',C2, "LineStyle","-", 'LineWidth', 2.5)
plot(diffHset, squeeze(flims(2,:,2)),'Color',C2, "LineStyle","-", 'LineWidth', 2.5)
fill(btwx, btwy, C2, 'FaceAlpha',0.05, 'EdgeColor', C2);
% add reference markers
plot(diffHset(1), fref, 'diamond','Color',lnCol,'MarkerSize',7,'MarkerEdgeColor',lnCol, 'LineWidth', 1.5)
plot(diffHset(end), fref, 'diamond','Color',lnCol,'MarkerSize',7.5, 'MarkerFaceColor',lnCol, 'MarkerEdgeColor',lnCol)
hold off

%% save everything

save('code output/Fig2.mat','flims1', 'flims2','mntx1')

%% load everything
%clear flims1
% load everything


load('code output/Fig2.mat','flims1', 'flims2','mntx1')

%load('code output/Fig3.mat','flims1', 'flims2','mntx1')
% load one thing
%load('code output/Fig3.mat','mntx1')

%% test bounds

% C0widths = parset2(2);  % step widths
% initC = stepfun(C0widths, xset); 
% 
% ii = 2;
% 
% taxisC = mntx1(1,ii,2) - errortol;
% 
% diffs = [0.05,0.05,diffHset3(ii), 0];
% 
% [soltest] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, dH, ftest,diffs,taxisM,taxisC, taxisT, diric,xset, tset,initC,C0low, C0high, M0low, M0high,rnsize, ampC0, ampM0, period0, icchoice); 
% 
% 
% figure(2)
% plot(xset, soltest(end,:,2), 'LineWidth',2, 'Color', [0.3020 0.7451 0.9333])
% hold on 
% plot(xset, soltest(end,:,1) + soltest(end,:,4), 'LineWidth',2, 'Color', [0.4667 0.6745 0.1882])
% hold off


