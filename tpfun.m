
% README: function for calculating the upper and lower tipping points (boundaries 
% of the region of bistability) for the nonspatial version of the Briggs
% pde model

% finit = intial set of fishing pressures to test

% input two finit sets, one for finding lower boundary and
% one for finding upper boundary 

% parset = model parameters (for the ODE/nonspatial version of the model)

function [lowtp, uptp] = tpfun(parset, finitL, finitU)
phiC = parset(1);
gTC= parset(2); 
gamma= parset(3);
gTI= parset(4); 
dC= parset(5); 
phiM = parset(6);
rM= parset(7);
gTV= parset(8); 
dv= parset(9);
omega= parset(10);
di= parset(11);
rH= parset(12);
dH= parset(13);

% turn off warnings
warning('off','symbolic:numeric:NumericalInstability')

% define the symbols
syms Mi C H Mv

fset = finitL;

% holding vector of eq values
Cstars = NaN(length(fset), 4);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 4);

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

%bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% store the fishing pressure at the tipping point and the equilibrium M and
% C covers just before the tipping point
%flow = fset(bstart);
% UPDATE: store the fishing pressure just before the tipping point
flow = fset(bstart-1);
Ceq = Cstars(bstart-1, 2);
Meq = Mstars(bstart-1, 1);

lowtp = [flow, Ceq, Meq];

% if isempty(bend) % same as if length(bend) == 0
%     fup = NaN;
%     flow = NaN;
% 
% else
%     % get more detailed locations

% repeat for the upper tipping point
fset = finitU;

% holding vector of eq values
Cstars = NaN(length(fset), 4);%not sure how many pos, real eq...maybe run a single 
% value in region of bistability to check how many solutions there were?
Mstars = NaN(length(fset), 4);

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

bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
%bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region

% store the fishing pressure at the tipping point and the equilibrium M and
% C covers just before the tipping point
%fup = fset(bend);
% UPDATE: store the fishing pressure just after the tipping point
fup = fset(bend+1);
Ceq = Cstars(bend+1, 2);
Meq = Mstars(bend+1, 1);

uptp = [fup, Ceq, Meq];

end