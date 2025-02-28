
% README: function for calculating the lower tipping points (boundary 
% of the region of bistability) for the nonspatial version of the Briggs
% pde model

% finitL = boundaries of the range of fishing pressures search over (vector
% where first element is the lower boundary and second element is the upper
% boundary)

% parset = model parameters (for the ODE/nonspatial version of the model)

function [lowtp] = tpfun(parset, finitL)
% output (lowtp) is the fishing pressure at the lower tipping point

% get the values of each parameter
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
phiH = parset(14);

% turn off warnings
warning('off','symbolic:numeric:NumericalInstability')

% define the symbols
syms Mi C H Mv

% set of fishing pressures to iterate over
fset = linspace(finitL(1), finitL(2), 50); 

% holding matrices for eq values
Cstars = NaN(length(fset), 4); % coral equilibria
Mstars = NaN(length(fset), 4); % macroalgal equilibria

for i = 1:length(fset)%for each element from 1 to length of fset
   
    fi = fset(i); % set fishing pressure to the ith element of fset

    % get the model eqns to solve
    eq1i = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi == 0;%Mi
    eq2i = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C ==0; %C
    eq3i = phiH + rH*H-dH*H*H-fi*H ==0; %H
    eq4i = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv ==0; % Mv
    % solve the eq values
    soli = vpasolve([eq1i, eq2i, eq3i, eq4i],[Mi,C, H, Mv], [0 Inf; 0 Inf; 0 Inf; 0 Inf]); % just pos and real equilibria
    % store the values of the eq C cover
    Cstars(i, 1:length(soli.C)) = sort(soli.C); % sort the equilibria from lowest to highest (or NA)
    Mstars(i, 1:length(soli.Mi)) = sort(soli.Mi + soli.Mv);
end

%bend = find(isnan(Cstars(:, 3))==0, 1, 'last' );% end of bistability region
bstart = find(isnan(Cstars(:, 3))==0, 1, 'first' );% start of bistability region = first row where 3rd column of Cstars is not NaN

% store the fishing pressure just before the tipping point and the
% equilibrium macroalgal cover at this fishing pressure
flow = fset(bstart-1);
Meq = Mstars(bstart-1, 1);

lowtp = [flow, Meq];


end


