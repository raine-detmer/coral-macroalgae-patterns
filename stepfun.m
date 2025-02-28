
% README: function for making the step function initial conditions
% C0high = high C cover (top of steps), C0low = low C cover (bottom of
% steps), C0widths = width of steps, xset = vector of x positions for
% PDE
function [initC] = stepfun(C0widths,xset) 
% output (initC) is vector with the locations of xset corresponding to C steps
% (where initial C over is high)

xlength = length(xset); % number of spatial locations

if C0widths >0 % if there are coral steps (have a width > 0)
nsteps = 1:floor(xlength/C0widths); % number of complete steps
% just want the odd ones (even = macroalgae steps)
% odds=T(mod(T,2)~=0);
nsteps = nsteps(mod(nsteps,2)~=0);
initC = NaN(max(nsteps)*C0widths,1); % holding vector to fill in

for i = nsteps % for each step
startindx = (i-1)*C0widths + 1; % beginning of ith step
stopindx = i*C0widths; % end of ith step
initC(startindx:stopindx) = startindx:stopindx; % indeces for the ith step

end

initC = xset(initC(isnan(initC)==0)); % locations in xset correponding to C steps

else
    initC = NaN;
end

end