
% README: function for making the step function initial conditions
% C0high = high C cover (top of steps), C0low = low C cover (bottom of
% steps), C0widths = width of steps, xset = vector of x positions for
% PDE
function [initC] = stepfun(C0widths,xset) %C0high, C0low,
xlength = length(xset);

if C0widths >0 
nsteps = 1:floor(xlength/C0widths); % number of complete steps
% just want the odd ones
% odds=T(mod(T,2)~=0);
nsteps = nsteps(mod(nsteps,2)~=0);
initC = NaN(max(nsteps)*C0widths,1);

for i = nsteps
startindx = (i-1)*C0widths + 1;
stopindx = i*C0widths;
initC(startindx:stopindx) = startindx:stopindx;

end

initC = xset(initC(isnan(initC)==0));

else
    initC = NaN;
end

end