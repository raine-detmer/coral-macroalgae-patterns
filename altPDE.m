% README: function for running the PDE version of the Mumby model with a
% type II grazing term and dynamic herbivores (like van de Leemput et al. 2016) 


function[sol] = altPDE(r,a,gamma,gz, rH, dH, f,h, d,alpha, beta, diff,taxisM,taxisC, taxisT, diric, ...
    x,t,initC,Clow, Chigh, Mlow, Mhigh, rnsize, ampC0, ampM0, period0,icchoice)

% Find solution by simulating with PDE solver
sol = pdepe(0,@pdefcn,@pdeic,@pdebc,x,t);
%Each row of sol represents a timepoint; each column a value of x
%Sol's third dimension follows assigment in the pdepe function.


%System of PDEs describing population dynamics
    function [c,g,s] = pdefcn(xi,t,y,DyDx)
        %k = kofx(xi); % fish carrying capacity is a function (called kofx) of space
        c = [1,1,1]'; % c= coefficient of the terms that get differentiated with respect to time
        %g = [diff(1)*DyDx(1),diff(2)*DyDx(2),diff(3)*DyDx(3)+taxisM*DyDx(1)+taxisC*DyDx(2)]';
        g = [diff(1)*DyDx(1),diff(2)*DyDx(2),diff(3)*DyDx(3)+taxisM*y(3)*DyDx(1)+taxisC*y(3)*DyDx(2)+ taxisT*y(3)*(-DyDx(1)-DyDx(2))]'; % things that get differentiated w/r/t space
        %third component describes changes in fish over space: have diffusion, could have taxis with respect to M and C
        % M with my alternative grazing formulation:
        %M = a*y(1)*y(2)-gz*(1-y(1)/(h+y(1)))*y(3)*y(1)+gamma*y(1)*(1-y(1)-y(2)) + alpha*(1-y(1)-y(2)); % macroalgae
        % M with van de Leemput grazing formulation:
        M = a*y(1)*y(2)-gz*y(3)*y(1)/(gz*h*y(1)+1)+gamma*y(1)*(1-y(1)-y(2)) + alpha*(1-y(1)-y(2)); % macroalgae
        C = r*(1-y(1)-y(2))*y(2)-d*y(2)-a*y(1)*y(2) + beta*(1-y(1)-y(2)); % coral
        H = rH*y(3)-dH*y(3)*y(3)-f*y(3); % herbivores
        %H = rH*y(3)*(1-y(3)/k); % herbivores
        s = [M,C,H]';
    end
%Initial conditions
    function [y0] = pdeic(xi)
        k = (rH-f)/dH; % initial herbivore abundance at each location
        
        % low coral
        if icchoice == 1
        y0 = [0.8,0.05,k]'; % if not random, specify initial M and C cover
        end

        % high coral
        if icchoice == 2
        y0 = [0.05,0.8,k]'; 
        end
        
        %RANDOMIZED
        if icchoice == 3
       % Mi = .5*rand;
        %y0 = [Mi,.8-Mi,k*(1+(2*rand-1))]'; %randomize macroalgae and fish
       % y0 = [Mi,0.8-Mi,k]'; %randomize macroalgae only
        % * rand generates random number between 0 and 1, multiple this by
        % 0.5 so max initial M is 0.5, say 0.2 of habitat is initially
        % turf, so initial C is 0.8-M. Initialize fish at carrying capacity
        % everywhere

         % update: make initial values and magnitude of randomness function
        % arguments
        Mi0 = Mhigh-Mhigh*rand*rnsize; % initial invul macroalgal cover
        C0 = Chigh-Chigh*rand*rnsize; % initial coral cover
        y0 = [Mi0,C0,k]'; % Minv, C, H, Mvuln

        % make sure Mhigh and Chigh don't sum to greater than 1 when using
        % random
        
        end
        
       % step function
        if icchoice == 4
        if ismember(xi, initC) ==1 % if xi is in initC
            %y0 = [0.05, 0.8, k]'; % high C
            y0 = [Mlow, Chigh, k]'; % high C
        else
            %y0 = [0.8, 0.05, k]'; % high M
            y0 = [Mhigh, Clow, k]'; % high M
        end
        end

        % sine wave
        if icchoice ==5

            C0i = ampC0*sin(period0*xi) + Chigh;
            M0i = ampM0*sin(period0*xi-pi) + Mhigh;
            y0 = [M0i, C0i, k]'; % need to transpose!!

        end
        
    end


%Dirichlet boundary conditions: biomass density goes to zero at the habitat
%edges; flux across the boundaries is allowed.

%Neumann (reflecting) boundary conditions: biomass density is non-zero at the habitat
%edges; flux across the boundaries is not allowed.
    function [pl,ql,pr,qr] = pdebc(xl,yl,xr,yr,t)
        if diric == 1 % Dirichlet
            pl = [yl(1),yl(2),yl(3)]';
            ql = [0,0,0]';
            pr = [yr(1),yr(2),yr(3)]';
            qr = [0,0,0]';
        else % Neumann
            pl = [0,0,0]';
            ql = [1,1,1]';
            pr = [0,0,0]';
            qr = [1,1,1]';
        end
    end

end

