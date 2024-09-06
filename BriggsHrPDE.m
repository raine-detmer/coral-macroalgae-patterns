% README: turning the Briggs model into a PDE (and adding dynamic fish to it)

% tried to add option for sin initial conditions but it seems to require a
% different method for initial condition specification

function[sol] = BriggsHrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, ...
    dH, f, diff,taxisM,taxisC, taxisT, diric,x,t,initC,Clow, Chigh, Mlow, Mhigh, rnsize, ...
    ampC0, ampM0, period0,icchoice)

% initial conditions for sin wave case
%if icchoice ==5
%global Mi00; global C00; global H00; global Mv00;
%Mi00 = (ampM0*sin(period0*x-pi) + Mhigh)*0.95;
%C00 = ampC0*sin(period0*x) + Chigh;
%Mv00 = (ampM0*sin(period0*x-pi) + Mhigh)*0.05;
%end

% Find solution by simulating with PDE solver
sol = pdepe(0,@pdefcn,@pdeic,@pdebc,x,t);
%Each row of sol represents a timepoint; each column a value of x
%Sol's third dimension follows assigment in the pdepe function.


%System of PDEs describing population dynamics
    function [c,g,s] = pdefcn(xi,t,y,DyDx)
        c = [1,1,1,1]'; % c= coefficient of the terms that get differentiated with respect to time
        g = [diff(1)*DyDx(1),diff(2)*DyDx(2),diff(3)*DyDx(3)+taxisM*y(3)*(DyDx(1)+ DyDx(4))+taxisC*y(3)*DyDx(2) + taxisT*y(3)*(-DyDx(1)-DyDx(2)-DyDx(4)), diff(4)*DyDx(4)]'; % things that get differentiated w/r/t space
        %third component describes changes in fish over space: have diffusion, could have taxis with respect to M and C
        %Mi = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi; % invulnerable macroalgae
        Mi = omega*y(4)+gTI*(1-y(1)-y(4)-y(2))*y(1)+gamma*gTI*y(1)*y(2)-di*y(3)*y(1);
        %C = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C; % coral
        C = phiC*(1-y(1)-y(4)-y(2))+gTC*(1-y(1)-y(4)-y(2))*y(2) -gamma*gTI*y(1)*y(2)-dC*y(2); % coral
        H = rH*y(3)-dH*y(3)*y(3)-f*y(3); % herbivores
        %H = rH*y(3)*(1-y(3)/k); % herbivores
        %Mv = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv;
        Mv = phiM*(1-y(1)-y(4)-y(2))+rM*(1-y(1)-y(4)-y(2))*y(1)+gTV*(1-y(1)-y(4)-y(2))*y(4)-dv*y(3)*y(4)-omega*y(4);
        s = [Mi,C,H, Mv]';
    end
%Initial conditions
    function [y0] = pdeic(xi)
        k = (rH-f)/dH; % initial herbivore abundance at each location
        
        % low coral
        if icchoice == 1
        y0 = [0.8,0.05,k, 0.05]'; % if not random, specify initial M and C cover
        end

        % high coral
        if icchoice == 2
        y0 = [0.05,0.8,k, 0.05]'; 
        end
        
        %RANDOMIZED
        if icchoice == 3
        %Mi = .5*rand;
        %y0 = [Mi,.8-Mi,k*(1+(2*rand-1))]'; %randomize macroalgae and fish
        %y0 = [Mi,0.8-Mi,k, 0]'; %randomize macroalgae only
        % * rand generates random number between 0 and 1, multiple this by
        % 0.5 so max initial M is 0.5, say 0.2 of habitat is initially
        % turf, so initial C is 0.8-M. Initialize fish at carrying capacity
        % everywhere

        % update: let M initially be higher
       % Mi = .9*rand; % initial invul macroalgal cover
       % Cprop = rand; % proportion of remaining cover that is coral
        % then say remaining cover that isn't Mi or C is 50% vuln M and 50%
        % turf
       % y0 = [Mi,(1-Mi)*Cprop,k, (1-Mi-(1-Mi)*Cprop)*0.5]'; % Minv, C, H, Mvuln

        % update: make initial values and magnitude of randomness function
        % arguments
        Mi0 = Mhigh-Mhigh*rand*rnsize; % initial invul macroalgal cover
        C0 = Chigh-Chigh*rand*rnsize; % initial coral cover
        % then say remaining cover that isn't Mi or C is 50% vuln M and 50%
        % turf
        %y0 = [Mi0,C0,k, (1-Mi0-C0)*0.5]'; % Minv, C, H, Mvuln
        % update: say that Mi is total macroalgae, and 95% of this is
        % invuln
        y0 = [Mi0*0.95,C0,k, Mi0*0.05]'; % Minv, C, H, Mvuln

        % make sure Mhigh and Chigh don't sum to greater than 1 when using
        % random

        end
        
       % Specific step wise distribution
        if icchoice == 4
        if ismember(xi, initC) ==1 % if xi is in initC
            %y0 = [0.05, 0.8, k]'; % high C
            %y0 = [Mlow, Chigh, k, 0]'; % high C
            y0 = [Mlow*0.95, Chigh, k, Mlow*0.05]'; % high C
        else
            %y0 = [0.8, 0.05, k]'; % high M
            %y0 = [Mhigh, Clow, k, 0]'; % high M
            y0 = [Mhigh*0.95, Clow, k, Mhigh*0.05]'; % high M
        end
        end

        % sine wave
        if icchoice ==5

            %xpos = xi;
            %C0i = ampC0*sin(period0*xpos) + Chigh;
            %M0i = ampM0*sin(period0*xpos-pi) + Mhigh;

            C0i = ampC0*sin(period0*xi) + Chigh;
            M0i = ampM0*sin(period0*xi-pi) + Mhigh;
            y0 = [M0i*0.95, C0i, k, M0i*0.05]'; % need to transpose!!

            % NOTE: amp0/2 + C0high + M0high-amp0/2 (high point + low
            % point) need to add up to <=1, so need C0high + M0high <=1
            %global Mi00; global C00; global H00; global Mv00;
            %C0i = C00(x==xi);
            %M0i = Mi00(x==xi);
            %M0v = Mv00(x==xi);
            %H0 = k;

            %y0 = [M0i, C0i, H0, M0v]'; % REMEMBER NEED TRANSPOSE (') to make this a column, could also separate with semicolons instead of commas

        end
        
    end


%Dirichlet boundary conditions: biomass density goes to zero at the habitat
%edges; flux across the boundaries is allowed.

%Neumann (reflecting) boundary conditions: biomass density is non-zero at the habitat
%edges; flux across the boundaries is not allowed.
    function [pl,ql,pr,qr] = pdebc(xl,yl,xr,yr,t)
        if diric == 1 % Dirichlet
            pl = [yl(1),yl(2),yl(3), yl(4)]';
            ql = [0,0,0,0]';
            pr = [yr(1),yr(2),yr(3), yr(4)]';
            qr = [0,0,0,0]';
        else % Neumann
            pl = [0,0,0,0]';
            ql = [1,1,1,1]';
            pr = [0,0,0,0]';
            qr = [1,1,1,1]';
        end
    end

end

