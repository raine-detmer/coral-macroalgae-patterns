% README: Briggs PDE model but with two herbivore populations that have separate
% taxis and diffusion and fishing pressures (used to make Figure 4)

function[sol] = Briggs2HrPDE(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, ...
    dH, f1,f2, diff,taxisM1,taxisC1, taxisT1,taxisM2,taxisC2, taxisT2, diric,x,t,initC,Clow, Chigh, Mlow, Mhigh, rnsize, ...
    ampC0, ampM0, period0,icchoice, phiH)


% Find solution by simulating with PDE solver
sol = pdepe(0,@pdefcn,@pdeic,@pdebc,x,t);
%Each row of sol represents a timepoint; each column a value of x
%Sol's third dimension follows assigment in the pdepe function.


%System of PDEs describing population dynamics
    function [c,g,s] = pdefcn(xi,t,y,DyDx)
        c = [1,1,1,1,1]'; % c= coefficient of the terms that get differentiated with respect to time
        g = [diff(1)*DyDx(1),diff(2)*DyDx(2),diff(3)*DyDx(3)+taxisM1*y(3)*(DyDx(1)+ DyDx(4)) + taxisC1*y(3)*DyDx(2) + taxisT1*y(3)*(-DyDx(1)-DyDx(2)-DyDx(4)), diff(4)*DyDx(4), diff(5)*DyDx(5)+taxisM2*y(5)*(DyDx(1)+ DyDx(4)) + taxisC2*y(5)*DyDx(2) + taxisT2*y(5)*(-DyDx(1)-DyDx(2)-DyDx(4))]'; % things that get differentiated w/r/t space
        %third component describes changes in fish over space: have diffusion, could have taxis with respect to M and C
        %Mi = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi; % invulnerable macroalgae
        Mi = omega*y(4)+gTI*(1-y(1)-y(4)-y(2))*y(1)+gamma*gTI*y(1)*y(2)-di*(y(3)+y(5))*y(1);
        %C = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C; % coral
        C = phiC*(1-y(1)-y(4)-y(2))+gTC*(1-y(1)-y(4)-y(2))*y(2) -gamma*gTI*y(1)*y(2)-dC*y(2); % coral
        H1 = 0.5*phiH + rH*y(3)-2*dH*y(3)*y(3)-f1*y(3); % first herbivore population
        H2 = 0.5*phiH + rH*y(5)-2*dH*y(5)*y(5)-f2*y(5); % second herbivore population
        Mv = phiM*(1-y(1)-y(4)-y(2))+rM*(1-y(1)-y(4)-y(2))*y(1)+gTV*(1-y(1)-y(4)-y(2))*y(4)-dv*(y(3)+y(5))*y(4)-omega*y(4);
        s = [Mi,C,H1, Mv, H2]';
    end
%Initial conditions
    function [y0] = pdeic(xi)
        % y0 = vector of initial conditions, 1st element = Mi, 2nd = coral,
        % 3rd = first herbivore population, 4th = Mv, 5th = second
        % herbivore population
        
        k1 = ((rH-f1) + sqrt((rH-f1)^2 + 4*2*dH*0.5*phiH))/(2*2*dH); % initial herbivore abundance at each location
        k2 = ((rH-f2) + sqrt((rH-f2)^2 + 4*2*dH*0.5*phiH))/(2*2*dH); % initial herbivore abundance at each location
        
        
        % low coral
        if icchoice == 1
        y0 = [0.8,0.05,k1, 0.05, k2]'; % uniform high M/low C
        end

        % high coral
        if icchoice == 2
        y0 = [0.05,0.8,k1, 0.05, k2]'; % uniform low M/high C
        end
        
        %RANDOMIZED
        if icchoice == 3

        
        % calculate random cover
        Mi0 = Mhigh-Mhigh*rand*rnsize; % initial invul macroalgal cover
        C0 = Chigh-Chigh*rand*rnsize; % initial coral cover
        % then say that Mi0 is total macroalgae, and 95% of this is
        % invuln and 5% is vulnerable
        y0 = [Mi0*0.95,C0,k1, Mi0*0.05, k2]'; % Minv, C, H1, Mvuln, H2

        % make sure Mhigh and Chigh don't sum to greater than 1 when using
        % random

        end
        
       % Specific step wise distribution
        if icchoice == 4
        if ismember(xi, initC) ==1 % if xi is in initC
           
            y0 = [Mlow*0.95, Chigh, k1, Mlow*0.05, k2]'; % high C
        else
            
            y0 = [Mhigh*0.95, Clow, k1, Mhigh*0.05, k2]'; % high M
        end
        end

        % sine wave
        if icchoice ==5


            C0i = ampC0*sin(period0*xi) + Chigh;
            M0i = ampM0*sin(period0*xi-pi) + Mhigh;
            y0 = [M0i*0.95, C0i, k1, M0i*0.05, k2]'; % need to transpose!!

         
        end
        
    end


%Dirichlet boundary conditions: biomass density goes to zero at the habitat
%edges; flux across the boundaries is allowed.

%Neumann (reflecting) boundary conditions: biomass density is non-zero at the habitat
%edges; flux across the boundaries is not allowed.
    function [pl,ql,pr,qr] = pdebc(xl,yl,xr,yr,t)
        if diric == 1 % Dirichlet
            pl = [yl(1),yl(2),yl(3), yl(4), yl(5)]';
            ql = [0,0,0,0,0]';
            pr = [yr(1),yr(2),yr(3), yr(4), yr(5)]';
            qr = [0,0,0,0,0]';
        else % Neumann
            pl = [0,0,0,0,0]';
            ql = [1,1,1,1,1]';
            pr = [0,0,0,0,0]';
            qr = [1,1,1,1,1]';
        end
    end

end

