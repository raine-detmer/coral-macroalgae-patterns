% README: function for simulating the PDE model used in the main text
% (local dynamics based on Briggs et al. 2018 model, with the addition of
% dynamic herbivores that have external recruitment)


function[sol] = BriggsHrPDEextH(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, ...
    dH, f, diff,taxisM,taxisC, taxisT, diric,x,t,initC,Clow, Chigh, Mlow, Mhigh, rnsize, ...
    ampC0, ampM0, period0,icchoice, phiH) % function of model parameters, boundary conditions, x, t, initial conditions


% Find solution by simulating with PDE solver
sol = pdepe(0,@pdefcn,@pdeic,@pdebc,x,t);
%Each row of sol represents a timepoint; each column a value of x
%Sol's third dimension follows assigment in the pdepe function.


%System of PDEs describing population dynamics
    function [c,g,s] = pdefcn(xi,t,y,DyDx)
        c = [1,1,1,1]'; % c= coefficient of the terms that get differentiated with respect to time
        g = [diff(1)*DyDx(1),diff(2)*DyDx(2),diff(3)*DyDx(3)+taxisM*y(3)*(DyDx(1)+ DyDx(4))+taxisC*y(3)*DyDx(2) + taxisT*y(3)*(-DyDx(1)-DyDx(2)-DyDx(4)), diff(4)*DyDx(4)]'; % things that get differentiated w/r/t space
        %third component describes changes in fish over space: have
        %diffusion, could have taxis with respect to C, M, and/or T
        %dMi/dt = omega*Mv+gTI*(1-Mi-Mv-C)*Mi+gamma*gTI*Mi*C-di*H*Mi; % invulnerable macroalgae
        Mi = omega*y(4)+gTI*(1-y(1)-y(4)-y(2))*y(1)+gamma*gTI*y(1)*y(2)-di*y(3)*y(1);
        %dC/dt = phiC*(1-Mi-Mv-C)+gTC*(1-Mi-Mv-C)*C -gamma*gTI*Mi*C-dC*C; % coral
        C = phiC*(1-y(1)-y(4)-y(2))+gTC*(1-y(1)-y(4)-y(2))*y(2) -gamma*gTI*y(1)*y(2)-dC*y(2); % coral
        % dH/d = phiH - rH*H - dH*H^2 - f*H % herbivores
        H = phiH + rH*y(3)-dH*y(3)*y(3)-f*y(3); 
        %Mv = phiM*(1-Mi-Mv-C)+rM*(1-Mi-Mv-C)*Mi+gTV*(1-Mi-Mv-C)*Mv-dv*H*Mv-omega*Mv;
        Mv = phiM*(1-y(1)-y(4)-y(2))+rM*(1-y(1)-y(4)-y(2))*y(1)+gTV*(1-y(1)-y(4)-y(2))*y(4)-dv*y(3)*y(4)-omega*y(4);
        s = [Mi,C,H, Mv]';
    end
%Initial conditions
    function [y0] = pdeic(xi)
    
        k = ((rH-f) + sqrt((rH-f)^2 + 4*dH*phiH))/(2*dH); % initial herbivore abundance at each location
        
        % low coral
        if icchoice == 1
        % y0 = vector of initial conditions, 1st element = Mi, 2nd = coral,
        % 3rd = herbivores, 4th = Mv
        y0 = [0.8,0.05,k, 0.05]'; % high uniform initial M cover, low uniform initial C cover
        end

        % high coral
        if icchoice == 2
        y0 = [0.05,0.8,k, 0.05]'; % low uniform initial M cover, high uniform initial C cover
        end
        
        %RANDOMIZED
        if icchoice == 3
        
        % calculate random cover
        Mi0 = Mhigh-Mhigh*rand*rnsize; % initial invul macroalgal cover
        C0 = Chigh-Chigh*rand*rnsize; % initial coral cover
        % assume that Mi0 is total macroalgae, and 95% of this is
        % invuln and 5% is vulnerable
        y0 = [Mi0*0.95,C0,k, Mi0*0.05]'; % Minv, C, H, Mvuln

        % NOTE: make sure Mhigh and Chigh don't sum to greater than 1 when using
        % random initial conditions

        end
        
       
        if icchoice == 4 % Specific step wise distribution
        if ismember(xi, initC) ==1 % if xi is in initC
            
            y0 = [Mlow*0.95, Chigh, k, Mlow*0.05]'; % low M and high C
        else
            y0 = [Mhigh*0.95, Clow, k, Mhigh*0.05]'; % high M and low C
        end
        end

        
        if icchoice ==5 % sine wave

            C0i = ampC0*sin(period0*xi) + Chigh; % coral cover is a sinusoidal function of spatial location xi
            M0i = ampM0*sin(period0*xi-pi) + Mhigh; % macroalgal cover is a sinusoidal function of spatial location xi
            y0 = [M0i*0.95, C0i, k, M0i*0.05]'; % need to transpose (') to make this a column vector; or could also separate with semicolons instead of commas

            % NOTE: amp0/2 + C0high + M0high-amp0/2 (high point + low
            % point) need to add up to <=1, so need C0high + M0high <=1
            

        end
        
    end


% function for boundary conditions
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

