% README: adding option for step-wise herbivore initial conditions to the
% Briggs PDE model


function[sol] = BriggsHrPDEHIC(phiC, gTC, gamma, gTI, dC, phiM, rM, gTV, dv, omega,di, rH, ...
    dH, f, diff,taxisM,taxisC, taxisT, diric,x,t,initC,Clow, Chigh, Mlow, Mhigh, rnsize, ...
    ampC0, ampM0, period0,icchoice, icchoiceH, initH)

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

        if icchoiceH ==1 % homogeneous fish
            H0 = k;
        end

        if icchoiceH == 2 % step-wise fish
            if ismember(xi, initH) ==1 % if xi is in initH
            H0 = k + 0.5*k; % high H
            else
            H0 = k - 0.5*k; % low H
            end

        end
        
        % low coral
        if icchoice == 1
        y0 = [0.8,0.05,H0, 0.05]'; % if not random, specify initial M and C cover
        end

        % high coral
        if icchoice == 2
        y0 = [0.05,0.8,H0, 0.05]'; 
        end
        
        %RANDOMIZED
        if icchoice == 3
       
        % update: make initial values and magnitude of randomness function
        % arguments
        Mi0 = Mhigh-Mhigh*rand*rnsize; % initial invul macroalgal cover
        C0 = Chigh-Chigh*rand*rnsize; % initial coral cover
        % then say remaining cover that isn't Mi or C is 50% vuln M and 50%
        % turf
        %y0 = [Mi0,C0,k, (1-Mi0-C0)*0.5]'; % Minv, C, H, Mvuln
        % update: say that Mi is total macroalgae, and 95% of this is
        % invuln
        y0 = [Mi0*0.95,C0,H0, Mi0*0.05]'; % Minv, C, H, Mvuln

        % make sure Mhigh and Chigh don't sum to greater than 1 when using
        % random

        end
        
       % Specific step wise distribution
        if icchoice == 4
        if ismember(xi, initC) ==1 % if xi is in initC
            %y0 = [0.05, 0.8, k]'; % high C
            %y0 = [Mlow, Chigh, k, 0]'; % high C
            y0 = [Mlow*0.95, Chigh, H0, Mlow*0.05]'; % high C
        else
            %y0 = [0.8, 0.05, k]'; % high M
            %y0 = [Mhigh, Clow, k, 0]'; % high M
            y0 = [Mhigh*0.95, Clow, H0, Mhigh*0.05]'; % high M
        end
        end

        % sine wave
        if icchoice ==5

            %xpos = xi;
            %C0i = ampC0*sin(period0*xpos) + Chigh;
            %M0i = ampM0*sin(period0*xpos-pi) + Mhigh;

            C0i = ampC0*sin(period0*xi) + Chigh;
            M0i = ampM0*sin(period0*xi-pi) + Mhigh;
            y0 = [M0i*0.95, C0i, H0, M0i*0.05]'; % need to transpose!!

            
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

