function bssn
    global v
    global n_var
    global r_min
    global r_max
    global diss_on
    global eta
    global puncture
    global sigma
    global precollapse
    v = 1; % Eulerian condition (v=0) or Lagrangian condition (v=1)
    eta = 1; % eta parameter in the Gamma-driver condition
    n_var = 12; % number of functions to evolve
    r_min = 0;
    r_max = 5;
    diss_on = 0; % dissipation (Kreiss-Oliger) on (=1) or off (=0)
    sigma = 0.1;
    puncture = 1;
    precollapse = 1;
    h = 1/10; % spatial grid
    N=round((r_max-r_min)/h);
    r=(r_min+h/2:h:r_max-h/2);

    
    t_end = 2;
    plot_RHS(20,r_min,r_max)
    [U_coarse, t_size_coarse]=state_result(1/10, t_end);
    [U_med, t_size_med]=state_result(1/20, t_end);
    [U_fine, t_size_fine]=state_result(1/40, t_end);
    [U_finest, t_size_finest]=state_result(1/80, t_end);
    
    [alpha_coarse, beta_r_coarse, B_r_coarse, chi_coarse, g_rr_coarse,...
    g_thth_coarse, A_rr_coarse, K_coarse, Gamma_r_coarse, H_coarse, M_r_coarse, G_coarse]=var_t(U_coarse);
    [alpha_med, beta_r_med, B_r_med, chi_med, g_rr_med, g_thth_med,...
    A_rr_med, K_med, Gamma_r_med, H_med, M_r_med, G_med]=var_t(U_med);
    [alpha_fine, beta_r_fine, B_r_fine, chi_fine, g_rr_fine, g_thth_fine,...
    A_rr_fine, K_fine, Gamma_r_fine, H_fine, M_r_fine, G_fine]=var_t(U_fine) ;
    [alpha_finest, beta_r_finest, B_r_finest, chi_finest, g_rr_finest, g_thth_finest,...
    A_rr_finest, K_finest, Gamma_r_finest, H_finest, M_r_finest, G_finest]=var_t(U_finest) ;
    
    error_medcoarse = error(beta_r_med(1:2:2*N, t_size_med), beta_r_coarse(:, t_size_coarse));
    error_finemed = 2*error(beta_r_fine(1:4:N*4, t_size_fine),beta_r_med(1:2:2*N, t_size_med));
    error_finestfine = 4*error(beta_r_finest(1:8:N*8, t_size_finest),beta_r_fine(1:4:4*N, t_size_fine));
    plot(r, error_medcoarse, '-r');
    hold on
    plot(r, error_finemed, 'b--');
    plot(r, error_finestfine, 'g-.');
    xlabel('$r$','Interpreter','latex');
    ylabel('$\Delta \chi$','Interpreter','latex') ;
    l=legend('$(\beta^r_{M/50}-\beta^r_{M/100})$', '$2(\beta^r_{M/100}-\beta^r_{M/200})$', '$4(\beta^r_{M/200}-\beta^r_{M/400})$');
    set(l, 'Interpreter', 'latex');
    set(l, 'FontSize', 14);
    %plotting
%     plot(r,chi_coarse(:, t_size_coarse), 'r')
%     hold on
%     plot(r,chi_med(:, t_size_med), 'b--')  
%     plot(r,chi_fine(:, t_size_fine), 'g-.')  
%     xlabel('$r$','Interpreter','latex')
%     ylabel('$\chi$','Interpreter','latex') 
%     legend('\eta = 0','\eta = 1','\eta = 2');
%     xlim([0,10]);
    
%     
%     %solving step
%     % initial condition
%      tspan = [0 t_end];
%      y0 = initial_cond(h,r_min,r_max);
%      [t,y] = ode45(@(t,y) dydt(t,y,h,N),tspan,y0);
%      %getting the output size
%      [t_size,y_size] = size(y)
%      % reshaping the array so it's nicer to unpackage
%      U = reshape(y, t_size, [], n_var);
%      % unpackaging into var(time, space)
%      [alpha, beta_r, B_r, chi, g_rr, g_thth, A_rr, K, Gamma_r, H, M_r, G] = var_t(U);
    
    

    %plotting results
%     for iter=1:t_size
%         [charac1, charac2]=CHARAC(chi(:,iter), g_rr(:,iter), g_thth(:,iter),...
%             A_rr(:,iter),K(:,iter), Gamma_r(:,iter), h, N);
%         plot(r, charac1)
%         ylim([-5 5])
%         pause(0.001)
%     end
%     xlabel('r')
%     ylabel('charac2')
% 
% Horizon computation

[horsize_r horsize_t] = size(var_t(U));
hor=zeros(horsize_t,1);
time=zeros(horsize_t,1);
    for iter=1:horsize_t
        time(iter) = t_end/horsize_t*iter;
        hor(iter) = Hor(g_rr(:,iter),g_thth(:,iter),K(:,iter),A_rr(:,iter),chi(:,iter),N,r,h);
    end
    
     plot(time,hor)
      
% 
%     %unpackaging
%     charac1=Q(:,:,1);
%     charac2=Q(:,:,2);
%     %plotting
%     for iter=1:t_size
%         plot(r,alpha(iter,:))
%         pause(0.001)
%     end
end

function [U, t_size]=state_result(h, t_end)
    global n_var
    global r_min
    global r_max
    N=round((r_max-r_min)/h);
    r=(r_min+h/2:h:r_max-h/2);
    
    % time at which we want to end the simulation
    % time to solve the equations
    tspan = [0 t_end];
    
    y0 = initial_cond(h,r_min,r_max);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [t,y] = ode45(@(t,y) dydt(t,y,h,N),tspan,y0, options);

    %getting the output size
    [t_size,y_size] = size(y);
    % reshaping the array so it's nicer to unpackage
    U = reshape(y, t_size, [], n_var);
end


function U=dydt(t, v_old, h, N)
    global n_var
    global r_min
    global r_max
    global v
    global diss_on
    global eta
    global puncture
    global sigma
    global precollapse
    U = zeros(n_var*N,1);
    % unpacking our previous state
    v_old=reshape(v_old,[],n_var);
    r=(r_min+h/2:h:r_max-h/2);
    [alpha, beta_r, B_r, chi, g_rr, g_thth, A_rr, K, Gamma_r, H, M_r, G] = unpackage(v_old);

    
    % radial derivatives
    
    % precollapsed alpha or unity
    if precollapse == 1
        alpha_p = f_prime(alpha,h,alpha(2),alpha(1),power(1+0.5/(r_max+h/2),-2),power(1+1/2/(r_max+3*h/2),-2),N);
        alpha_pp = f_pprime(alpha,h,alpha(2),alpha(1),power(1+0.5/(r_max+h/2),-2),power(1+1/2/(r_max+3*h/2),-2),N);
    else
        alpha_p = f_prime(alpha,h,alpha(2),alpha(1),1,1,N);
        alpha_pp = f_pprime(alpha,h,alpha(2),alpha(1),1,1,N);
    end
    
    beta_r_p = f_prime(beta_r,h,-beta_r(2),-beta_r(1),0,0,N);
    beta_r_pp = f_pprime(beta_r,h,-beta_r(2),-beta_r(1),0,0,N);
    B_r_p= f_prime(B_r,h,-B_r(2),-B_r(1),0,0,N);
    % appropriate boundary conditions for flat space
    if puncture == 0
        chi_p = f_prime(chi,h,chi(2),chi(1),1,1,N);
        chi_pp = f_pprime(chi,h,chi(2),chi(1),1,1,N);
    else
        chi_p = f_prime(chi,h,chi(2),chi(1),power(1 + 1/2/(r_max+h/2), -4),power(1 + 1/2/(r_max+3*h/2), -4),N);
        chi_pp = f_pprime(chi,h,chi(2),chi(1),power(1 + 1/2/(r_max+h/2), -4),power(1 + 1/2/(r_max+3*h/2),-4),N);
    end
    g_rr_p = f_prime(g_rr,h,g_rr(2),g_rr(1),1,1,N);
    g_rr_pp = f_pprime(g_rr,h,g_rr(2),g_rr(1),1,1,N);
    g_thth_p = f_prime(g_thth,h,g_thth(2), g_thth(1), (r_max+h/2)^2,(r_max+3*h/2)^2,N);
    g_thth_pp = f_pprime(g_thth,h,g_thth(2), g_thth(1),(r_max+h/2)^2,(r_max+3*h/2)^2,N);
    A_rr_p =  f_prime(A_rr,h,A_rr(2),A_rr(1),0,0,N);
    K_p = f_prime(K,h,K(2),K(1),0,0,N);
    Gamma_r_p = f_prime(Gamma_r,h,-Gamma_r(2),-Gamma_r(1),-2/(r_max+h/2),-2/(r_max+3*h/2),N);
    % The initial conditions below are what we used before for r_min = 1.
    % at r_min = 0 some stuff go wrong
    if r_min ~= 0
        g_thth_p = f_prime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        g_thth_pp = f_pprime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        Gamma_r_p = f_prime(Gamma_r,h,-2/(r_min-3*h/2),-2/(r_min-h/2),-2/(r_max+h/2),-2/(r_max+3*h/2) ,N);
    end


    % constraint evolution equations
    H_p = f_prime(H,h,H(2),H(1),0,0,N);
    M_r_p = f_prime(M_r,h,M_r(2),M_r(1),0,0,N);
    M_r_pp = f_pprime(M_r,h,M_r(2),M_r(1),0,0,N);
    G_p = f_prime(G,h,G(2),G(1),0,0,N);
    G_pp = f_pprime(G,h,G(2),G(1),0,0,N);

    % time derivatives for each state variable
    if precollapse == 0
        alpha_t = beta_r.*alpha_p-2*alpha.*K...
                  +diss_on*KreissOliger(alpha,sigma,h,alpha(3),alpha(2),alpha(1),1,1,1,N); % eqn 1
    else
        alpha_t = beta_r.*alpha_p-2*alpha.*K...
                  +diss_on*KreissOliger(alpha,sigma,h,alpha(3),alpha(2),alpha(1),...
                  power(1+1/2./(r_max+h/2),-2),power(1+1/2/(r_max+3*h/2),-2),power(1+1/2/(r_max+5*h/2),-2),N); % eqn 1
    end
    beta_r_t = 3/4*B_r+beta_r.*beta_r_p...
               +diss_on*KreissOliger(beta_r,sigma,h,-beta_r(3),-beta_r(2),-beta_r(1),0,0,0,N); % note this is only 2a), what is B_r?
    chi_t = 2/3*K.*alpha.*chi - v.*beta_r.*g_rr_p.*chi./(3*g_rr)...
            -2*v.*beta_r.*g_thth_p.*chi./(3*g_thth)-2/3*v.*beta_r_p.*chi...
            +beta_r.*chi_p...
            +diss_on*KreissOliger(chi,sigma,h,chi(3),chi(2),chi(1),...
            power(1 + 1/2/(r_max+h/2), -4),power(1 + 1/2/(r_max+3*h/2), -4),power(1 + 1/2/(r_max+5*h/2), -4),N);
    g_rr_t = -2*A_rr.*alpha-v.*beta_r.*g_rr_p./3+beta_r.*g_rr_p...
             -2*g_rr.*v.*beta_r.*g_thth_p./(3*g_thth)+2*g_rr.*beta_r_p...
             -2/3*g_rr.*v.*beta_r_p...
             +diss_on*KreissOliger(g_rr,sigma,h,g_rr(3),g_rr(2),g_rr(1),1,1,1,N);
    g_thth_t = A_rr.*g_thth.*alpha./g_rr-g_thth.*v.*beta_r.*g_rr_p./(3*g_rr)...
               -2/3*v.*beta_r.*g_thth_p+beta_r.*g_thth_p...
               -2/3*g_thth.*v.*beta_r_p...
               +diss_on*KreissOliger(g_thth,sigma,h,g_thth(3),g_thth(2),g_thth(1),(r_max+h/2)^2,(r_max+3*h/2)^2,(r_max+5*h/2)^2,N);
    A_rr_t = -2*alpha.*A_rr.^2./g_rr+K.*alpha.*A_rr-v.*beta_r.*g_rr_p.*A_rr./(3*g_rr)...
             -2*v.*beta_r.*g_thth_p.*A_rr./(3*g_thth)-2/3*v.*beta_r_p.*A_rr...
             +2*beta_r_p.*A_rr+2*alpha.*chi.*g_rr_p.^2./(3*g_rr.^2)... % end of first line
             -alpha.*chi.*g_thth_p.^2./(3*g_thth.^2)-alpha.*chi_p.^2./(6*chi)...
             -2*g_rr.*alpha.*chi./(3*g_thth)+beta_r.*A_rr_p+2/3*g_rr.*alpha.*chi.*Gamma_r_p...
             -alpha.*chi.*g_rr_p.*g_thth_p./(2*g_rr.*g_thth)+chi.*g_rr_p.*alpha_p./(3*g_rr)... % end of second line
             +chi.*g_thth_p.*alpha_p./(3*g_thth)-alpha.*g_rr_p.*chi_p./(6*g_rr)...
             -alpha.*g_thth_p.*chi_p./(6*g_thth)-2/3*alpha_p.*chi_p...
             -alpha.*chi.*g_rr_pp./(3*g_rr)+alpha.*chi.*g_thth_pp./(3*g_thth)...
             -2/3*chi.*alpha_pp+alpha.*chi_pp/3....
             +diss_on*KreissOliger(A_rr,sigma,h,A_rr(3),A_rr(2),A_rr(1),0,0,0,N);
    K_t = 3*alpha.*A_rr.^2./(2*g_rr.^2)+K.^2.*alpha./3+beta_r.*K_p...
          +chi.*g_rr_p.*alpha_p./(2*g_rr.^2)-chi.*g_thth_p.*alpha_p./(g_rr.*g_thth)...
          +alpha_p.*chi_p./(2*g_rr)-chi.*alpha_pp./g_rr...
          +diss_on*KreissOliger(K,0.05,h,K(3),K(2),K(1),0,0,0,N);
    Gamma_r_t = -v.*beta_r.*g_thth_p.^2./(g_rr.*g_thth.^2)...
                +A_rr.*alpha.*g_thth_p./(g_rr.^2.*g_thth)...
               -v.*beta_r_p.*g_thth_p./(3*g_rr.*g_thth)+beta_r_p.*g_thth_p./(g_rr.*g_thth)...
               +beta_r.*Gamma_r_p+A_rr.*alpha.*g_rr_p./(g_rr.^3)...
               -4/3*alpha.*K_p./g_rr-2*A_rr.*alpha_p./g_rr.^2.... % end of first line
               +v.*g_rr_p.*beta_r_p./(2*g_rr.^2)-g_rr_p.*beta_r_p./(2*g_rr.^2)...
               -3*A_rr.*alpha.*chi_p./(g_rr.^2.*chi)+v.*beta_r.*g_rr_pp./(6*g_rr.^2)...
               +v.*beta_r.*g_thth_pp./(3*g_rr.*g_thth)+v.*beta_r_pp./(3*g_rr)...
               +beta_r_pp./g_rr...
               +diss_on*KreissOliger(Gamma_r,sigma,h,-Gamma_r(3),-Gamma_r(2),-Gamma_r(1),-2/(r_max+h/2),-2/(r_max+3*h/2),-2/(r_max+5*h/2),N);
    % This has to be here since it uses Gamma_r_t
    B_r_t = Gamma_r_t- eta.*B_r+beta_r.*B_r_p - beta_r.*Gamma_r_p ...
           +diss_on*KreissOliger(K,sigma,h,-B_r(3),-B_r(2),-B_r(1),0,0,0,N);
    % This is from 2b). Not sure what the parameter \eta should be but I've set it to 0 for now

    % Constraint evolution system
    H_t = beta_r.*H + 2/3.*alpha.*K.*H-2.*alpha.*A_rr.*chi./g_rr.*G_p-2.*alpha./g_rr.*chi.*M_r_p...
        + (alpha.*chi_p./g_rr+alpha.*chi.*g_rr_p./g_rr./g_rr-4.*alpha_p.*chi./g_rr...
        -2.*alpha.*chi.*g_thth_p./g_thth).*M_r;
    M_r_t = beta_r.*M_r_p + beta_r_p.*M_r - alpha.*K.*M_r - alpha_p./3.*H + alpha./3.*H_p + 2/3.*alpha.*chi.*G_pp...
        +(2/3.*alpha_p.*chi - alpha.*chi_p./3 + alpha.*chi.*g_thth_p./g_thth).*G_p;
    G_t = beta_r.*G_p+ 2.*alpha./g_rr.*M_r;

    U = [alpha_t; beta_r_t; B_r_t; chi_t; g_rr_t; g_thth_t; A_rr_t; K_t; Gamma_r_t; H_t; M_r_t; G_t];
    

end

% this requires the fine grid points to be at the same points
function y=error(fine, coarse)
    y=coarse-fine;
end

function [alpha, beta_r, B_r, chi, g_rr, g_thth, A_rr, K, Gamma_r, H, M_r, G]=var_t(U)
    alpha=U(:,:,1).';
    beta_r=U(:,:,2).';
    B_r = U(:,:,3).';
    chi=U(:,:,4).';
	g_rr=U(:,:,5).';
	g_thth=U(:,:,6).';
	A_rr=U(:,:,7).';
	K=U(:,:,8).';
	Gamma_r=U(:,:,9).';
    H=U(:,:,10).';
    M_r=U(:,:,11).';
    G=U(:,:,12).';
end

function [charac1, charac2]=CHARAC(chi, g_rr, g_thth, A_rr, K, Gamma_r, h, N)
    global r_min
    global r_max

    chi_p = f_prime(chi,h,chi(2),chi(1),power(1 + 1/2/(r_max+h/2), -4),power(1 + 1/2/(r_max+3*h/2), -4),N);
    g_rr_p = f_prime(g_rr,h,g_rr(2),g_rr(1),1,1,N);
    g_thth_p = f_prime(g_thth,h,g_thth(2), g_thth(1),(r_max+h/2)^2,(r_max+3*h/2)^2,N);

    charac1=Gamma_r - 3/2.*A_rr./((g_rr.^3).*chi).^0.5+1/2*chi_p./(2*g_rr.*chi)-g_rr_p./(2.*g_rr.^2)...
        +g_thth_p./(2.*g_rr.*g_thth)+ K./(g_rr.*chi).^0.5;

    charac2=Gamma_r + 3/2.*A_rr./((g_rr.^3).*chi).^0.5+1/2*chi_p./(2*g_rr.*chi)-g_rr_p./(2.*g_rr.^2)...
        +g_thth_p./(2.*g_rr.*g_thth)- K./(g_rr.*chi).^0.5;
end

function v_old=initial_cond(h,r_min,r_max)
    % puncture == 1 means we have a puncture, otherwise we have flat space
    global n_var
    global puncture
    global precollapse
    M =1;
    N=round((r_max-r_min)/h);
    r=(r_min+h/2:h:r_max-h/2);
    v_old=zeros(N,n_var);
    % flat initial conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = ones(N,1);
    beta_r = zeros(N,1);
    B_r = zeros(N,1);
    chi = ones(N,1);
    g_rr = ones(N,1); 
    g_thth = (r.^2).';
    %g_thth = max(g_thth, 0.01);
    A_rr = zeros(N,1);
    K = zeros(N,1);
    Gamma_r = (-2./r).';
    %Gamma_r = min(-2/0.1, Gamma_r);
    % parameter s determines whether we use cap for chi or inv_r
    s = 0;
    if precollapse == 1
        alpha = power(1+M/2./(r.'),-2);
    end
    
    if puncture == 1
    	% punctured Schwarzchild BH ICs
        if s==1
 	        % tune r0 near 2 but not equal or above
         	chi = cap(1.8,r.',M,N);
        else
 		% tune r0 near 0 but not equal or below
        chi = power(1 + M./(2*r.'), -4);
 		%chi = power(1.+M/2.*inv_r(r.',.25,N),-4.);
        end
    end
    
    % radial derivatives
    if precollapse == 0
        alpha_p = f_prime(alpha,h,1,1,1,1,N);
        alpha_pp = f_pprime(alpha,h,1,1,1,1,N);
    else
        alpha_p = f_prime(alpha,h,0,0,power(1+M/2./(r_max+h/2),-2),power(1+M/2./(r_max+h/2),-2),N);
        alpha_pp = f_pprime(alpha,h,power(1+M/2./(r_max+h/2),-2),power(1+M/2./(r_max+h/2),-2),1,1,N);
    end
    beta_r_p = f_prime(beta_r,h,0,0,0,0,N);
    beta_r_pp = f_pprime(beta_r,h,0,0,0,0,N);
    B_r_p= f_prime(B_r,h,0,0,0,0,N);
    chi_p = f_prime(chi,h,1,1,1,1,N);
    chi_pp = f_pprime(chi,h,1,1,1,1,N);
    g_rr_p = f_prime(g_rr,h,1,1,1,1,N);
    g_rr_pp = f_pprime(g_rr,h,1,1,1,1,N);
    g_thth_p = f_prime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
    g_thth_pp = f_pprime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
    A_rr_p =  f_prime(A_rr,h,0,0,0,0,N);
    K_p = f_prime(K,h,0,0,0,0,N);
    Gamma_r_p = f_prime(Gamma_r,h, -2/(r_min-3*h/2), -2/(r_min-h/2),-2/(r_max+h/2),-2/(r_max+3*h/2),N);
    
    % Building constraints
    H = -3/2*A_rr.*A_rr./g_rr./g_rr + 2/3.*K.*K - 5/2.*chi_p.*chi_p./chi./g_rr...
            +2*chi_pp./g_rr + 2*chi./g_thth-2*chi.*g_thth_pp./g_rr./g_thth...
            +2*chi_p.*g_thth_p./g_rr./g_thth+chi.*g_rr_p.*g_thth_p./g_rr./g_rr./g_thth...
            -chi_p.*g_rr_p./g_rr./g_rr+chi.*g_thth_p.*g_thth_p./2./g_rr./g_thth./g_thth;
    M_r = A_rr_p./g_rr - 2/3.*K_p...
        - 3/2.*A_rr./g_rr.*(chi_p./chi - g_thth_p./g_thth+ g_rr_p./g_rr);
    G = -g_rr_p./2./g_rr./g_rr+Gamma_r+g_thth_p./g_thth./g_rr;
        
    % initializing the initial state with ICs
    v_old(:,1) = alpha;
    v_old(:,2) = beta_r;
    v_old(:,3) = B_r;
    v_old(:,4) = chi;
    v_old(:,5) = g_rr;
    v_old(:,6) = g_thth;
    v_old(:,7) = A_rr;
    v_old(:,8) = K;
    v_old(:,9) = Gamma_r;
    v_old(:,10) = H;
    v_old(:,11) = M_r;
    v_old(:,12) = G;
end

function plot_RHS(n_iter,r_min,r_max)
     L=zeros(1,10);
     H=zeros(1,10);
     for j=1:n_iter
         h=0.5^j;
         N = round((r_max-r_min)/h);
         v_old = initial_cond(h,r_min,r_max);
         U = dydt(1, v_old, h, N);
         U=reshape(U.',[],12);
         Uj=U(:,7);
         L(1,j)=sqrt(h)*sqrt(sum((abs(Uj).^2)));%sqrt(mean((Uj).^2));
         H(1,j)=h;
     end
     %plot loglog plot
     loglog(H,L)
     xlabel('h')
     ylabel('RHS')
     %plot slope of loglog plot    
     %plot(log10(H(1,1:19)),(log10(L(1,2:20)) - log10(L(1,1:19)))./( log10(H(1,2:20)) - log10(H(1,1:19)) ) )
end

% this function unpackages the *compressed* state vector v
% N is the length of the grid
function [alpha, beta_r, B_r, chi, g_rr, g_thth, A_rr, K, Gamma_r, H, M_r, G]=unpackage(v)
    alpha = v(:,1);
    beta_r = v(:,2);
    B_r = v(:,3);
    chi=v(:,4);
    g_rr=v(:,5);
    g_thth=v(:,6);
    A_rr=v(:,7);
    K=v(:,8);
    Gamma_r=v(:,9);
    H=v(:,10);
    M_r=v(:,11);
    G=v(:,12);
end

% This function returns f'(x) where f is one of the state variables
function y=f_prime(f,h,a,b,c,d,N)
    y=zeros(N,1);
    % These are for the two-level boundary conditions at each end
    y(1) = (-f(3) + 8*f(2) - 8*b + a)./(12*h);
    y(2) = (-f(4) + 8*f(3) - 8*f(1) + b)./(12*h);
    y(N-1) = (-c + 8*f(N) - 8*f(N-2) + f(N-3))./(12*h);
    y(N) = (-d + 8*c - 8*f(N-1) + f(N-2))./(12*h);
    % Computing the middle parts
    y(3:N-2) = (-f(5:N) + 8*f(4:N-1) - 8*f(2:N-3) + f(1:N-4))./(12*h);
end

% This function returns f''(x) where f is one of the state variables
function y=f_pprime(f,h,a,b,c,d,N)
    y=zeros(N,1);
    % These are for the two-level boundary conditions at each end
    y(1) = (-f(3) + 16*f(2) - 30*f(1) + 16*b - a)./(12*h^2);
    y(2) = (-f(4) + 16*f(3) - 30*f(2) + 16*f(1) - b)./(12*h^2);
    y(N-1) = (-c + 16*f(N) - 30*f(N-1) + 16*f(N-2) - f(N-3))./(12*h^2);
    y(N) = (-d + 16*c - 30*f(N) + 16*f(N-1) - f(N-2))./(12*h^2);
    % Computing the middle parts
    y(3:N-2) = (-f(5:N) + 16*f(4:N-1) - 30*f(3:N-2) + 16*f(2:N-3) - f(1:N-4))./(12*h^2);
end

% capping function for chi at t=0
function chi0=cap(r0,r,M,N)
    f=zeros(N,1);
    a=(2-r0)./2./r0;
    f= r0.*r + (1./(1+M/2./r)-r0.*r).*heaviside(r-a);
    chi0 = power(f,4).';
end

% Kreiss-Oliger dissipation term
% since our finite differencing is fourth order we have 6th order third
% derivative, which requires 4 points to each side
function v=KreissOliger(f,sigma,h,a,b,c,x,y,z,N)
    v=zeros(N,1);
    prefactor = (sigma*(-1)^6*h^5/2^6);
    v(1) = prefactor*(f(4)-6*f(3)+15*f(2)-20*f(1)+15*c-6*b+1*a)./h^6;
    v(2) = prefactor*(f(5)-6*f(4)+15*f(3)-20*f(2)+15*f(1)-6*c+b)./h^6;
    v(3) = prefactor*(f(6)-6*f(5)+15*f(4)-20*f(3)+15*f(2)-6*f(1)+c)./h^6;
    v(N-2) = prefactor*(x-6*f(N)+15*f(N-1)-20*f(N-2)+15*f(N-3)-6*f(N-4)+f(N-5))./h^6;
    v(N-1) = prefactor*(y-6*x+15*f(N)-20*f(N-1)+15*f(N-2)-6*f(N-3)+f(N-4))./h^6;
    v(N) = prefactor*(z-6*y+15*x-20*f(N)+15*f(N-1)-6*f(N-2)+f(N-3))./h^6;
    % Computing the middle parts
    v(4:N-3) = prefactor*(f(7:N)-6*f(6:N-1)+15*f(5:N-2)-20*f(4:N-3)...
        +15*f(3:N-4)-6*f(2:N-5)+f(1:N-6))./h^6;
end

% Horizon and expansion 
function [hum] = Hor(g_rr,g_thth,K,A_rr,chi,N,r,h)
    global r_max
    hum=0.;
    j=1;
    chi_inv = 1./chi;
    g_thth_p = f_prime(g_thth,h,g_thth(2), g_thth(1),(r_max+h/2)^2,(r_max+3*h/2)^2,N);
    % reconstruct thth-entry of full extrinsic curvature
    if j==0
        % without dividing by g_thth prior to Theta 
        K_thth = chi_inv.*(-A_rr.*g_thth./g_rr + 1/3.*g_thth.*K);
          % Leo: is it safe to invert g_rr/g_thth and take sqrt of it?
        % expansion Theta (Eq. (3.3) in the apparent horizon finding paper)
        Theta = g_thth_p./g_thth./sqrt(g_rr)-2.*K_thth./g_thth;
    else
        % dividing g_thth prior to Theta
        K_thth = chi_inv.*(-A_rr./g_rr + 1/3.*K);
        % expansion Theta (Eq. (3.3) in the apparent horizon finding paper)
        Theta = g_thth_p./g_thth./sqrt(g_rr)-2.*K_thth;
    end

    % Simple root-finder for horizon
    % threshold "0" up to machine precision
    hold = 1.e-15;
    
    for i=1:size(Theta)-1
        s = Theta(i)*Theta(i+1);
        % negative successive multiples => Theta crosses 0 at index i
        if s < hold
            % average between radii i and i+1
            hum = 1/2*(r(i)+r(i+1));
            % terminate loop
            break;
        end
    end
end
