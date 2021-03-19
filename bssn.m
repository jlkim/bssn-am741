function bssn
    % call the other function here
    % example:
    %[A,B] = compute_bssn(1, 2, 0.000001)
    % This will return a matrix A of time derivatives of state variables
    % and a matrix B of time derivatives of constraints
    %defining L(containing RMS of A_rr_t) and H(the values of h)
%     L=zeros(1,10);
%     H=zeros(1,10);
%     for j=1:20
%         h=0.5^j;
%         [U,C] = compute_bssn(1,10,h);
%         Uj=U(:,7);
%         L(1,j)=sqrt(mean((Uj).^2));
%         H(1,j)=h;
%         
%     end
    % plot loglog plot
%     loglog(H,L)
    % plot slope of loglog plot
    %plot(log10(H(1,1:19)),(log10(L(1,2:20)) - log10(L(1,1:19)))./( log10(H(1,2:20)) - log10(H(1,1:19)) ) )
    r_min = 1;
    r_max = 10;
    h = 0.1;
    N=round((r_max-r_min)/h);
    % time at which we want to end the simulation
    t_end=1;
    % number of timesteps to be taken
    dt=0.1;
    %n_it=t_end/dt;
    n_it=1;
    % number of variables to evolve
    n_variables=9;
    n_constraints=3;
    % Eulerian condition (v=0) or Lagrangian condition (v=1)
    v=0;
    % eta parameter in the Gamma-driver condition
    eta = 0;
    
    % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % radial grid
    r=(r_min+h/2:h:r_max-h/2); % might change this later
    % temporal grid
    t=(0:dt:n_it*dt);
    % arrays for the numerical approximations at times n and n+1
    v_new=zeros(N,n_variables);
    v_old=zeros(N,n_variables);
    v_oldold=zeros(N,n_variables);
    % constraint arrays
    c_new = zeros(N,n_constraints);
    c_old = zeros(N,n_constraints);
    
    U = zeros(N,n_variables);
    C = zeros(N,n_constraints);
    
    % flat initial conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = ones(N,1);
    beta_r = zeros(N,1);
    B = zeros(N,1);
    chi = ones(N,1);
    g_rr = ones(N,1); 
    g_thth = r.^2;
    A_rr = zeros(N,1);
    K = zeros(N,1);
    Gamma_r = -2./r;
    % punctured Schwarzchild BH ICs
    g_rr = ones(N,1);
    g_thth = r.*r;
    % tune r0 near 2 but not above.
    chi = cap(1.8,r,1,N);
    
    % initializing the initial state with ICs
    v_old(:,1) = alpha;
    v_old(:,2) = beta_r;
    v_old(:,3) = B;
    v_old(:,4) = chi;
    v_old(:,5) = g_rr;
    v_old(:,6) = g_thth;
    v_old(:,7) = A_rr;
    v_old(:,8) = K;
    v_old(:,9) = Gamma_r;
    
    % time to solve the equations
    tspan = [0 1];
    % setting IC
    y0 = v_old;
    % solving step
    [t,y] = ode45(@(t,y) dydt(t,y,h,N),tspan,y0);
    % getting the output size
    [t_size,y_size] = size(y);
    % unpackaging
    alpha=y(:,1:N);
    beta_r=y(:,N+1:2*N);
    B = y(:,2*N+1:3*N);
    chi=y(:,3*N+1:4*N);
	g_rr=y(:,4*N+1:5*N);
	g_thth=y(:,5*N+1:6*N);
	A_rr=y(:,6*N+1:7*N);
	K=y(:,7*N+1:8*N);
	Gamma_r=y(:,8*N+1:9*N);
    
    % plotting results
    for iter = 1:t_size
        plot(r, beta_r(iter,:))
        pause(0.05)
    end
end


% temporary name, could change it later 
function [U,C]=compute_bssn(r_min, r_max, h)
    % U is the (N,9) matrix with the time derivatives of state variables
    % V is the (N,3) matrix with the time derivatives of constraints
    
    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N=round((r_max-r_min)/h);
    % time at which we want to end the simulation
    t_end=1;
    % number of timesteps to be taken
    dt=0.1;
    %n_it=t_end/dt;
    n_it=1;
    % number of variables to evolve
    n_variables=9;
    n_constraints=3;
    % Eulerian condition (v=0) or Lagrangian condition (v=1)
    v=0;
    % eta parameter in the Gamma-driver condition
    eta = 0;
    % mass parameter for the Schwarzchild BH
    M = 1.;
    
    % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % radial grid
    r=(r_min+h/2:h:r_max-h/2); % might change this later
    % temporal grid
    t=(0:dt:n_it*dt);
    % arrays for the numerical approximations at times n and n+1
    v_new=zeros(N,n_variables);
    v_old=zeros(N,n_variables);
    v_oldold=zeros(N,n_variables);
    % constraint arrays
    c_new = zeros(N,n_constraints);
    c_old = zeros(N,n_constraints);
    
    U = zeros(N,n_variables);
    C = zeros(N,n_constraints);

    
    % the initial conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flat space ICs
    alpha = ones(N,1);
    beta_r = zeros(N,1);
    B = zeros(N,1);
    chi = ones(N,1);
    g_rr = ones(N,1); 
    g_thth = r.^2;
    A_rr = zeros(N,1);
    K = zeros(N,1);
    Gamma_r = -2./r;
    % punctured Schwarzchild BH ICs
    g_rr = ones(N,1);
    g_thth = r.*r;
    % tune r0 near 2 but not above.
    chi = cap(1.8,r,M);
    
    % initializing the initial state with ICs
    v_old(:,1) = alpha;
    v_old(:,2) = beta_r;
    v_old(:,3) = B;
    v_old(:,4) = chi;
    v_old(:,5) = g_rr;
    v_old(:,6) = g_thth;
    v_old(:,7) = A_rr;
    v_old(:,8) = K;
    v_old(:,9) = Gamma_r;
    
    % initial constraint values --- On-shell initial configurations (all constraints vanish at t=0)
    H = c_old(:,1);
    M_r = c_old(:,2);
    G = c_old(:,3);

    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:n_it
        % unpacking our previous state
        alpha=v_old(:,1);
        beta_r=v_old(:,2);
        B = v_old(:,3);
        chi=v_old(:,4);
        g_rr=v_old(:,5);
        g_thth=v_old(:,6);
        A_rr=v_old(:,7);
        K=v_old(:,8);
        Gamma_r=v_old(:,9);
        
        % radial derivatives
        alpha_p = f_prime(alpha,h,1,1,1,1,N);
        alpha_pp = f_pprime(alpha,h,1,1,1,1,N);
        beta_r_p = f_prime(beta_r,h,0,0,0,0,N);
        beta_r_pp = f_pprime(beta_r,h,0,0,0,0,N);
        B_p= f_prime(B,h,0,0,0,0,N);
        chi_p = f_prime(chi,h,1,1,1,1,N);
        chi_pp = f_pprime(chi,h,1,1,1,1,N);
        g_rr_p = f_prime(g_rr,h,1,1,1,1,N);
        g_rr_pp = f_pprime(g_rr,h,1,1,1,1,N);
        g_thth_p = f_prime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        %g_thth_p = f_prime(g_thth,h,(r_min-2*h)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        g_thth_pp = f_pprime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        A_rr_p =  f_prime(A_rr,h,0,0,0,0,N);
        K_p = f_prime(K,h,0,0,0,0,N);
        Gamma_r_p = f_prime(Gamma_r,h,-2/(r_min-3*h/2),-2/(r_min-h/2),-2/(r_max+h/2),-2/(r_max+3*h/2),N);
        H_p = f_prime(H,h,0,0,0,0,N);
        M_r_p = f_prime(M_r,h,0,0,0,0,N);
        M_r_pp = f_pprime(M_r,h,0,0,0,0,N);
        G_p = f_prime(G,h,0,0,0,0,N);
        G_pp = f_pprime(G,h,0,0,0,0,N);
        
        % Building constraints
        H = -3/2*A_rr.*A_rr./g_rr./g_rr + 2/3.*K.*K - 5/2.*chi_p.*chi_p./chi./g_rr...
                +2.*chi_pp./g_rr + 2.*chi./g_thth-2.*chi.*g_thth_pp./g_rr./g_thth...
                +2.*chi_p.*g_thth_p./g_rr./g_thth+chi.*g_rr_p.*g_thth_p./g_rr./g_rr./g_thth...
                -chi_p.*g_rr_p./g_rr./g_rr+chi.*g_thth_p.*g_thth_p./2./g_rr./g_thth./g_thth;
        M_r = A_rr_p./g_rr - 2/3.*K_p...
            - 3/2.*A_rr./g_rr.*(chi_p./chi - g_thth_p./g_thth+ g_rr_p./g_rr);
        G = -g_rr_p./2./g_rr./g_rr+Gamma_r+g_thth_p./g_thth./g_rr;
        
        
        % time derivatives for each state variable
        alpha_t = beta_r.*alpha_p-2*alpha.*K; % eqn 1
        beta_r_t = 3/4*B+beta_r.*beta_r_p; % note this is only 2a), what is B?
        chi_t = 2/3*K.*alpha.*chi - v.*beta_r.*g_rr_p.*chi./(3*g_rr)...
                -2*v.*beta_r.*g_thth_p.*chi./(3*g_thth)-2/3*v.*beta_r_p.*chi...
                +beta_r.*chi_p;
        g_rr_t = -2*A_rr.*alpha-v.*beta_r.*g_rr_p./3+beta_r.*g_rr_p...
                 -2*g_rr.*v.*beta_r.*g_thth_p./(3*g_thth)+2*g_rr.*beta_r_p...
                 -2/3*g_rr.*v.*beta_r_p;
        g_thth_t = A_rr.*g_thth.*alpha./g_rr-g_thth.*v.*beta_r.*g_rr_p./(3*g_rr)...
                   -2/3*v.*beta_r.*g_thth_p+beta_r.*g_thth_p...
                   -2/3*g_thth.*v.*beta_r_p;
        A_rr_t = -2*alpha.*A_rr.^2./g_rr+K.*alpha.*A_rr-v.*beta_r.*g_rr_p.*A_rr./(3*g_rr)...
                 -2*v.*beta_r.*g_thth_p.*A_rr./(3*g_thth)-2/3*v.*beta_r_p.*A_rr...
                 +2*beta_r_p.*A_rr+2*alpha.*chi.*g_rr_p.^2./(3*g_rr.^2)... % end of first line
                 -alpha.*chi.*g_thth_p.^2./(3*g_thth.^2)-alpha.*chi_p.^2./(6*chi)...
                 -2*g_rr.*alpha.*chi./(3*g_thth)+beta_r.*A_rr_p+2/3*g_rr.*alpha.*chi.*Gamma_r_p...
                 -alpha.*chi.*g_rr_p.*g_thth_p./(2*g_rr.*g_thth)+chi.*g_rr_p.*alpha_p./(3*g_rr)... % end of second line
                 +chi.*g_thth_p.*alpha_p./(3*g_thth)-alpha.*g_rr_p.*chi_p./(6*g_rr)...
                 -alpha.*g_thth_p.*chi_p./(6*g_thth)-2/3*alpha_p.*chi_p...
                 -alpha.*chi.*g_rr_pp./(3*g_rr)+alpha.*chi.*g_thth_pp./(3*g_thth)...
                 -2/3*chi.*alpha_pp+alpha.*chi_pp/3;
        K_t = 3*alpha.*A_rr.^2./(2*g_rr.^2)+K.^2.*alpha./3+beta_r.*K_p...
              +chi.*g_rr_p.*alpha_p./(2*g_rr.^2)-chi.*g_thth_p.*alpha_p./(g_rr.*g_thth)...
              +alpha_p.*chi_p./(2*g_rr)-chi.*alpha_pp./g_rr;
        Gamma_r_t = -v.*beta_r.*g_thth_p.^2./(g_rr.*g_thth.^2)...
                    +A_rr.*alpha.*g_thth_p./(g_rr.^2.*g_thth)...
                   -v.*beta_r_p.*g_thth_p./(3*g_rr.*g_thth)+beta_r_p.*g_thth_p./(g_rr.*g_thth)...
                   +beta_r.*Gamma_r_p+A_rr.*alpha.*g_rr_p./(g_rr.^3)...
                   -4/3*alpha.*K_p./g_rr-2*A_rr.*alpha_p./g_rr.^2.... % end of first line
                   +v.*g_rr_p.*beta_r_p./(2*g_rr.^2)-g_rr_p.*beta_r_p./(2*g_rr.^2)...
                   -3*A_rr.*alpha.*chi_p./(g_rr.^2.*chi)+v.*beta_r.*g_rr_pp./(6*g_rr.^2)...
                   +v.*beta_r.*g_thth_pp./(3*g_rr.*g_thth)+v.*beta_r_pp./(3*g_rr)...
                   +beta_r_pp./g_rr;
        % This has to be here since it uses Gamma_r_t
        B_t = Gamma_r_t+beta_r.*(B_p - Gamma_r_p) - eta.*B;
        % This is from 2b). Not sure what the parameter \eta should be but I've set it to 0 for now above.
             
        % Constraint evolution system
        H_t = beta_r.*H + 2/3.*alpha.*K.*H-2.*alpha.*A_rr.*chi./g_rr.*G_p-2.*alpha./g_rr.*chi.*M_r_p...
               + (alpha.*chi_p./g_rr+alpha.*chi.*g_rr_p./g_rr./g_rr-4.*alpha_p.*chi./g_rr...
                    -2.*alpha.*chi.*g_thth_p./g_thth).*M_r;
        M_t = beta_r.*M_r_p + beta_r_p.*M_r - alpha.*K.*M_r - alpha_p./3.*H + alpha./3.*H_p + 2/3.*alpha.*chi.*G_pp...
                +(2/3.*alpha_p.*chi - alpha.*chi_p./3 + alpha.*chi.*g_thth_p./g_thth).*G_p;
        G_t = beta_r.*G_p+ 2.*alpha./g_rr.*M_r;
        
        % this part is for time evolution, not implemented yet
        v_old=v_new;
        c_old=c_new;
    end
    % function outputs, A(N,9) is the time derivatives
    % B(N,3) are the constraint time derivatives
    U(:,1) = alpha_t;
    U(:,2) = beta_r_t;
    U(:,3) = B_t;
    U(:,4) = chi_t;
    U(:,5) = g_rr_t;
    U(:,6) = g_thth_t;
    U(:,7) = A_rr_t;
    U(:,8) = K_t;
    U(:,9) = Gamma_r_t;
    C(:,1) = H_t;
    C(:,2) = M_t;
    C(:,3) = G_t;
    
   
end

function U=dydt(t, v_old, h, N)
    eta = 1;
    v = 0;
    U = zeros(9*N,1);
    % unpacking our previous state
    [alpha, beta_r, B, chi, g_rr, g_thth, A_rr, K, Gamma_r] = unpackage(v_old, N);

    % radial derivatives
    alpha_p = f_prime(alpha,h,alpha(2),alpha(1),alpha(N),alpha(N-1),N);
    alpha_pp = f_pprime(alpha,h,alpha(2),alpha(1),alpha(N),alpha(N-1),N);
    beta_r_p = f_prime(beta_r,h,-beta_r(2),-beta_r(1),-beta_r(N),-beta_r(N-1),N);
    beta_r_pp = f_pprime(beta_r,h,-beta_r(2),-beta_r(1),-beta_r(N),-beta_r(N-1),N);
    B_p= f_prime(B,h,-B(2),-B(1),-B(N),-B(N-1),N);
    chi_p = f_prime(chi,h,chi(2),chi(1),chi(N),chi(N-1),N);
    chi_pp = f_pprime(chi,h,chi(2),chi(1),chi(N),chi(N-1),N);
    g_rr_p = f_prime(g_rr,h,g_rr(2),g_rr(1),g_rr(N),g_rr(N-1),N);
    g_rr_pp = f_pprime(g_rr,h,g_rr(2),g_rr(1),g_rr(N),g_rr(N-1),N);
    g_thth_p = f_prime(g_thth,h,g_thth(2), g_thth(1), g_thth(N), g_thth(N-1),N);
    g_thth_pp = f_pprime(g_thth,h,g_thth(2), g_thth(1), g_thth(N), g_thth(N-1),N);
    A_rr_p =  f_prime(A_rr,h,A_rr(2),A_rr(1),A_rr(N),A_rr(N-1),N);
    K_p = f_prime(K,h,K(2),K(1),K(N),K(N-1),N);
    Gamma_r_p = f_prime(Gamma_r,h,-Gamma_r(2),-Gamma_r(1),-Gamma_r(N), -Gamma_r(N-1) ,N);
    
    % Building constraints
    H = -3/2*A_rr.*A_rr./g_rr./g_rr + 2/3.*K.*K - 5/2.*chi_p.*chi_p./chi./g_rr...
            +2.*chi_pp./g_rr + 2.*chi./g_thth-2.*chi.*g_thth_pp./g_rr./g_thth...
            +2.*chi_p.*g_thth_p./g_rr./g_thth+chi.*g_rr_p.*g_thth_p./g_rr./g_rr./g_thth...
            -chi_p.*g_rr_p./g_rr./g_rr+chi.*g_thth_p.*g_thth_p./2./g_rr./g_thth./g_thth;
    M_r = A_rr_p./g_rr - 2/3.*K_p...
        - 3/2.*A_rr./g_rr.*(chi_p./chi - g_thth_p./g_thth+ g_rr_p./g_rr);
    G = -g_rr_p./2./g_rr./g_rr+Gamma_r+g_thth_p./g_thth./g_rr;
    
    H_p = f_prime(H,h,0,0,0,0,N);
    M_r_p = f_prime(M_r,h,0,0,0,0,N);
    M_r_pp = f_pprime(M_r,h,0,0,0,0,N);
    G_p = f_prime(G,h,0,0,0,0,N);
    G_pp = f_pprime(G,h,0,0,0,0,N);
    
    % time derivatives for each state variable
    alpha_t = beta_r.*alpha_p-2*alpha.*K; % eqn 1
    beta_r_t = 3/4*B+beta_r.*beta_r_p; % note this is only 2a), what is B?
    chi_t = 2/3*K.*alpha.*chi - v.*beta_r.*g_rr_p.*chi./(3*g_rr)...
            -2*v.*beta_r.*g_thth_p.*chi./(3*g_thth)-2/3*v.*beta_r_p.*chi...
            +beta_r.*chi_p;
    g_rr_t = -2*A_rr.*alpha-v.*beta_r.*g_rr_p./3+beta_r.*g_rr_p...
             -2*g_rr.*v.*beta_r.*g_thth_p./(3*g_thth)+2*g_rr.*beta_r_p...
             -2/3*g_rr.*v.*beta_r_p;
    g_thth_t = A_rr.*g_thth.*alpha./g_rr-g_thth.*v.*beta_r.*g_rr_p./(3*g_rr)...
               -2/3*v.*beta_r.*g_thth_p+beta_r.*g_thth_p...
               -2/3*g_thth.*v.*beta_r_p;
    A_rr_t = -2*alpha.*A_rr.^2./g_rr+K.*alpha.*A_rr-v.*beta_r.*g_rr_p.*A_rr./(3*g_rr)...
             -2*v.*beta_r.*g_thth_p.*A_rr./(3*g_thth)-2/3*v.*beta_r_p.*A_rr...
             +2*beta_r_p.*A_rr+2*alpha.*chi.*g_rr_p.^2./(3*g_rr.^2)... % end of first line
             -alpha.*chi.*g_thth_p.^2./(3*g_thth.^2)-alpha.*chi_p.^2./(6*chi)...
             -2*g_rr.*alpha.*chi./(3*g_thth)+beta_r.*A_rr_p+2/3*g_rr.*alpha.*chi.*Gamma_r_p...
             -alpha.*chi.*g_rr_p.*g_thth_p./(2*g_rr.*g_thth)+chi.*g_rr_p.*alpha_p./(3*g_rr)... % end of second line
             +chi.*g_thth_p.*alpha_p./(3*g_thth)-alpha.*g_rr_p.*chi_p./(6*g_rr)...
             -alpha.*g_thth_p.*chi_p./(6*g_thth)-2/3*alpha_p.*chi_p...
             -alpha.*chi.*g_rr_pp./(3*g_rr)+alpha.*chi.*g_thth_pp./(3*g_thth)...
             -2/3*chi.*alpha_pp+alpha.*chi_pp/3;
    K_t = 3*alpha.*A_rr.^2./(2*g_rr.^2)+K.^2.*alpha./3+beta_r.*K_p...
          +chi.*g_rr_p.*alpha_p./(2*g_rr.^2)-chi.*g_thth_p.*alpha_p./(g_rr.*g_thth)...
          +alpha_p.*chi_p./(2*g_rr)-chi.*alpha_pp./g_rr;
    Gamma_r_t = -v.*beta_r.*g_thth_p.^2./(g_rr.*g_thth.^2)...
                +A_rr.*alpha.*g_thth_p./(g_rr.^2.*g_thth)...
               -v.*beta_r_p.*g_thth_p./(3*g_rr.*g_thth)+beta_r_p.*g_thth_p./(g_rr.*g_thth)...
               +beta_r.*Gamma_r_p+A_rr.*alpha.*g_rr_p./(g_rr.^3)...
               -4/3*alpha.*K_p./g_rr-2*A_rr.*alpha_p./g_rr.^2.... % end of first line
               +v.*g_rr_p.*beta_r_p./(2*g_rr.^2)-g_rr_p.*beta_r_p./(2*g_rr.^2)...
               -3*A_rr.*alpha.*chi_p./(g_rr.^2.*chi)+v.*beta_r.*g_rr_pp./(6*g_rr.^2)...
               +v.*beta_r.*g_thth_pp./(3*g_rr.*g_thth)+v.*beta_r_pp./(3*g_rr)...
               +beta_r_pp./g_rr;
    % This has to be here since it uses Gamma_r_t
    B_t = Gamma_r_t+beta_r.*(B_p - Gamma_r_p) - eta.*B;
    % This is from 2b). Not sure what the parameter \eta should be but I've set it to 0 for now above.
    U = [alpha_t; beta_r_t; B_t; chi_t; g_rr_t; g_thth_t; A_rr_t; K_t; Gamma_r_t];
    
end

% this function unpackages the *compressed* state vector v
% N is the length of the grid
function [alpha, beta_r, B, chi, g_rr, g_thth, A_rr, K, Gamma_r]=unpackage(v,N)
    alpha = v(1:N);
    beta_r = v(N+1:2*N);
    B = v(2*N+1:3*N);
    chi=v(3*N+1:4*N);
    g_rr=v(4*N+1:5*N);
    g_thth=v(5*N+1:6*N);
    A_rr=v(6*N+1:7*N);
    K=v(7*N+1:8*N);
    Gamma_r=v(8*N+1:9*N);
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
    chi0 = power(f,4);
end

% Kreiss-Oliger dissipation term
% since our finite differencing is fourth order we have 6th order third
% derivative, which requires 4 points to each side
function v=KreissOliger(f,sigma,h,a,b,c,d,w,x,y,z,N)
    v=zeros(N,1);
    prefactor = (sigma*(-1)^6*h^5/2^6);
    v(1) = prefactor*(7/240*f(5)-3/10*f(4)+169/120*f(3)-61/30*f(2)...
           +61/30*d-169/120*c+3/10*b-7/240*a)./h^3;
    v(2) = prefactor*(7/240*f(6)-3/10*f(5)+169/120*f(4)-61/30*f(3)...
           +61/30*f(1)-169/120*d+3/10*c-7/240*b)./h^3;       
    v(3) = prefactor*(7/240*f(7)-3/10*f(6)+169/120*f(5)-61/30*f(4)...
           +61/30*f(2)-169/120*f(1)+3/10*d-7/240*c)./h^3;
    v(4) = prefactor*(7/240*f(8)-3/10*f(7)+169/120*f(6)-61/30*f(5)...
           +61/30*f(3)-169/120*f(2)+3/10*f(1)-7/240*d)./h^3;    
    v(N-3) = prefactor*(7/240*w-3/10*f(N)+169/120*f(N-1)-61/30*f(N-2)...
             +61/30*f(N-4)-169/120*f(N-5)+3/10*f(N-6)-7/240*f(N-7))./h^3;  
    v(N-2) = prefactor*(7/240*x-3/10*w+169/120*f(N)-61/30*f(N-1)...
             +61/30*f(N-3)-169/120*f(N-4)+3/10*f(N-5)-7/240*f(N-6))./h^3;         
    v(N-1) = prefactor*(7/240*y-3/10*x+169/120*w-61/30*f(N)...
             +61/30*f(N-2)-169/120*f(N-3)+3/10*f(N-4)-7/240*f(N-5))./h^3;    
    v(N) = prefactor*(7/240*z-3/10*y+169/120*x-61/30*w...
           +61/30*f(N-1)-169/120*f(N-2)+3/10*f(N-3)-7/240*f(N-4))./h^3;          
    % Computing the middle parts
    v(5:N-4) = prefactor*(7/240*f(9:N)-3/10*f(8:N-1)+169/120*f(7:N-2)-61/30*f(6:N-3)...
               +61/30*f(4:N-5)-169/120*f(3:N-6)+3/10*f(2:N-7)-7/240*f(1:N-8))./h^3;
end
