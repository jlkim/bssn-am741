function bssn

    % Winter 2021
    % Assignment C1

    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of grid points in space (not counting the point at x=1, which is the same as the point at x=0, due to the periodicity)
    r_max=2;
    r_min=1;
    % spatial step
    h =0.1;
    N=(r_max-r_min)/h;
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
    eta = 0.;

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

    
    % the initial condition (need to implement)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_old(:,1) = ones(N,1);
    v_old(:,2) = zeros(N,1);
    v_old(:,3) = zeros(N,1);
    v_old(:,4) = ones(N,1);
    v_old(:,5) = ones(N,1); 
    v_old(:,6) = r.^2;
    v_old(:,7) = zeros(N,1);
    v_old(:,8) = zeros(N,1);
    v_old(:,9) = -2./r;
    % unpacking the values so it's more readable
    alpha=v_old(:,1);
    beta=v_old(:,2);
    B = v_old(:,3);
    chi=v_old(:,4);
    g_rr=v_old(:,5);
    g_thth=v_old(:,6);
    A_rr=v_old(:,7);
    K=v_old(:,8);
    Gamma_r=v_old(:,9);
    
    % initial constraint values --- On-shell initial configurations (all constraints vanish at t=0)
    H = c_old(:,1);
    M = c_old(:,2);
    G = c_old(:,3);

    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:n_it
        % unpacking our previous state
        alpha=v_old(:,1);
        beta=v_old(:,2);
        B = v_old(:,3);
        chi=v_old(:,4);
        g_rr=v_old(:,5);
        g_thth=v_old(:,6);
        A_rr=v_old(:,7);
        K=v_old(:,8);
        Gamma_r=v_old(:,9);
        
        % radial derivatives (the boundary conditions are incomplete!)
        alpha_p = f_prime(alpha,h,1,1,1,1,N);
        alpha_pp = f_pprime(alpha,h,1,1,1,1,N);
        % Could also compute with alpha_pp = f_prime(alpha_p,h,1,1,1,1,N)
        beta_p = f_prime(beta,h,0,0,0,0,N);
        beta_pp = f_pprime(beta,h,0,0,0,0,N);
        B_p= f_prime(B,h,0,0,0,0,N);
        % I'm also unsure about what B exactly is (a gauge fixing quantity?), so I'm not sure what would be the right boundary condition. I just set them all to 0 for now. 
        chi_p = f_prime(chi,h,1,1,1,1,N);
        chi_pp = f_pprime(chi,h,1,1,1,1,N);
        g_rr_p = f_prime(g_rr,h,1,1,1,1,N);
        g_rr_pp = f_pprime(g_rr,h,1,1,1,1,N);
        g_thth_p = f_prime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        g_thth_pp = f_pprime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
        A_rr_p =  f_pprime(A_rr,h,0,0,0,0,N);
        K_p = f_prime(K,h,0,0,0,0,N);
        Gamma_r_p = f_prime(Gamma_r,h,1/(r_min-3*h/2),1/(r_min-h/2),1/(r_max-h/2),1/(r_max+3*h/2),N);
        % BCs. for constraints depends on the Cauchy surface; pick finite energy density at inner BH boundary?
        H_p = f_prime(H,h,1,1,0,0,N);
        M_p = f_prime(M,h,0,0,0,0,N);
        M_pp = f_prime(M_p,h,0,0,0,0,N);
        G_p = f_prime(G,h,0,0,0,0,N);
        G_pp = f_prime(G_p,h,0,0,0,0,N);
        
        % Building constraints
        H = -3/2*A_rr.*A_rr./g_rr./g_rr + 2/3.*K.*K - 5/2.*chi_p.*chi_p./chi./g_rr...
                +2.*chi_pp./g_rr + 2.*chi./g_thth-2.*chi.*g_thth_pp./g_rr./g_thth...
                +2.*chi_p.*g_thth_p./g_rr./g_thth+chi.*g_rr_p.*g_thth_p./g_rr./g_rr./g_thth...
                -chi_p.*g_rr_p./g_rr./g_rr+chi.*g_thth_p.*g_thth_p./2./g_rr./g_thth./g_thth;
        M = A_rr_p./g_rr - 2/3.*K_p...
            - 3/2.*A_rr./g_rr.*(chi_p./chi - g_thth_p./g_thth+ g_rr_p./g_rr);
        G = -g_rr_p./2./g_rr./g_rr+Gamma_r+g_thth_p./g_thth./g_rr;
        
        
        % time derivatives for each state variable (this part needs dbl
        % checking)
        alpha_t = beta.*alpha_p-2*alpha.*K % eqn 1
        beta_t = 3/4*B+beta.*beta_p % note this is only 2a), what is B?
        chi_t = 2/3*K.*alpha.*chi - v.*beta.*g_rr_p.*chi./(3*g_rr)...
                -2*v.*beta.*g_thth_p.*chi./(3*g_thth)-2/3*v.*beta_p.*chi...
                +beta.*chi_p
        g_rr_t = -2*A_rr.*alpha-v.*beta.*g_rr_p./3+beta.*g_rr_p...
                 -2*g_rr.*v.*beta.*g_thth_p./(3*g_thth)+2*g_rr.*beta_p...
                 -2/3*g_rr.*v.*beta_p
        g_thth_t = A_rr.*g_thth.*alpha./g_rr-g_thth.*v.*beta.*g_rr_p./(3*g_rr)...
                   -2/3*v.*beta.*g_thth_p+beta.*g_thth_p...
                   -2/3*g_thth.*v.*beta_p
        A_rr_t = -2*alpha.*A_rr.^2./g_rr+K.*alpha.*A_rr-v.*beta.*g_rr_p.*A_rr./(3*g_rr)...
                 -2*v.*beta.*g_thth_p.*A_rr./(3*g_thth)-3/2*v.*beta_p.*A_rr...
                 +2*alpha.*chi.*g_rr_p.^2./(3*g_rr.^2)... % end of first line
                 -alpha.*chi.*g_thth_p.^2./(3*g_thth.^2)-alpha.*chi_p.^2./(6*chi)...
                 -2*g_rr.*alpha.*chi./(3*g_thth)+beta.*A_rr_p+2/3*g_rr.*alpha.*chi.*Gamma_r_p...
                 -alpha.*chi.*g_rr_p.*g_thth_p./(2*g_rr.*g_thth)+chi.*g_rr_p.*alpha_p./(3*g_rr)... % end of second line
                 +chi.*g_thth_p.*alpha_p./(3*g_thth)-alpha.*g_rr_p.*chi_p./(6*g_rr)...
                 -alpha.*g_thth_p.*chi_p./(6*g_thth)-2/3*alpha_p.*chi_p...
                 -alpha.*chi.*g_rr_pp./(3*g_rr)+alpha.*chi.*g_thth_pp./(3*g_thth)...
                 -2/3*chi.*alpha_pp+alpha.*chi_pp/3
        K_t = 3*alpha.*A_rr.^2./(2*g_rr.^2)+K.^2.*alpha./3+beta.*K_p...
              +chi.*g_rr_p.*alpha_p./(2*g_rr.^2)-chi.*g_thth_p.*alpha_p./g_rr./g_thth...
              +alpha_p.*chi_p./(2*g_rr)-chi.*alpha_pp./g_rr
        Gamma_r_t = -v.*beta.*g_thth_p.^2./(g_rr.*g_thth.^2)...
                    +A_rr.*alpha.*g_thth_p./(g_rr.^2.*g_thth)...
                   -v.*beta_p.*g_thth_p./(3*g_rr.*g_thth)+beta_p.*g_thth_p./(g_rr.*g_thth)...
                   +beta.*Gamma_r_p+A_rr.*alpha.*g_rr_p./(g_rr.^3)...
                   -4/3*alpha.*K_p./g_rr-2*A_rr.*alpha_p./g_rr.^2.... % end of first line
                   +v.*g_rr_p.*beta_p./(2*g_rr.^2)-g_rr_p.*beta_p./(2*g_rr.^2)...
                   -3*A_rr.*alpha.*chi_p./(g_rr.^2.*chi)+v.*beta.*g_rr_pp./(6*g_rr.^2)...
                   +v.*beta.*g_thth_pp./(3*g_rr.*g_thth)+v.*beta_pp./(3*g_rr)...
                   +beta_pp./g_rr
        % This has to be here since it uses Gamma_r_t
        B_t = Gamma_r_t+beta.*(B_p - Gamma_r_p) - eta.*B;
        % This is from 2b). Not sure what the parameter \eta should be but I've set it to 0 for now above.
             
        % Constraint evolution system
        H_t = beta.*H + 2/3.*alpha.*K.*H-2.*alpha.*A_rr.*chi./g_rr.*G_p-2.*alpha./g_rr.*chi.*M_p...
               + (alpha.*chi_p./g_rr+alpha.*chi.*g_rr_p./g_rr./g_rr-4.*alpha_p.*chi./g_rr...
                    -2.*alpha.*chi.*g_thth_p./g_thth).*M
        M_t = beta.*M_p + beta_p.*M - alpha.*K.*M - alpha_p./3.*H + alpha./3.*H_p + 2/3.*alpha.*chi.*G_pp...
                +(2/3.*alpha_p.*chi - alpha.*chi_p./3 + alpha.*chi.*g_thth_p./g_thth).*G_p
        G_t = beta.*G_p+ 2.*alpha./g_rr.*M
        
        v_old=v_new;
        c_old = c_new;
    end
end

% Do we need both fprime and fpprime? Couldn't we compute f'' by putting f' into f_prime?

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
