function bssn

    % Winter 2021
    % Assignment C1

    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of grid points in space (not counting the point at x=1, which is the same as the point at x=0, due to the periodicity)
    r_max=10;
    r_min=1;
    % spatial step
    h =0.1;
    N=(r_max-r_min)/h;
    % time at which we want to end the simulation
    t_end=1;
    % number of timesteps to be taken
    dt=0.01;
    n_it=t_end/dt;
    % number of variables to evolve
    n_variables=2;
    % Eulerian condition (v=0) or Lagrangian condition (v=1)
    v=0;
    % Theta?
    theta=pi/2;

    % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % radial grid
    r=(r_min:h:r_max); % might change this later
    % temporal grid
    t=(0:dt:n_it*dt);
    % arrays for the numerical approximations at times n and n+1
    v_new=zeros(N,n_variables);
    v_old=zeros(N,n_variables);
    v_oldold=zeros(N,n_variables);

    
    % the initial condition (need to implement)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_old(:,1) = 1; 
    v_old(:,2) = 1;
    v_old(:,3) = 1;
    v_old(:,4) = 1;
    v_old(:,5) = 1; 
    v_old(:,6) = r.^2;
    v_old(:,7) = r.^2 * sin(theta);
    v_old(:,8) = 1;
    v_old(:,9) = 1;
    


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
        alpha_p = f_prime(alpha,h,a,b,c,d);
        alpha_pp = f_pprime(alpha,h,a,b,c,d);
        beta_p = f_prime(beta,h,a,b,c,d);
        beta_pp = f_pprime(beta,h,a,b,c,d);
        chi_p = f_prime(chi,h,a,b,c,d);
        chi_pp = f_pprime(chi,h,a,b,c,d);
        g_rr_p = f_prime(g_rr,h,a,b,c,d);
        g_rr_pp = f_pprime(g_rr,h,a,b,c,d);
        g_thth_p = f_prime(g_thth,h,a,b,c,d);
        g_thth_pp = f_pprime(g_thth,h,a,b,c,d);
        K_p = f_prime(K,h,a,b,c,d);
        Gamma_r_p = f_prime(Gamma_r,h,a,b,c,d);
       
        
        % time derivatives for each state variable (this part needs dbl
        % checking)
        alpha_t = beta.*alpha_p-2*alpha.*K; % eqn 1
        beta_t = 3/4*B+beta.*beta_p; % note this is only 2a), what is B?
        chi_t = 2/3*K.*alpha.*chi - v.*beta.*g_rr_p.*chi./(3*g_rr)...
                -2*v.*beta.*g_thth_p.*chi./(3*g_thth)-2/3*v.*beta_p.*chi...
                +beta.*chi_p;
        g_rr_t = -2*A_rr.*alpha-v.*beta.*g_rr_p./3+beta.*g_rr_p...
                 -2*g_rr.*v.*beta.*g_thth_p./(3*g_thth)+2*g_rr.*beta_p...
                 -2/3*g_rr.*v.*beta_p;
        g_thth_t = A_rr.*g_thth.*alpha./g_rr-g_thth.*v.*beta.*g_rr_p./(3*g_rr)...
                   -2/3*v.*beta.*g_thth_p+beta.*g_thth_p...
                   -2/3*g_thth.*v.*beta_p;
        A_rr_t = -2*alpha.*A_rr.^2./g_rr+K.*alpha.*A_rr-v.*beta.*g_rr_p.*A_rr./(3*g_rr)...
                 -2*v.*beta.*g_thth_p.*A_rr./(3*g_thth)-3/2*v.*beta_p.*A_rr...
                 +2*alpha.*chi.*g_rr_p.^2./(3*g_rr.^2)... % end of first line
                 -alpha.*chi.*g_thth_p.^2./(3*g_thth.^2)-alpha.*chi_p.^2/(6*chi)...
                 -2*g_rr.*alpha.*chi./(3*g_thth)+beta.*A_rr_p+2/3*g_rr.*alpha.*chi.*Gamma_r_p...
                 -alpha.*chi.*g_rr_p.*g_thth_p./(2*g_rr.*g_thth)+chi.*g_rr_p.*alpha_p./(3*g_rr)... % end of second line
                 +chi.*g_thth_p.*alpha_p./(3*g_thth)-alpha.*g_rr_p.*chi_p./(6*g_rr)...
                 -alpha.*g_thth_p.*chi_p./(6*g_thth)-2/3*alpha_p.*chi_p...
                 -alpha.*chi.*g_rr_pp./(3*g_rr)+alpha.*chi.*g_thth_pp./(3*g_thth)...
                 -2/3*chi.*alpha_pp+alpha.*chi_pp/3;
        K_t = 3*alpha.*A_rr.^2/(2*g_rr.^2)+K.^2*alpha./3+beta.*K_p...
              +chi.*g_rr_p.*alpha_p./(2*g_rr.^2)-chi.*g_thth_p.*alpha_p./g_rr./g_thth...
              +alpha_p.*chi_p./(2*g_rr)-chi.*alpha_pp./g_rr;
        Gamma_r_t = -v.*beta.*g_thth_p.^2/(g_rr.*g_thth.^2)...
                    +A_rr.*alpha.*g_thth_p./(g_rr.^2.*g_thth)...
                   -v.*beta_p.*g_thth_p./(3*g_rr.*g_thth)+beta_p.*g_thth_p./(g_rr.*g_thth)...
                   +beta.*Gamma_r_p+A_rr.*alpha.*g_rr_p./(g_rr.^3)...
                   -4/3*alpha.*K_p./g_rr-2*A_rr.*alpha_p./g_rr.^2.... % end of first line
                   +v.*g_rr_p.*beta_p./(2*g_rr.^2)-g_rr_p.*beta_p./(2*g_rr.^2)...
                   -3*A_rr.*alpha.*chi_p./(g_rr.^2.*chi)+v.*beta.*g_rr_pp./(6*g_rr.^2)...
                   +v.*beta.*g_thth_pp./(3*g_rr.*g_thth)+v.*beta_pp./(3*g_rr)...
                   +beta_pp./g_rr;
             
        v_oldold=v_old;
        v_old=v_new;
    end
end

% This function returns f'(x) where f is one of the state variables
function y=f_prime(f,h,a,b,c,d)
    y=zeros(N);
    % These are for the two-level boundary conditions at each end
    y(1) = a;
    y(2) = b;
    y(N-1) = c;
    y(N) = d;
    % Computing the middle parts
    y(3:N-2) = (-f(5:N) + 8*f(4:N-1) - 8*f(2:N-3) + f(1:N-4))./(12*h);
end

% This function returns f''(x) where f is one of the state variables
function y=f_pprime(f,h,a,b,c,d)
    y=zeros(N);
    % These are for the two-level boundary conditions at each end
    y(1) = a;
    y(2) = b;
    y(N-1) = c;
    y(N) = d;
    % Computing the middle parts
    y(3:N-2) = (-f(5:N) + 16*f(4:N-1) - 30*f(3:N-2) + 16*f(2:N-3) - f(1:N-4))./(12*h^2);
end