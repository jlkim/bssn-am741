function bssn
    % call the other function here
    % example:
    %[A,B] = compute_bssn(1, 2, 0.000001)
    % This will return a matrix A of time derivatives of state variables
    % and a matrix B of time derivatives of constraints
    %defining L(containing RMS of A_rr_t) and H(the values of h)
    L=zeros(1,10)
    H=zeros(1,10)
    for j=1:5
        h=0.1^j
        [A,B] = compute_bssn(1,2,h);
        Aj=A(:,7)
        L(1,j)=sqrt(mean((Aj).^2));
        H(1,j)=h
        
    end
    L
    H
    %plotting
    loglog(H,L)



end

 

% temporary name, could change it later 
function [A,B]=compute_bssn(r_min, r_max, h)
    % A is the (N,9) matrix with the time derivatives of state variables
    % B is the (N,3) matrix with the time derivatives of constraints
    
    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N=ceil((r_max-r_min)/h);
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
    
    A = zeros(N,n_variables);
    B = zeros(N,n_constraints);

    
    % the initial conditions
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
              +chi.*g_rr_p.*alpha_p./(2*g_rr.^2)-chi.*g_thth_p.*alpha_p./g_rr./g_thth...
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
    A(:,1) = alpha_t;
    A(:,2) = beta_r_t;
    A(:,3) = B_t;
    A(:,4) = chi_t;
    A(:,5) = g_rr_t;
    A(:,6) = g_thth_t;
    A(:,7) = A_rr_t;
    A(:,8) = K_t;
    A(:,9) = Gamma_r_t;
    B(:,1) = H_t;
    B(:,2) = M_t;
    B(:,3) = G_t;
    
    
    
   
end

function y=time_derivative(v, h, N)
    % TODO: move time derivatives here
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
