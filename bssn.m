function bssn
global v
global eta
global n_var
global r_min
global r_max
v = 0; % Eulerian condition (v=0) or Lagrangian condition (v=1)
eta = 0; % eta parameter in the Gamma-driver condition
n_var = 12; % number of functions to evolve
r_min = 0;
r_max = 10;
h = 0.1; % spatial grid
N=round((r_max-r_min)/h);
r=(r_min+h/2:h:r_max-h/2);
% time at which we want to end the simulation
t_end=1;

% time to solve the equations
tspan = [0 t_end];
% initial condition
y0 = initial_cond(h,r_min,r_max,0)
% if we want to plot the RHS we call the function below where the
% inputs are plot_RHS(num_iterations, r_min, r_max, puncture (0 or 1))
%plot_RHS(20,1,10,0)
% solving step
[t,y] = ode45(@(t,y) dydt(t,y,h,N),tspan,y0);
%getting the output size
[t_size,y_size] = size(y)
% reshaping the array so it's nicer to unpackage
for iter=1:t_size
    U(iter,:,:)=reshape(y(iter,:),[],n_var);
end
%unpackaging
alpha=U(:,:,1);
beta_r=U(:,:,2);
B_r = U(:,:,3);
chi=U(:,:,4);
g_rr=U(:,:,5);
g_thth=U(:,:,6);
A_rr=U(:,:,7);
K=U(:,:,8);
Gamma_r=U(:,:,9);
H=U(:,:,10);
M_r=U(:,:,11);
G=U(:,:,12);
%plotting results
%     for iter = 1:t_size
%        plot(r, alpha(iter,:))
%        pause(0.005)
%     end
for iter=1:t_size
    Q(iter,:,:)=reshape(y(iter,:),[],n_var);
end

%unpackaging
charac1=Q(:,:,1);
charac2=Q(:,:,2);
%plotting
for iter=1:t_size
    plot(r,charac1(iter,:))
    pause(0.005)
end
end

function U=dydt(t, v_old, h, N)
global n_var
global r_min
global r_max
eta = 0;
v = 0;
U = zeros(n_var*N,1);
% unpacking our previous state
v_old=reshape(v_old,[],n_var);
r=(r_min+h/2:h:r_max-h/2);
[alpha, beta_r, B_r, chi, g_rr, g_thth, A_rr, K, Gamma_r, H, M_r, G] = unpackage(v_old);



% radial derivatives
alpha_p = f_prime(alpha,h,alpha(2),alpha(1),1,1,N);
alpha_pp = f_pprime(alpha,h,alpha(2),alpha(1),1,1,N);
beta_r_p = f_prime(beta_r,h,-beta_r(2),-beta_r(1),0,0,N);
beta_r_pp = f_pprime(beta_r,h,-beta_r(2),-beta_r(1),0,0,N);
B_r_p= f_prime(B_r,h,-B_r(2),-B_r(1),0,0,N);
chi_p = f_prime(chi,h,chi(2),chi(1),1,1,N);
chi_pp = f_pprime(chi,h,chi(2),chi(1),1,1,N);
g_rr_p = f_prime(g_rr,h,g_rr(2),g_rr(1),1,1,N);
g_rr_pp = f_pprime(g_rr,h,g_rr(2),g_rr(1),1,1,N);
g_thth_p = f_prime(g_thth,h,g_thth(2), g_thth(1), (r_max+h/2)^2,(r_max+3*h/2)^2,N);
g_thth_pp = f_pprime(g_thth,h,g_thth(2), g_thth(1),(r_max+h/2)^2,(r_max+3*h/2)^2,N);
A_rr_p =  f_prime(A_rr,h,A_rr(2),A_rr(1),0,0,N);
K_p = f_prime(K,h,K(2),K(1),0,0,N);
Gamma_r_p = f_prime(Gamma_r,h,-2/(r_min-3*h/2),-2/(r_min-h/2),-2/(r_max+h/2),-2/(r_max+3*h/2),N);
% The initial conditions below are what we used before for r_min = 1.
% at r_min = 0 some stuff go wrong
%g_thth_p = f_prime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
%g_thth_pp = f_pprime(g_thth,h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
%Gamma_r_p = f_prime(Gamma_r,h,-Gamma_r(2),-Gamma_r(1),-2/(r_max+h/2),-2/(r_max+3*h/2) ,N);

% constraint evolution equations
H_p = f_prime(H,h,H(2),H(1),0,0,N);
M_r_p = f_prime(M_r,h,M_r(2),M_r(1),0,0,N);
M_r_pp = f_pprime(M_r,h,M_r(2),M_r(1),0,0,N);
G_p = f_prime(G,h,G(2),G(1),0,0,N);
G_pp = f_pprime(G,h,G(2),G(1),0,0,N);

% time derivatives for each state variable
alpha_t = beta_r.*alpha_p-2*alpha.*K; % eqn 1
beta_r_t = 3/4*B_r+beta_r.*beta_r_p; % note this is only 2a), what is B_r?
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
B_r_t = Gamma_r_t+beta_r.*(B_r_p - Gamma_r_p) - eta.*B_r;
% This is from 2b). Not sure what the parameter \eta should be but I've set it to 0 for now above.

% Constraint evolution system
H_t = beta_r.*H + 2/3.*alpha.*K.*H-2.*alpha.*A_rr.*chi./g_rr.*G_p-2.*alpha./g_rr.*chi.*M_r_p...
    + (alpha.*chi_p./g_rr+alpha.*chi.*g_rr_p./g_rr./g_rr-4.*alpha_p.*chi./g_rr...
    -2.*alpha.*chi.*g_thth_p./g_thth).*M_r;
M_r_t = beta_r.*M_r_p + beta_r_p.*M_r - alpha.*K.*M_r - alpha_p./3.*H + alpha./3.*H_p + 2/3.*alpha.*chi.*G_pp...
    +(2/3.*alpha_p.*chi - alpha.*chi_p./3 + alpha.*chi.*g_thth_p./g_thth).*G_p;
G_t = beta_r.*G_p+ 2.*alpha./g_rr.*M_r;

U = [alpha_t; beta_r_t; B_r_t; chi_t; g_rr_t; g_thth_t; A_rr_t; K_t; Gamma_r_t; H_t; M_r_t; G_t];

end


function Q=CHARAC(t, Z_old, h, N)
n_varr=2
global r_min
global r_max
Q = zeros(n_varr*N,1);
% unpacking our previous state
Z_old=reshape(Z_old,[],n_varr);
r=(r_min+h/2:h:r_max-h/2);
[Gamma_r,A_rr,g_rr,chi,g_rr,g_thth] = unpackage1(Z_old);
for iter=1:t_size
    
    chi_p = f_prime(chi(iter,:),h,1,1,1,1,N);
    g_rr_p = f_prime(g_rr(iter,:),h,1,1,1,1,N);
    g_thth_p = f_prime(g_thth(iter,:),h,(r_min-3*h/2)^2,(r_min-h/2)^2,(r_max+h/2)^2,(r_max+3*h/2)^2,N);
    
    
    
    charac1= Gamma_r - 3/2.*A_rr./((g_rr.^3).*chi).^0.5 +1/2.*chi_p./(2.*g_rr.*chi)- 1/2.*g_rr_p./(2.*g_rr.^2)...
        +g_thth_p./(2.*g_rr.*g_thth)+ K./(g_rr.*chi).^0.5;
    
    charac2=Gamma_r + 3/2.*A_rr./((g_rr.^3).*chi).^0.5 +1/2.*chi_p./(2.*g_rr.*chi)- 1/2.*g_rr_p./(2.*g_rr.^2)...
        +g_thth_p./(2.*g_rr.*g_thth)- K./(g_rr.*chi).^0.5;
    
    Q=[charac1,charac2] ;
end


end
function [charac1,charac2]=unpackage1(Z)
charac1 = Z(:,1);
charac2 = Z(:,2);

end






function v_old=initial_cond(h,r_min,r_max,puncture)
% puncture == 1 means we have a puncture, otherwise we have flat space
global n_var
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
A_rr = zeros(N,1);
K = zeros(N,1);
Gamma_r = (-2./r).';

if puncture == 1
    % punctured Schwarzchild BH ICs
    % tune r0 near 2 but not above.
    chi = cap(1.8,r,M,N);
end

% radial derivatives
alpha_p = f_prime(alpha,h,1,1,1,1,N);
alpha_pp = f_pprime(alpha,h,1,1,1,1,N);
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
Gamma_r_p = f_prime(Gamma_r,h,2/(r_min-3*h/2),2/(r_min-h/2),-2/(r_max+h/2),-2/(r_max+3*h/2),N);

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

function plot_RHS(n_iter,r_min,r_max,puncture)
L=zeros(1,10);
H=zeros(1,10);
for j=1:n_iter
    h=0.5^j;
    N = round((r_max-r_min)/h);
    v_old = initial_cond(h,r_min,r_max,puncture);
    U = dydt(1, v_old, h, N);
    U=reshape(U.',[],12);
    Uj=U(:,7);
    L(1,j)=sqrt(mean((Uj).^2));
    H(1,j)=h;
end
%plot loglog plot
loglog(H,L)
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