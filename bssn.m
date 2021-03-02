function linadv

    % Winter 2021
    % Assignment C1

    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the advection speed
    v=16;
    x_0=0;
    % number of grid points in space (not counting the point at x=1, which is the same as the point at x=0, due to the periodicity)
    x_max=8;
    x_min=-8;
    % spatial step
    h =0.1;
    N=(x_max-x_min)/h;
    % safety constant for stability (should be smaller than 1)
    cfl=.8;
    % CFl timestep limit for the explicit methods
    x=(x_min:h:x_max-h);
    %dt=cfl*h/(4/h^2+6*abs(max(one_soliton_IC(x))));
    dt=cfl*h/(4/h^2+6*abs(max(two_soliton_IC(x))));
    %dt=cfl*h/(4/h^2+6*abs(max(Gauss_IC(x))));
    % time at which we want to end the simulation
    t_end=1;
    % number of timesteps to be taken
    n_it=t_end/dt;
    % number of different methods we want to try
    n_methods=2;

    % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spatial grid (not including the point at x=1, which is the same as the point at x=0, due to the periodicity)
    x=(x_min:h:x_max-h);
    % temporal grid
    t=(0:dt:n_it*dt);
    % arrays for the numerical approximations at times n and n+1
    v_new=zeros(N,n_methods);
    v_old=zeros(N,n_methods);
    v_oldold=zeros(N,n_methods);
    % array for the exact solution
    v_exact=zeros(N,1);

    % the initial condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:n_methods
       % One soliton IC
       v_old(:,i)=-8./(cosh(2*x)).^2;
       % Two soliton IC
       %v_old(:,i)=two_soliton_IC(x);
       % Gaussian initial condition
       %v_old(:,i)=Gauss_IC(x);
    end


    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:n_it

        % method 1: LF
        % This is a 3 level problem, need to explicitly compute first
        % time step first
        
        if iter == 1
           v_new(3:N-2,1)= v_old(3:N-2,1) + 6/3*dt/h*(v_old(2:N-3,1)+v_old(3:N-2,1)+v_old(4:N-1,1))...
                         .*(v_old(4:N-1,1)-v_old(2:N-3,1))-dt/h^3*(v_old(5:N,1)-2*v_old(4:N-1,1)+2*v_old(2:N-3,1)-v_old(1:N-4,1));
            v_new(1,1) = v_old(1,1) + 6/3*dt/h*(v_old(N,1)+v_old(1,1)+v_old(2,1))...
                         .*(v_old(2,1)-v_old(N,1))-dt/h^3*(v_old(3,1)-2*v_old(2,1)+2*v_old(N,1)-v_old(N-1,1));
            v_new(2,1) = v_old(2,1) + 6/3*dt/h*(v_old(1,1)+v_old(2,1)+v_old(3,1))...
                         .*(v_old(3,1)-v_old(1,1))-dt/h^3*(v_old(4,1)-2*v_old(3,1)+2*v_old(1,1)-v_old(N,1))      ;      
            v_new(N-1,1) = v_old(N-1,1) + 6/3*dt/h*(v_old(N-2,1)+v_old(N-1,1)+v_old(N,1))...
                         .*(v_old(N,1)-v_old(N-2,1))-dt/h^3*(v_old(1,1)-2*v_old(N,1)+2*v_old(N-2,1)-v_old(N-3,1))  ;
            v_new(N,1) = v_old(N,1) + 6/3*dt/h*(v_old(N-1,1)+v_old(N,1)+v_old(1,1))...
                         .*(v_old(1,1)-v_old(N-1,1))-dt/h^3*(v_old(2,1)-2*v_old(1,1)+2*v_old(N-1,1)-v_old(N-2,1)) ;   
        else    
            v_new(3:N-2,1) = v_oldold(3:N-2,1) + 6/3*dt/h*(v_old(2:N-3,1)+v_old(3:N-2,1)+v_old(4:N-1,1))...
                             .*(v_old(4:N-1,1)-v_old(2:N-3,1))-dt/h^3*(v_old(5:N,1)-2*v_old(4:N-1,1)+2*v_old(2:N-3,1)-v_old(1:N-4,1));
            v_new(1,1) = v_oldold(1,1) + 6/3*dt/h*(v_old(N,1)+v_old(1,1)+v_old(2,1))...
                            .*(v_old(2,1)-v_old(N,1))-dt/h^3*(v_old(3,1)-2*v_old(2,1)+2*v_old(N,1)-v_old(N-1,1));
            v_new(2,1) = v_oldold(2,1) + 6/3*dt/h*(v_old(1,1)+v_old(2,1)+v_old(3,1))...
                            .*(v_old(3,1)-v_old(1,1))-dt/h^3*(v_old(4,1)-2*v_old(3,1)+2*v_old(1,1)-v_old(N,1))      ;      
            v_new(N-1,1) = v_oldold(N-1,1) + 6/3*dt/h*(v_old(N-2,1)+v_old(N-1,1)+v_old(N,1))...
                            .*(v_old(N,1)-v_old(N-2,1))-dt/h^3*(v_old(1,1)-2*v_old(N,1)+2*v_old(N-2,1)-v_old(N-3,1))  ;
            v_new(N,1) = v_oldold(N,1) + 6/3*dt/h*(v_old(N-1,1)+v_old(N,1)+v_old(1,1))...
                            .*(v_old(1,1)-v_old(N-1,1))-dt/h^3*(v_old(2,1)-2*v_old(1,1)+2*v_old(N-1,1)-v_old(N-2,1)) ;             
        end
        
        % method 2: Goda - calculate A matrix each loop
        [A,B]=mat_linadv_Goda(N,h,dt,v_old(:,2));
        v_new(:,2)=A\(B*v_old(:,2));
        
        % Initial condition #1
        %v_exact(:) = -v./(2* (cosh(sqrt(v)/2*(x-v*t(iter+1)-x_0))).^2 )...
        %             -v./(2* (cosh(sqrt(v)/2*(x-v*t(iter+1)+16))).^2 )...
        %             -v./(2* (cosh(sqrt(v)/2*(x-v*t(iter+1)+32))).^2 );
        v_exact(:) = -v./(2* (cosh(sqrt(v)/2*(x-v*t(iter+1)-x_0))).^2 )...
        %             -v./(2* (cosh(sqrt(v)/2*(x-v*t(iter+1)+16))).^2 )...
        %             -v./(2* (cosh(sqrt(v)/2*(x-v*t(iter+1)+32))).^2 );
        

                
        % graphical output
        %plot(x,v_exact(:),'^r-')
        plot(x,v_new(:,1),'*b-')
        hold on
        %plot(x,v_new(:,1),'*b-')
        plot(x,v_new(:,2),'+g-')
        axis([-8 8 -9 1.2])
        xlabel('x')
        ylabel('v')
        title('KdV simulation')
        pause(0.001)
        hold off
        
		%legend('Exact','LF', 'Goda')
        
        % prepare for the next iteration
        v_oldold=v_old;
        v_old=v_new;
    end
    legend('LF', 'Goda')
end

function y=one_soliton_IC(x)
    y=-8./((cosh(2*x)).^2);
end

function y=two_soliton_IC(x)
    y=-8./((cosh(2*x+8)).^2)-2./((cosh(x)).^2);
end

function y=Gauss_IC(x)
    y=-5*exp(-((x/2).^2));
end

% the matrices for the implicit method
%-------------------------
% implicit Goda scheme
function [A,B]=mat_linadv_Goda(N,h,dt,v_old)
  e = ones(N,1); 
  A = spdiags([-0.5*dt/h^3*e dt/h^3*e+dt/h*([v_old(2:N); v_old(1)]+v_old)...
               e -dt/h^3*e-dt/h*([v_old(N); v_old(1:N-1)]+v_old) ...
               0.5*dt/h^3*e], [-2:2], N, N);
  B=speye(N);
  % Necessary boundary conditions for first and last couple rows
  A(1,N) = dt/h^3+dt/h*(v_old(1)+v_old(N));
  A(1,N-1) = -0.5*dt/h^3;
  A(2,N) = -0.5*dt/h^3;
  A(N-1,1) = 0.5*dt/h^3;
  A(N,1) = -dt/h^3-dt/h*(v_old(N)+v_old(1));
  A(N,2) = 0.5*dt/h^3;
end