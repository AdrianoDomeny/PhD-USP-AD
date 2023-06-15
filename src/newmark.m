
% -----------------------------------------------------------------
%  newmark__ode_solver.m
%  
%  This function computes the time response of the following
%  initial value problem:
%  
%  [M]*Qacce + [C]*Qvelo + [K]*Qdisp = F + FNL   for t > 0
%                              Qdisp = Qdisp0    for t = 0
%                              Qvelo = Qvelo0    for t = 0
%  where
%  
%  [M]    - mass      matrix
%  [C]    - damping   matrix
%  [K]    - stiffness matrix
%  
%  F      - linear    force vector in R^{Ndofs} at time t
%  FNL    - nonlinear force vector in R^{Ndofs} at time t
%  
%  Qacce  - acceleration vector in R^{Ndofs} at time t
%  Qvelo  - velocity     vector in R^{Ndofs} at time t
%  Qdisp  - displacement vector in R^{Ndofs} at time t
%  
%  Qvelo0 - initial velocity     vector in R^{Ndofs} at time t
%  Qdisp0 - initial displacement vector in R^{Ndofs} at time t
%  
%  
%  Algorithm:
%  This function uses an implicit Newmark sheme for time integration 
%  of the equation of motion and a fixed point iteration, combined 
%  with the Successive Over-Relaxation (SOR) technique, to solve the 
%  nonlinear system of algebraic equations associated.
%  
%  
%  References:
%  M. Geradin and A. Cardona
%  Flexible Multibody Dynamics: A Finite Element Approach
%  John Willey & Sons, 2001, page 29
%  
%  T. J. R. Hughes
%  The Finite Element Method: 
%  Linear Static and Dynamic Finite Element Analysis 
%  Dover Civil and Mechanical Engineering, 2000, page 490
%  
%  C. Soize
%  Dynamique des structures, elements de base et
%  concepts fondamentaux
%  Ellipses Marketing, 2001, pages 55--56
%  
%  N. Newmark
%  A method of computation for structural dynamics
%  Journal of the Engineering Mechanics Division 
%  vol. 85 pp. 67--94, 1959
%  
%  
%  Input:
%  Nelem         - number of elements
%  Ndofs         - number of degress of freedon
%  Nldofs        - number of local degrees of freedon
%  Nbc           - number of geometric boundary conditions
%  Ndt           - number of time steps
%  Nx            - number of points in the spatial mesh
%  Ngp           - number of points for Gaussian quadrature
%  dt            - time step
%  M             - (Ndofs  x Ndofs) mass matrix
%  C             - (Ndofs  x Ndofs) damping matrix
%  K             - (Ndofs  x Ndofs) stiffness matrix
%  F             - (Ndofs  x 1    ) linear force vector
%  Qdisp0        - (Ndofs  x 1    ) generalized initial displacment
%  Qvelo0        - (Ndofs  x 1    ) generalized initial velocity
%  B             - (Nbc    x Ndofs) constraints matrix
%  h             - (Ndofs  x 1    ) constraints vector
%  h_dot         - (Ndofs  x 1    ) constraints vector time derivative
%  gauss_points  - Gauss points vector
%  gauss_weights - Gauss weights vector
%  shape_tab     - shape functions table
%                  (Nldfos*Ngp x Ndofs_node*Nelem)
%  dshape_tab    - shape functions derivatives table
%                  (Nldfos*Ngp x Ndofs_node*Nelem)
%  GlobalDoF     - (Nldofs x Nelem) global DoF identification matrix
%  xmesh         - (Nx  x 1) spatial mesh vector
%  time          - (Ndt x 1) temporal mesh vector
%  PHI           - (Ndofs x Nred) modal projection matrix
%  phys_param    - physical parameters vector
%  num_param     - numerical parameters vector
%  
%  Output:
%  Qdisp      - (Ndofs x Ndt) generalized displacement
%  Qvelo      - (Ndofs x Ndt) generalized velocity
%  Qacce      - (Ndofs x Ndt) generalized accelaration
%  lambda     - (Nbc   x Ndt) Lagrande multipliers matrix
%  error_flag - error flag
%  res_vec    - residual vector
%  iter_vec   - iteration counter vector
%  time_flag  - time advance flag
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Sep 22, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [Qdisp,Qvelo,Qacce,lambda,...
          error_flag,res_vec,iter_vec,time_flag] = ...
          newmark(Nelem,Ndofs,Nldofs,...
                              Nbc,Ndt,Nx,Ngp,dt,...
                              M,C,K,F,Qdisp0,Qvelo0,...
                              B,h,h_dot,gauss_points,gauss_weights,...
                              shape_tab,dshape_tab,GlobalDoF,...
                              xmesh,time,PHI,...
                              phys_param,num_param)
	
    % error flag
	error_flag = 1;
	
    % check number of arguments
    if nargin < 26
        error(' Too few inputs.')
    elseif nargin > 27
        error(' Too many inputs.')
    end
    
    % check input arguments
    if Ndofs <= 0
        error('Ndofs must be a positive integer.')
    end
    
    if Nbc <= 0
        error('Nbc must be a positive integer.')
    end
    
    if Ndt <= 0
        error('Ndt must be a positive integer.')
    end
    
    if dt <= 0.0
        error('time step must be positive.')
    end
    
    
    % convert time to a column vector (if necessary)
    if find( size(time) == max(size(time)) ) > 1
        time = time';
    end
    
    
    % define numerical parameters
    if nargin == 27
        
        % maximum number of iterations
        max_iter = num_param(1);
        
        if max_iter <= 0
            error('max_iter must be a positive integer.')
        end
        
        % time step refinement parameter
        Ndt_ref = num_param(2);
        
        if Ndt_ref <= 1
            error('Ndt_ref must be integer greater than 1.')
        end
        
        % tolerance
        tol = num_param(3);
        
        if tol <= 0.0
            error('tol must be a positive.')
        end
        
        % SOR relaxation parameter
        w_sor = num_param(4);
        
        if w_sor <= 0.0 || w_sor >= 2.0
            error('w_sor must geater than 0 and less than 2.')
        end
        
        % reduced model dimension
        Nred = num_param(5);
        
        if Nred <= 0 || Nred > Ndofs
            error('Nred must be a positie integer less than or equal to Ndofs.')
        end
        
        % number of iteration before show current time on the screen
        show_time_flag = num_param(6);
        
	else
        
        % maximum number of iterations
        max_iter = 20;
        
        % time step refinement parameter
        Ndt_ref = 10;
        
        % tolerance
        tol = 1.0e-3;
        
        % SOR relaxation parameter
        w_sor = 1.0;
        
        % reduced model dimension
        Nred = Ndofs;
        
        % number of iteration before show current time on the screen
        show_time_flag = 0;

    end
    
    
    % preallocate memory for residual vector
    res_vec = zeros(Ndt,1);
    
	
    % preallocate memory for iteration counter vector
    iter_vec = zeros(Ndt,1);
    
    
    % preallocate memory for displacement
    Qdisp = zeros(Nred,Ndt);
    
    
    % preallocate memory for velocity
    Qvelo = zeros(Nred,Ndt);
    
    
    % preallocate memory for accelation
    Qacce = zeros(Nred,Ndt);
    
    
    % preallocate memory for lagrange multipliers
    lambda = zeros(Nbc,Ndt);
        
    
    % The Newmark scheme has two parameters, gamma and beta, and they
    % control its stability and accuracy characteristics. The method 
    % is unconditionnaly stable if
    %
    %          gamma >= 0.5 and beta >= 0.25(gamma + 0.5)^2
    %
    % If gamma = 0.5, then the method is stable, but there is no 
    % damping of the high frequencies introduced by the discretization
    %
    % If gamma > 0.5, then the high frequencies are damped. This damping
    % is introduced through a parameter alpha > 0.
    
    %alpha = 0.3;
    %alpha = 0.1;
    %alpha = 0.01;
    %alpha = 0.001;
    %alpha = 0.0001;
    alpha = 0.0;
    gamma = 0.5 + alpha;
    beta  = 0.25*(0.5 + gamma)^2;
    
    
    % Newmark method constants
    a0 = 1/(beta*dt*dt);
    a1 = gamma/(beta*dt);
    a2 = 1/(beta*dt);
    a3 = 1/(2*beta) - 1;
    a4 = gamma/beta - 1;
    a5 = 0.5*dt*(gamma/beta - 2);
    a6 = dt*(1 - gamma);
    a7 = gamma*dt;
    a8 = a0/a1;
    a9 = a0*(a4/a1);
    
    
    % caractheristic time
    %tau_st = phys_param(52);
    
    
    % save operations for future use
    one_minus_w_sor = 1.0 - w_sor;
    
    
    % compute the effective stiffness matrix
    Keff = a0*M + a1*C + K;
    
    
    % scale the matrix B
    a0B = a0*B;
    
    % transpose matrix of the matrix Ba0
    a0BT = a0B';
    
    
    % Cholesky decomposition of Keff = LKeff*LKeff^T
    [LKeff,Keff_posdef_flag] = chol(Keff,'lower');
    
    
    % check Cholesky decomposition of Keff for error
    if Keff_posdef_flag == 0
        
        
        % define uppper triangular matrix  for Keff
        UKeff = LKeff';
        
        
        % constrained stiffness matrix
        Kbar = a0B*(UKeff\(LKeff\a0BT));
       
        
        % Cholesky decomposition of Kbar = LKbar*LKbar^T
        [LKbar,Kbar_posdef_flag] = chol(Kbar,'lower');
        
        
        % check the Cholesky decomposition of Kbar for error
        if Kbar_posdef_flag == 0
            
            % define uppper triangular matrix for Kbar
            UKbar = LKbar';
        else
        
            % LU decomposition of Kbar = LKbar*UKbar^T
            [LKbar,UKbar] = lu(Kbar);
        end
        
    else
        
        % define extended stiffness matrix
        Khat = [Keff a0BT; a0B zeros(Nbc,Nbc)];
        
        
        % LU decomposition of Khat = LKhat*UKhat^T
        [LKhat,UKhat] = lu(Khat);
    end
    
    
	% time t = t0
	t0 = time(1);
	
    
	% time t = t1
	%t1 = time(2);
    
    
    % f(t) = 1
    % temporal function to be multiplied by F
    f_t0 = sin(t0);
    
    
    % compute the linear force at t = t0
    FL_t0 = M*(a0*Qdisp0 + a2*Qvelo0 + a3*zeros(Nred,1)) ...
          + C*(a1*Qdisp0 + a4*Qvelo0 + a5*zeros(Nred,1)) ...
          + f_t0*F;
    
	
    % compute the nonlinear terms force at t = t0
    %FNL_t0 = fem__assembly_nonlinear(Nelem,Ndofs,Nldofs,Nx,Ngp,...
    %                                 Qdisp0,Qvelo0,zeros(Nred,1),...
    %                                 gauss_points,gauss_weights,...
    %                                 shape_tab,dshape_tab,...
    %                                 GlobalDoF,xmesh,t0,...
    %                                 phys_param,PHI);
	
	
	% compute the effective force at t = t0
	Feff_t0 = FL_t0;% + FNL_t0;
    
    
    % g(t) = t and g_dot(t) = 1
    % temporal functions to be multiplied by h and h_dot
    %g_t0 = t0;
    %g_t1 = t1;
    
    %g_dot_t0 = 1.0;
    %g_dot_t1 = 1.0;
    
    
	% compute the effective constraint vector t = t0
	%a0Heff_t0 = a8*h_dot(:,1)*g_dot_t1 + ...
    %            a9*h_dot(:,1)*g_t0 + ...
    %            a0*h(:,1)*g_t0;
    a0Heff_t0 = a0*h(:,1);
	
    
    % compute the Lagrange multiplier at t = t0
    if Keff_posdef_flag == 0
        
        % Lagrange multipliers at t = t0
        lambda0 = UKbar\(LKbar\(a0B*(UKeff\(LKeff\Feff_t0)) - a0Heff_t0));

    else

        % define the extended effective force at t = t0
        Fhat_t0 = [Feff_t0; a0Heff_t0];

        % compute the extended displacement at t = t0
        Qhat_t0 = UKhat\(LKhat\Fhat_t0);
        
        % Lagrange multipliers at t = t0
        lambda0 = Qhat_t0(Nred+1:Nred+Nbc,1);
    end
    
    
    % procedure to avoid lambda0 be infinity or NaN
    lambda0(isnan(lambda0)) = 0.0;
    lambda0(isinf(lambda0)) = 0.0;
    
    
	% compute the acceleration at t = t0
	Qacce0 = M\(f_t0*F(:,1) + FNL_t0 - C*Qvelo0 - K*Qdisp0 - B'*lambda0);
    
    
    % procedure to avoid Qacce0 be infinity or NaN
    Qacce0(isnan(Qacce0)) = 0.0;
    Qacce0(isinf(Qacce0)) = 0.0;
    
    
    % set the initial conditions
     Qdisp(:,1) =  Qdisp0;
     Qvelo(:,1) =  Qvelo0;
     Qacce(:,1) =  Qacce0;
    lambda(:,1) = lambda0;

    
    % loop on time-samplings
    for n=1:(Ndt-1)
        
        % initialize the generalized vectors at t = t_{n}
         Qdisp_tn =  Qdisp(:,n);
         Qvelo_tn =  Qvelo(:,n);
         Qacce_tn =  Qacce(:,n);
        lambda_tn = lambda(:,n);
        
        
        % initialize the generalized vectors at t = t_{n+1}
         Qdisp_tn1 =  Qdisp_tn;
         Qvelo_tn1 =  Qvelo_tn;
         Qacce_tn1 =  Qacce_tn;
        lambda_tn1 = lambda_tn;
        
        
        % time t = t_{n}
        tn = time(n);
        
        
        % time t = t_{n+1}
        tn1 = time(n+1);
        
        
        % f(t) = 1
        % temporal function at t = t_{n+1} to be multiplied by F
        %f_tn1 = min(tn1/tau_st,1.0);
        f_tn1 = sin(tn1);
        
        
        % compute the linear force at t = t_{n+1}
        FL_tn1 = M*(a0*Qdisp_tn + a2*Qvelo_tn + a3*Qacce_tn) ...
               + C*(a1*Qdisp_tn + a4*Qvelo_tn + a5*Qacce_tn) ...
               + f_tn1*F;
        
        
        % compute the nonlinear force at t = t_{n+1}
        %FNL_tn1 = fem__assembly_nonlinear(Nelem,Ndofs,Nldofs,Nx,Ngp,...
        %                                  Qdisp_tn1,Qvelo_tn1,Qacce_tn1,...
        %                                  gauss_points,gauss_weights,...
        %                                  shape_tab,dshape_tab,...
        %                                  GlobalDoF,xmesh,tn1,...
        %                                  phys_param,PHI);
        
        
        % compute the effective force at t = t_{n+1}
        Feff_tn1 = FL_tn1;% + FNL_tn1;
        
        
        % g(t) = t and g_dot(t) = 1
        % temporal functions to be multiplied by h and h_dot
        %g_tn  = tn;
        %g_tn1 = tn1;
        
        %g_dot_tn = 1.0;
        %g_dot_tn1 = 1.0;
        
        
        % compute the effective constraint vector t = t_{n+1}
        %a0Heff_tn1 = a8*h_dot(:,1)*g_dot_tn1 + ...
        %             a9*h_dot(:,1)*g_tn + ...
        %             a0*h(:,1)*g_tn;
        a0Heff_tn1 = a0*h(:,1);
        
        
        % initial residual at t = t_{n+1}
        res_vec(n+1) = tol + 1.0;
        
        
        % every "show_time_flag" iterations plot the time on the screen
        if mod(n,show_time_flag) == 0
            disp(['time = ',num2str(tn)])
            disp(' ')
        end
        
        
        % fixed point iteration
        while res_vec(n+1) > tol && iter_vec(n+1) < max_iter
            
            
            % update iteration counter
            iter_vec(n+1) = iter_vec(n+1) + 1;
            
            
            % define the old extended displacement at t = t_{n+1}
            Qhat_tn1_old = [Qdisp_tn1; lambda_tn1];
            
            
            % solve the nonlinear system of algebraic equations
            if Keff_posdef_flag == 0
                
                % compute the new Lagrange multipliers at t = t_{n+1}
                lambda_tn1 = ...
                  UKbar\(LKbar\(a0B*(UKeff\(LKeff\Feff_tn1))-a0Heff_tn1));
                
                
                % compute the new generalized displacement at t = t_{n+1}
                Qdisp_tn1 = UKeff\(LKeff\(Feff_tn1 - a0BT*lambda_tn1));
                
                
                % define the new extended displacement at t = t_{n+1}
                Qhat_tn1 = [Qdisp_tn1; lambda_tn1];

            else
                
                % define the extended effective force at t = t_{n+1}
                Fhat = [Feff_tn1; a0Heff_tn1];
                
                
                % compute the new extended displacement at t = t_{n+1}
                Qhat_tn1 = UKhat\(LKhat\Fhat);
            end
            
            
            % successive over relaxation
            Qhat_tn1 = w_sor*Qhat_tn1 + one_minus_w_sor*Qhat_tn1_old;
            
            
            % generalized displacement at t = t_{n+1}
            Qdisp_tn1 = Qhat_tn1(1:Nred,1);
            
            
            % Lagrange multipliers at t = t_{n+1}
            lambda_tn1 = Qhat_tn1(Nred+1:Nred+Nbc,1);
            
            
            % compute the new generalized acceleration at t = t_{n+1}
            Qacce_tn1 =   a0*(Qdisp_tn1 - Qdisp_tn) ...
                        - a2*Qvelo_tn ...
                        - a3*Qacce_tn;
            
            
            % compute the new generalized velocity at t = t_{n+1}
            Qvelo_tn1 = Qvelo_tn + a6*Qacce_tn + a7*Qacce_tn1;
            
            
            % compute the new nonlinear force at t = t_{n+1}
            FNL_tn1 = fem__assembly_nonlinear(Nelem,Ndofs,Nldofs,Nx,Ngp,...
                                              Qdisp_tn1,Qvelo_tn1,Qacce_tn1,...
                                              gauss_points,gauss_weights,...
                                              shape_tab,dshape_tab,...
                                              GlobalDoF,xmesh,tn1,...
                                              phys_param,PHI);
            
            
            % compute the new effective force at t = t_{n+1}
            Feff_tn1 = FL_tn1 + FNL_tn1;
            
            
            % compute the new residual at t = t_{n+1}
            res_num      =     norm(Qhat_tn1 - Qhat_tn1_old);
            res_den      = 0.5*norm(Qhat_tn1 + Qhat_tn1_old);
            res_vec(n+1) = res_num/res_den;
        end
        
        
        % time advance flag
        time_flag = n + 1;
        
            
        % check the convergence
        if res_vec(n+1) < tol
            
            % update generalized vectors at t = t_{n+1}
             Qdisp(:,n+1) =  Qdisp_tn1;
             Qvelo(:,n+1) =  Qvelo_tn1;
             Qacce(:,n+1) =  Qacce_tn1;
            lambda(:,n+1) = lambda_tn1;
            
            % error flag
            error_flag = 0;
        else
            
            % define the refined time step
            dt_ref = (tn1-tn)/(Ndt_ref-1);
            
            
            % define refined temporal mesh
            time_ref = (tn:dt_ref:tn1)';
            

            % define refined constaints vectors
            %h_ref     = (h(:,n+1)/tn1)*(time_ref');
            %h_dot_ref = h_dot(:,n+1)*ones(1,Ndt_ref);
            
            
            % recursive call to Newmark routine
            [Qdisp_ref,Qvelo_ref,Qacce_ref,lambda_ref,...
                                 error_flag,res_vec_ref] = ...
             newmark(Nelem,Ndofs,Nldofs,...
                                 Nbc,Ndt_ref,Nx,Ngp,dt_ref,...
                                 M,C,K,F,Qdisp(:,n),Qvelo(:,n),...
                                 B,h,h_dot,...
                                 gauss_points,gauss_weights,...
                                 shape_tab,dshape_tab,GlobalDoF,...
                                 xmesh,time_ref,PHI,...
                                 phys_param,num_param);
            
            
            % check if refined Newmark integration converged
            if error_flag == 0
                
                % update generalized vectors at t = t_{n+1}
                 Qdisp(:,n+1) =  Qdisp_ref(:,Ndt_ref);
                 Qvelo(:,n+1) =  Qvelo_ref(:,Ndt_ref);
                 Qacce(:,n+1) =  Qacce_ref(:,Ndt_ref);
                lambda(:,n+1) = lambda_ref(:,Ndt_ref);
                
                % update residual at t = t_{n+1}
                res_vec(n+1) = res_vec_ref(Ndt_ref);
            else
                
                % break the loop
                break
            end
        end
    end
    
return
% -----------------------------------------------------------------
