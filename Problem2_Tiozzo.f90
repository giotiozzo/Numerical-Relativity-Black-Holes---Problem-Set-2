MODULE ADM_1plus1
IMPLICIT NONE
    
CONTAINS
    !Computes spatial derivative using second-order centered finite differences with PERIODIC boundary conditions
    !Input: u(N)=function to differentiate, N=number of grid points, dx=grid spacing
    !Output: du(N)=derivative du/dx
    
    SUBROUTINE derivative_periodic(u,du,N,dx)
        INTEGER,INTENT(IN)::N
        REAL*8,INTENT(IN)::u(N),dx
        REAL*8,INTENT(OUT)::du(N)
        INTEGER::i
        
        !Interior points: standard centered difference
        DO i=2,N-1
            du(i)=(u(i+1)-u(i-1))/(2.d0*dx)
        END DO
        
        ! Periodic boundary conditions
        ! Left boundary: wraps to right
        du(1)=(u(2)-u(N))/(2.d0*dx)
        
        ! Right boundary: wraps to left
        du(N)=(u(1)-u(N-1))/(2.d0*dx)
        
    END SUBROUTINE derivative_periodic
    
    !Computes logarithmic derivatives: D_g=dx ln(g),D_alpha=dx ln(alpha)
    
    SUBROUTINE compute_Dg_Dalpha(g,alpha,Dg,Dalpha,N,dx)
        INTEGER,INTENT(IN)::N
        REAL*8,INTENT(IN)::g(N),alpha(N),dx
        REAL*8,INTENT(OUT)::Dg(N),Dalpha(N)
        REAL*8::g_deriv(N),alpha_deriv(N)
        INTEGER::i
        
        !Compute derivatives
        CALL derivative_periodic(g,g_deriv,N,dx)
        CALL derivative_periodic(alpha,alpha_deriv,N,dx)
        
        !Compute logarithmic derivatives
        DO i=1,N
            Dg(i)=g_deriv(i)/g(i)
            Dalpha(i)=alpha_deriv(i)/alpha(i)
        END DO
        
    END SUBROUTINE compute_Dg_Dalpha
    

    !Computes initial K_tilde=SQRT g K from second derivative of h(x)
    !Formula: K_xx=-h''/SQRTg,so K_tilde=SQRTg*K=-h''
    
    SUBROUTINE compute_initial_K(h,K_tilde,N,dx)
        INTEGER,INTENT(IN)::N
        REAL*8,INTENT(IN)::h(N),dx
        REAL*8,INTENT(OUT)::K_tilde(N)
        REAL*8::hp(N),hpp(N)
        INTEGER::i
        
        !First derivative
        CALL derivative_periodic(h,hp,N,dx)
        
        !Second derivative
        CALL derivative_periodic(hp,hpp,N,dx)
        
        DO i=1,N
            K_tilde(i)=-hpp(i)
        END DO
        
    END SUBROUTINE compute_initial_K

    !Computes right-hand side of evolution equations
    
    SUBROUTINE compute_rhs(g,K_tilde,alpha,Dalpha,&
                          rhs_g,rhs_Dg,rhs_K,rhs_alpha,rhs_Dalpha,N,dx,f_slicing)
        INTEGER,INTENT(IN)::N
        REAL*8,INTENT(IN)::g(N),K_tilde(N),alpha(N),Dalpha(N),dx,f_slicing
        REAL*8,INTENT(OUT)::rhs_g(N),rhs_Dg(N),rhs_K(N),rhs_alpha(N),rhs_Dalpha(N)
        REAL*8::K(N),temp(N),temp_deriv(N)
        INTEGER::i
        
        ! Compute K from K_tilde
        DO i=1,N
            K(i)=K_tilde(i)/SQRT(g(i))
        END DO
        
        ! Evolution equation for g: dg/dt = -2*alpha*g*K
        DO i=1,N
            rhs_g(i)=-2.d0*alpha(i)*g(i)*K(i)
        END DO
        
        ! Evolution equation for D_g: dD_g/dt = -d/dx(2*alpha*K)
        DO i=1,N
            temp(i)=2.d0*alpha(i)*K(i)
        END DO
        CALL derivative_periodic(temp,temp_deriv,N,dx)
        DO i=1,N
            rhs_Dg(i)=-temp_deriv(i)
        END DO
        
        ! Evolution equation for K_tilde: dK_tilde/dt = -d/dx(alpha*D_alpha*g^(-1/2))
        DO i=1,N
            temp(i)=alpha(i)*Dalpha(i)/SQRT(g(i))
        END DO
        CALL derivative_periodic(temp,temp_deriv,N,dx)
        DO i=1,N
            rhs_K(i)=-temp_deriv(i)
        END DO
        
        ! Evolution equation for alpha: dalpha/dt = -alpha^2*f*K
        DO i=1,N
            rhs_alpha(i)=-alpha(i)**2*f_slicing*K(i)
        END DO
        
        ! Evolution equation for D_alpha: dD_alpha/dt = -d/dx(alpha*f*K)
        DO i=1,N
            temp(i)=alpha(i)*f_slicing*K(i)
        END DO
        CALL derivative_periodic(temp,temp_deriv,N,dx)
        DO i=1,N
            rhs_Dalpha(i)=-temp_deriv(i)
        END DO
        
    END SUBROUTINE compute_rhs
    
    !Performs one time step using fourth-order Runge-Kutta method
    
    SUBROUTINE rk4_step(g,Dg,K_tilde,alpha,Dalpha,N,dx,dt,f_slicing)
        INTEGER,INTENT(IN)::N
        REAL*8,INTENT(INOUT)::g(N),Dg(N),K_tilde(N),alpha(N),Dalpha(N)
        REAL*8,INTENT(IN)::dx,dt,f_slicing
        
        !RK4 intermediate variables
        REAL*8::k1_g(N),k1_Dg(N),k1_K(N),k1_alpha(N),k1_Dalpha(N)
        REAL*8::k2_g(N),k2_Dg(N),k2_K(N),k2_alpha(N),k2_Dalpha(N)
        REAL*8::k3_g(N),k3_Dg(N),k3_K(N),k3_alpha(N),k3_Dalpha(N)
        REAL*8::k4_g(N),k4_Dg(N),k4_K(N),k4_alpha(N),k4_Dalpha(N)
        REAL*8::g_temp(N),Dg_temp(N),K_temp(N),alpha_temp(N),Dalpha_temp(N)
        
        !Stage 1: k1=f(y)
        CALL compute_rhs(g,K_tilde,alpha,Dalpha,&
                        k1_g,k1_Dg,k1_K,k1_alpha,k1_Dalpha,N,dx,f_slicing)
        
        !Stage 2: k2=f(y+dt*k1/2)
        g_temp=g+0.5d0*dt*k1_g
        Dg_temp=Dg+0.5d0*dt*k1_Dg
        K_temp=K_tilde+0.5d0*dt*k1_K
        alpha_temp=alpha+0.5d0*dt*k1_alpha
        Dalpha_temp=Dalpha+0.5d0*dt*k1_Dalpha
        
        CALL compute_rhs(g_temp,K_temp,alpha_temp,Dalpha_temp,&
                        k2_g,k2_Dg,k2_K,k2_alpha,k2_Dalpha,N,dx,f_slicing)
        
        !Stage 3: k3=f(y+dt*k2/2)
        g_temp=g+0.5d0*dt*k2_g
        Dg_temp=Dg+0.5d0*dt*k2_Dg
        K_temp=K_tilde+0.5d0*dt*k2_K
        alpha_temp=alpha+0.5d0*dt*k2_alpha
        Dalpha_temp=Dalpha+0.5d0*dt*k2_Dalpha
        
        CALL compute_rhs(g_temp,K_temp,alpha_temp,Dalpha_temp,&
                        k3_g,k3_Dg,k3_K,k3_alpha,k3_Dalpha,N,dx,f_slicing)
        
        !Stage 4: k4=f(y+dt*k3)
        g_temp=g+dt*k3_g
        Dg_temp=Dg+dt*k3_Dg
        K_temp=K_tilde+dt*k3_K
        alpha_temp=alpha+dt*k3_alpha
        Dalpha_temp=Dalpha+dt*k3_Dalpha
        
        CALL compute_rhs(g_temp,K_temp,alpha_temp,Dalpha_temp,&
                        k4_g,k4_Dg,k4_K,k4_alpha,k4_Dalpha,N,dx,f_slicing)
        
        !Update: y_new=y+dt*(k1+2*k2+2*k3+k4)/6
        g=g+dt*(k1_g+2.d0*k2_g+2.d0*k3_g+k4_g)/6.d0
        Dg=Dg+dt*(k1_Dg+2.d0*k2_Dg+2.d0*k3_Dg+k4_Dg)/6.d0
        K_tilde=K_tilde+dt*(k1_K+2.d0*k2_K+2.d0*k3_K+k4_K)/6.d0
        alpha=alpha+dt*(k1_alpha+2.d0*k2_alpha+2.d0*k3_alpha+k4_alpha)/6.d0
        Dalpha=Dalpha+dt*(k1_Dalpha+2.d0*k2_Dalpha+2.d0*k3_Dalpha+k4_Dalpha)/6.d0
        
    END SUBROUTINE rk4_step
    
END MODULE ADM_1plus1

!Main Program

PROGRAM Main
    USE ADM_1plus1
    IMPLICIT NONE
    
    !Simulation parameters (CORRECTED to match problem requirements)
    INTEGER,PARAMETER::N=400              !Number of spatial grid points
    REAL*8,PARAMETER::x_min=-10.d0        !Left boundary
    REAL*8,PARAMETER::x_max=40.d0         !Right boundary
    REAL*8,PARAMETER::t_final=10.d0       !Final simulation time
    REAL*8,PARAMETER::CFL=0.05d0          !CFL stability factor
    REAL*8,PARAMETER::f_slicing=1.d0      !Bona-Masso parameter (f=1)
    
    !Grid and time stepping
    REAL*8::dx,dt,t
    REAL*8::x(N)
    INTEGER::i,step,max_steps
    
    !ADM variables
    REAL*8::g(N),Dg(N),K_tilde(N),alpha(N),Dalpha(N)
    
    !Initial data: h(x)=0.1*exp(-x^2/6) as specified in problem
    REAL*8::h(N),hp(N)
    REAL*8,PARAMETER::h0=0.1d0  !Amplitude
    
    !Output storage, save at t=0,2.5,5,7.5,10
    REAL*8::alpha_data(N,5)
    REAL*8::g_data(N,5)
    REAL*8::K_data(N,5)
    REAL*8::save_times(5)
    INTEGER::save_counter
    INTEGER::output_unit
    REAL*8::t_output
    
    save_times=(/0.d0, 2.5d0, 5.d0, 7.5d0, 10.d0/)
    
    PRINT *,'Numerical Relativity: 1+1 ADM Evolution'
    PRINT *,'Gauge Wave Propagation with Periodic BC'
    PRINT *,''
    
    !Set up spatial grid
    dx=(x_max-x_min)/DBLE(N)
    DO i=1,N
        x(i)=x_min+(i-1)*dx
    END DO
    
    !Calculate time step from CFL condition
    dt=CFL*dx
    max_steps=INT(t_final/dt)+1
    
    PRINT *,'Grid parameters:'
    PRINT *,'Spatial points: N =',N
    PRINT *,'Domain: x in [',x_min,',',x_max,']'
    PRINT *,'Grid spacing: dx =',dx
    PRINT *,'Time step: dt =',dt
    PRINT *,'CFL number:',CFL
    PRINT *,'Total steps:',max_steps
    PRINT *,'Boundary conditions: PERIODIC'
    PRINT *,'Slicing: Bona-Masso with f =',f_slicing

    !INITIAL DATA: h(x)=0.1*exp(-x^2/6)
    PRINT *,'Initial data:'
    PRINT *,'Profile: h(x)=0.1*exp(-x^2/6)'
    
    !Calculate h(x)=0.1*exp(-x^2/6)
    DO i=1,N
        h(i)=h0*EXP(-x(i)**2/6.d0)
    END DO
    
    !Calculate h'(x) using periodic finite differences
    CALL derivative_periodic(h,hp,N,dx)
    
    !Initialize ADM variables from initial slicing
    !g=1-h'^2
    !alpha=sqrt((1-h')/(1+h'))
    DO i=1,N
        g(i)=1.d0-hp(i)**2
        alpha(i)=SQRT((1.d0-hp(i))/(1.d0+hp(i)))
    END DO
    
    !Calculate D_g=dxln(g) and D_alpha=dxln(alpha)
    CALL compute_Dg_Dalpha(g,alpha,Dg,Dalpha,N,dx)
    
    !Calculate K_tilde=sqrt(g)*K from second derivative of h
    CALL compute_initial_K(h,K_tilde,N,dx)
    
    PRINT *,'Initial data set successfully.'
    PRINT *,'max(alpha) =',MAXVAL(alpha)
    PRINT *,'min(alpha) =',MINVAL(alpha)
    PRINT *,'max(g) =',MAXVAL(g)
    PRINT *,'min(g) =',MINVAL(g)
    PRINT *,'max(K_tilde) =',MAXVAL(K_tilde)
    PRINT *,'min(K_tilde) =',MINVAL(K_tilde)
    
    !TIME EVOLUTION LOOP
    
    PRINT *,'Starting time evolution'
    
    t=0.d0
    
    !Store initial data (t=0)
    alpha_data(:,1)=alpha
    g_data(:,1)=g
    K_data(:,1)=K_tilde/SQRT(g)  ! Store K
    save_counter=2
    
    !Main evolution loop
    DO step=1,max_steps
        !Evolve one time step using RK4
        CALL rk4_step(g,Dg,K_tilde,alpha,Dalpha,N,dx,dt,f_slicing)
        
        t=t+dt
        
        !Progress indicator every 100 steps
        IF (MOD(step,100)==0) THEN
            PRINT *,'Progress: t =',t,'/',t_final,'  (step',step,')'
        END IF
        
        !Save data at specified times
        IF (save_counter<=5) THEN
            IF (ABS(t-save_times(save_counter))<dt/2.d0 .OR. &
                t>save_times(save_counter)) THEN
                alpha_data(:,save_counter)=alpha
                g_data(:,save_counter)=g
                K_data(:,save_counter)=K_tilde/SQRT(g)
                PRINT *,'Saved snapshot at t=',save_times(save_counter)
                save_counter=save_counter+1
            END IF
        END IF
        
        IF (t>=t_final) EXIT
    END DO
    
    PRINT *,'Time evolution completed successfully.'
    PRINT *,'Final time: t =',t
    PRINT *,'  Final max(alpha)=',MAXVAL(alpha)
    PRINT *,'  Final min(alpha)=',MINVAL(alpha)
    PRINT *,'  Final max(g)=',MAXVAL(g)
    PRINT *,'  Final min(g)=',MINVAL(g)
    
    !OUTPUT RESULTS
    
    PRINT *,'Writing output files'
    
    OPEN(NEWUNIT=output_unit,FILE='evolution_data.dat',STATUS='REPLACE')
    
    WRITE(output_unit,'(A)') '# 1+1 ADM Evolution with Periodic Boundary Conditions'
    WRITE(output_unit,'(A)') '# h(x)=0.1*exp(-x^2/6)'
    WRITE(output_unit,'(A)') '# Domain: x in [-10, 40], N=400 points'
    WRITE(output_unit,'(A)') '# CFL=0.05, f=1 (Bona-Masso slicing)'
    WRITE(output_unit,'(A)') '# Format: x t alpha g K'
    WRITE(output_unit,'(A)') '# Time snapshots: t=0,2.5,5,7.5,10'
    WRITE(output_unit,'(A)') '#'
    
    DO step=1,5
        t_output=save_times(step)
        DO i=1,N
            WRITE(output_unit,'(5(ES24.15E3))') x(i),t_output,&
                alpha_data(i,step),g_data(i,step),K_data(i,step)
        END DO
        !Blank line between time snapshots for easier plotting
        WRITE(output_unit,'(A)') ''
    END DO
    
    CLOSE(output_unit)
    
    PRINT *,'Data written to: evolution_data.dat'
    PRINT *,'Simulation completed'
    
END PROGRAM Main
