%close all;
clear all
for pcoefB=1:101
    pcoefB
    coefB=0.5+(pcoefB-1)/100*(3-0.5);    
    gamma=0.1;beta=0.3*coefB;R0=beta/gamma;alpha=0.1;temps=200;dt=1;
    S(1)=0.999;I(1)=0.001;R(1)=0;u(1)=0;
    popsize=67000000;ICUobj=3000;
    error(1)=(0.02/coefB*I(1)*popsize-ICUobj)/ICUobj;
    Perror=error(1);Ierror=error(1)*dt;
    npid=500;
    gammaP=4;gammaI=2;
    t=1;critere=0;
    while critere==0
        t=t+1;
        R0=beta/gamma*(1-u(t-1));
        S(t)=S(t-1)-(R0*gamma*I(t-1)*S(t-1)-alpha*R(t-1))*dt;
        I(t)=I(t-1)+(R0*S(t-1)-1)*gamma*I(t-1)*dt;
        R(t)=R(t-1)+(gamma*I(t-1)-alpha*R(t-1))*dt;
        if S(t)<0
            S(t)=0;
        end    
        if I(t)<0
            I(t)=0;
        end 
        if R(t)<0
            I(t)=0;
        end 
        error(t)=(0.02/coefB*I(t)*popsize-ICUobj)/ICUobj;
        Perror=error(t);
        Ierror(t)=Ierror(t-1)+error(t)*dt;
        u(t)=(gammaP*Perror+gammaI*Ierror(t)); 
        if u(t)<0
            u(t)=0;
        elseif u(t)>1
            u(t)=1;
        end
        if t>20
            if sum(abs(error(t-20:t)))<0.01 
              opt(pcoefB)=t;
              critere=1;
            end
            if t==temps             
               critere=1;
               opt(pcoefB)=NaN;
            end
        end
    end

end
%figure
hold on
plot([0.5:(1/100*(3-0.5)):3],opt,'g')
