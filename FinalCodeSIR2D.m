close all;clear all
variant=1;%1: low contagious/more dangerous variant;2: high contagious/less dangerous variant;
if variant==1
    coefB=1;
elseif variant==2
    coefB=2;
end
gamma=0.1;beta=0.3*coefB;R0=beta/gamma;alpha=0.1;temps=100;dt=1;
S(1)=0.999;I(1)=0.001;R(1)=0;u(1)=0;
popsize=67000000;ICUobj=3000;
error(1)=(0.02/coefB*I(1)*popsize-ICUobj)/ICUobj;
Perror=error(1);Ierror=error(1)*dt;
npid=500;
for p=1:npid
    p
    for i=1:npid
        gammaP=1+(p-1)*5/npid;gammaI=(i-1)*2/npid;
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
                  opt(p,i)=t;
                  critere=1;
                end
                if t==temps             
                   critere=1;
                   opt(p,i)=NaN;
                end
            end
        end
    end
end

%%%%%%%%%%% Figures
figure
[v,loc] = min(opt(:));
[ii,jj,k] = ind2sub(size(opt),loc);
gammaP=1+(ii-1)*5/npid;gammaI=(jj-1)*2/npid;
contourf(([1:npid]-1)*2/npid,1+([1:npid]-1)*5/npid,opt);