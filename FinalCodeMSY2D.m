close all;clear all
%initialisation
r=0.25; K=4; alpha=0.25; temps=2000;dt=1;% Model parameters
beta=0.05;price=4.5;cost=1.5;uref=4.5;
xMSY=(K+alpha)/2;yMSY=r*(K-xMSY).*(xMSY-alpha)./xMSY;MSY=xMSY*yMSY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaP=0.2;gammaI=0.2;gammaD=0.2% Values are good for the biomass error
ErrorType=1;%1: biomass error; 2: effort error; 3: harvest error
npasX=200;npasY=200;
pasX=4/npasX;pasY=1/npasY;
[X,Y] = meshgrid(pasX/2:pasX:(4-pasX/2), pasY/2:pasY:(1-pasY/2));
%%%%%%%%% loop for each point of the grid
for ix=1:npasX % loop on biomass values
    ix
    for jy=1:npasY % loop on effort values
        x(1)=X(jy,ix);y(1)=Y(jy,ix);
        t=1;critere=0;Ierror(t)=0;u(1)=uref;
        while critere==0 %%%%%% time loop
            t=t+1;
            Xp=x(t-1);Yp=y(t-1);
            for tat=1:10
                Xp=Xp+(-Yp*Xp+r*(Xp-alpha)*(K-Xp))*dt/10;%small time step for taking into account non linearities of the dynamics
            end
            x(t)=Xp;
            if x(t)<0
                x(t)=0;
            end    
            if ErrorType==1
                error(t)=(x(t-1)-xMSY)/xMSY;%biomass error
                if t>5
                    Ierror(t)=sum(error(t-5:t));
                else
                   Ierror(t)=sum(error(1:t));
                end
                if t>2
                    Derror(t)=error(t)-error(t-1);
                else
                    Derror(t)=0;
                end
                
            elseif ErrorType==2
                error(t)=-(y(t-1)-yMSY)/yMSY;%effort error
                if t>5
                   Ierror(t)=sum(error(t-5:t));
                else
                   Ierror(t)=sum(error(1:t));
                end
                if t>2
                    Derror(t)=error(t)-error(t-1);
                else
                    Derror(t)=0;
                end
            elseif ErrorType==3
                error(t)=-(x(t-1)*y(t-1)-MSY)/MSY;%harvest error
                if t>5
                   Ierror(t)=sum(error(t-5:t));
                else
                   Ierror(t)=sum(error(1:t));
                end
                if t>2
                    Derror(t)=error(t)-error(t-1);
                else
                    Derror(t)=0;
                end
            end
            u(t)=u(t-1)-(gammaP*error(t)+gammaI*Ierror(t)+gammaD*Derror(t));
            if u(t)<0
                u(t)=0;
            end
            if u(t)>uref*2
                u(t)=uref*2;
            end
            y(t)=y(t-1)+beta*y(t-1)*(1-y(t-1))*(price*x(t-1)-cost-u(t-1));
            if y(t)>1
                y(t)=1;
            elseif y(t)<0
                y(t)=0;
            end
            if t==temps
                critere=1;
                to(ix,jy)=NaN;
            elseif t>4 && abs(x(t)-xMSY)<0.1 && sum(abs(error(t-4:t)))<0.0001 
                critere=1;
                to(ix,jy)=t;
            end
        end
    end
end

%%%%%%%%%%%%%% Figure
figure
hold on
xeq=0:0.01:5;
yeq=r*(K-xeq).*(xeq-alpha)./xeq;
ymsye=MSY./xeq;
contourf(Y(:,1),X(1,:),to,100)
plot(yeq,xeq,'r','Linewidth',3)
plot(yMSY,xMSY,'p','Linewidth',9,'color','red')
ylim([0;4]);xlim([0;1])