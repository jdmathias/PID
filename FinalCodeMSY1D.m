close all;clear all
%initialisation
r=0.25; K=4; alpha=0.25; temps=2000;dt=1;% Model parameters
coef=1;%coef for increasing/decreasing MSY target. coef=1.4 for reproducing figure d
xMSY=coef*(K+alpha)/2;
yMSY=r*(K-xMSY).*(xMSY-alpha)./xMSY;
MSY=xMSY*yMSY;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaP=0.05;gammaI=0.02;% gain gammaI=0.02 for PI figure
ErrorType=3;%1: biomass error; 2: effort error; 3: harvest error
npasX=200;npasY=200;
pasX=4/npasX;pasY=1/npasY;
[X,Y] = meshgrid(pasX/2:pasX:(4-pasX/2), pasY/2:pasY:(1-pasY/2));
%%%%%%%%% loop for each point of the grid
for ix=1:npasX % loop on biomass values
    ix
    for jy=1:npasY % loop on effort values
        x(1)=X(jy,ix);y(1)=Y(jy,ix);
        t=1;critere=0;Ierror(1)=0;
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
                Ierror(t)=Ierror(t-1)+error(t);
            elseif ErrorType==2
                error(t)=-(y(t-1)-yMSY)/yMSY;%effort error
                Ierror(t)=Ierror(t-1)+error(t);
            elseif ErrorType==3
                error(t)=-(x(t-1)*y(t-1)-MSY)/MSY;%harvest error
                Ierror(t)=Ierror(t-1)+error(t);
            end
            y(t)=y(t-1)+gammaP*error(t)+gammaI*Ierror(t);
            if y(t)>1
                y(t)=1;
            elseif y(t)<0
                y(t)=0;
            end
            if t==temps
                critere=1;
                to(ix,jy)=NaN;
            elseif t>4 && sum(abs(error(t-4:t)))<0.0001 && abs(x(t)-xMSY)<0.1 
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
contourf(Y(:,1),X(1,:),to)
plot(yeq,xeq,'r',ymsye,xeq,'g','Linewidth',3)
plot(yeq,xeq,'r','Linewidth',3)
plot(yMSY,xMSY,'p','Linewidth',9,'color','red')
ylim([0;4]);xlim([0;1])