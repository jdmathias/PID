%clear all;%close all;
phi11=0.912;phi21=0.03833;I=0.25;C(1)=818.985;K(1)=135;mu(1)=0.039;
gamma=0.3;npas=19;

%% total factor productivity
a(1)=3.79762;
ga(1)=0.079;
for t=2:npas
    ga(t)=ga(1)*exp(-0.006*(t-1)*5);
    a(t)=a(t-1)/(1-ga(t-1));
end
%% baseline carbon intensity
s(1)=0.489;gs(2)=-0.01;
for t=2:npas
    gs(t+1)=gs(t)*(1-0.001)^5;
    s(t)=s(t-1)*exp(gs(t)*5);
end
%% population UN medium trajectory
l=[6838 7349 7758 8142 8501 8839 9157 9454 9725 9969 10184 10376 10548 10701 10837 10954 11055 11142 11213];
l=l/1000;
%%%%%%%%%%%%%%%%%%%%% capital stock K
deltak=0.1;phi1=0;phi2=0.0027;teta2=2.8;deltab=0.025;
%% backstop price
bp(1)=344;
for t=2:npas
    bp(t)=bp(t-1)*(1-deltab);
end
%% teta1 function
teta1=bp.*s/1000/teta2;
Cd(1)=teta1(1)*mu(1).^teta2;
Tat=(C(1)-300*2.13)/213;
Od(1)=phi2*Tat.^2;
Ierror(1)=0;
coe=0.01;
gammaPi=1/0.1930*coe*1;gammaIi=1/1.8233*coe*1;gammaDi=1/0.0034*coe*1;
Vmin(2:npas)=100000000;
Vmax(2:npas)=-1;

for p=1:101
    p
    for i=1:101
        for d=1:101
            gammaP=gammaPi/3+(p-1)/100*(gammaPi*3-gammaPi/3);
            gammaI=gammaIi/3+(i-1)/100*(gammaIi*3-gammaIi/3);
            gammaD=gammaDi/3+(d-1)/100*(gammaDi*3-gammaDi/3);
            for t=2:npas
                error(t)=(C(t-1)-350*2.13)/(350*2.13);
                if t>2
                    Perror(t)=error(t);
                    Ierror(t)=Ierror(t-1)+error(t);
                    Derror(t)=(error(t)-error(t-1));
                else
                    Perror(t)=error(t);
                    Ierror(t)=Ierror(t-1)+error(t);
                    Derror(t)=0;
                end
                mu(t)=mu(t-1)+(gammaP*Perror(t)+gammaI*Ierror(t)+gammaD*Derror(t));
                if mu(t)<0
                    mu(t)=0;
                elseif mu(t)>1
                    mu(t)=1;
                end

                Cd(t)=teta1(t)*mu(t).^teta2;
                Tat=(C(t-1)-300*2.13)/213;
                Od(t)=phi2*Tat.^2;
                K(t)=5*(I*(1-Cd(t))./(1+Od(t))*a(t-1).*K(t-1).^gamma*l(t-1)^(1-gamma))+(1-deltak)^5*K(t-1);
                E(t)=(1-mu(t))*s(t)*a(t)*K(t)^gamma*l(t)^(1-gamma);
                C(t)=E(t)+phi11*C(t-1)+70;
                if C(t)<Vmin(t)
                    Vmin(t)=C(t);
                end
                if C(t)>Vmax(t)
                    Vmax(t)=C(t);
                end
            end

        end
    end
end
figure
x1=2010:5:(2010+(npas-1)*5);x2=x1;Vmin(1)=C(1);Vmax(1)=C(1);
hae=fill([x1 x2(19:-1:1)],[Vmin Vmax(19:-1:1)]/2.13,'r')
set(hae,'facealapha',0.9)