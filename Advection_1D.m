% u_t+c*u_x=0
clear all;
close all;
nraf=5;
order=1;
for n=1:nraf
    nx=2^(n+6);
    aa=-20;bb=20;
    dx=(bb-aa)/nx;
    x=aa+dx:dx:bb;
    c=.5;x0=-5;
    uex=@(t,x) exp(-(x-x0-c*t).^2/4.);
    CFL=.4;
    Tf=5.;
    dt=CFL*dx/c;
    u0=uex(0,x);
    naff=1;
    it=0;
    for t=dt:dt:Tf
%         if mod(it,naff)==0
%             plot(x,u0,'x-b');%,x,uex(t,x),'x-r');
%             drawnow;
%             pause(.1)
%         end
        if order==1
            for i=2:nx
                u(i)=u0(i)-c*dt*(u0(i)-u0(i-1))/dx;
            end
            u(1)=u(nx);
        end
        if order==2
            for i=3:nx
                u(i)=u0(i)-c*dt*(3*u0(i)-4*u0(i-1)+u0(i-2))/(2*dx)+c^2*dt^2*(u0(i)-2*u0(i-1)+u0(i-2))/(2*dx^2);
            end
            u(1)=u(nx-1);u(2)=u(nx);
        end
        u0=u;
        it=it+1;
    end
    err_two(n)=sqrt(dx*sum((u-uex(t,x)).^2));
    Dt(n)=dt;
    Dx(n)=dx;
end
ratedt=log(err_two(1:nraf-1)./err_two(2:nraf))./log(Dt(1:nraf-1)./Dt(2:nraf))
ratedx=log(err_two(1:nraf-1)./err_two(2:nraf))./log(Dx(1:nraf-1)./Dx(2:nraf))