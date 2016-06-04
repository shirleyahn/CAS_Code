%x=linspace(-1.5,1.5);
x=linspace(-1.0,1.0);
%y=linspace(-0.5,1.25);
y=linspace(-1.0,1.0);
[X,Y]=meshgrid(x,y);
Z=exp(-X.^2)+Y.^2;
% three hole potential
%Z=4*(X.^2+Y.^2-1).^2.*Y.^2-exp(-4*((X-1).^2+Y.^2))-exp(-4*((X+1).^2+Y.^2))+exp(8*(X-1.5))+exp(-8*(X+1.5))+exp(-4*(Y+0.25))+0.2*exp(-8*X.^2);
%Z=2.25*(X-0.5*sin(2*pi*Y)).^2+2.25*cos(2*pi*Y);
%Z=3.0*(sqrt(X.^2+Y.^2)-3.0).^2+2.25*cos(2*tan(Y./X))-4.5*cos(2*tan(Y./X));
contour3(X,Y,Z,50),colorbar
xlabel('x')
ylabel('y')
%%
Z=zeros(length(x),length(y));
x_index=0;
for m=-1.5:0.01:1.5
    x_index=x_index+1;
    y_index=0;
    for n=-0.5:0.01:1.25
        y_index=y_index+1;
        Z(x_index,y_index)=4*(m^2+n^2-1)^2*n^2-exp(-4*((m-1)^2+n^2))-exp(-4*((m+1)^2+n^2))+exp(8*(m-1.5))+exp(-8*(m+1.5))+exp(-4*(n+0.25))+0.2*exp(-8*m^2);
    end
end
