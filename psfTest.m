lambda = 0.193;
xmin = -5;
xmax = 5;
res = 512;
NA = 0.01;
Zern = [0, 0, 0, 0.1, 0.2, -0.2,0.2, 0.2, 0.5,0.6];
PSF=Zernike2PSF(Zern,lambda, NA, xmin, xmax, res);hold on
axis equal tight;
xlabel('um'),ylabel('um');
set(gca,'FontSize',14);
% axis off; axis on
th= linspace(0,2*pi,100);
x=3.65*cos(th);
y=3.65*sin(th);
h=plot(x,y,'r');
set(h,'LineWidth',2);