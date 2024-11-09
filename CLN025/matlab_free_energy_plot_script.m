close all
clear all
hold off
NGRID=130;
xx=load('../best3_scaled_rgyr_pc1_pc2_1429pt.txt');
ff=load('./best3_free_energy_1429pt.txt');
x=xx(:,1);
y=xx(:,2);
z=xx(:,3);
[minx,maxx] = bounds(x);
[miny,maxy] = bounds(y);
[minz,maxz] = bounds(z);
[xq,yq,zq] = meshgrid(linspace(minx,maxx,NGRID),linspace(miny,maxy,NGRID),linspace(minz,maxz,NGRID));
f=ff(:,1)*2.479;
fq = griddata(x,y,z,f,xq,yq,zq,'natural');
patch(isosurface(xq,zq,yq,fq, 9.),'FaceColor','red','FaceAlpha',0.1,'EdgeColor','none');
patch(isosurface(xq,zq,yq,fq, 4.),'FaceColor','blue','FaceAlpha',0.3,'EdgeColor','none');
patch(isosurface(xq,zq,yq,fq, -0.25),'FaceColor','blue','FaceAlpha',1.,'EdgeColor','none');
camlight;
campos([-15.2,-17.8,10.8]);
view(-52.3,22.7);
camtarget([7.8,-0.049,0.104]);
camup([0 0 1]);
camva(9.66);
light("Style","local","Position",[-17,-29,15]);
xlabel('RGYR'); ylabel('PC2'); zlabel('PC1');
grid on;
set(gca,'FontName','Helvetica','FontSize',12);
lgd = legend(' 9.0 kJ/mol',' 4.0 kJ/mol','-0.25 kJ/mol');
title(lgd,'Free energy')
lighting gouraud
