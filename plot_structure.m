function plot_structure(eta_1,eta_2,theta_1,theta_2,theta_3)
% INPUTS: Degrees of freedom as defined in Project
% -------------------------------------------------------------------------
% Geometry data (modify accordingly)
% -------------------------------------------------------------------------
c1 = 0.65;
c2 = 0.5;
h1 = 1.25;
h2 = 1.75;
x1 = 0.22;
x2 = 0.2;
h0 = 1.5;
% -------------------------------------------------------------------------
% Build geometry
% -------------------------------------------------------------------------
% Scale amplitudes of displacements
x = [eta_1,eta_2,theta_1*c1,theta_2*c1,theta_3*c1];
x = 0.5*c1*x/max(abs(x)); % Scale displacements to chordlength
maxangle = max(abs([(x(1)-x(2))/h0,x(3)/c1,x(4)/c1,x(5)/c1]));
if maxangle > 0.4 
    scale = 0.4/maxangle;
else
    scale = 1;
end
% Panels definition
% (1) Streamwise position of hinge point
% (2) Panel shortest distance from vertical wall
% (3) Vertical position of hinge point (in reference configuration)
% (4) Panel chord
% (5) Panel span
% (7) Hinge distance from panel leading edge
% (8) Vertical displacement of hinge point
% (9) Rotation angle of panel from horizontal plane
panels = [   
x1,     h2, 0,  c1, h1, x1, x(1)*scale, x(3)/c1*scale
x1,     0,  0,  c1, h2, x1, x(1)*scale, x(4)/c1*scale
x1+h0,  0,  0,  c2, h2, x2, x(2)*scale, x(5)/c1*scale
];
% Bars definition
% (1) Streamwise position of point 1
% (2) Distance from vertical wall of point 1
% (3) Vertical position of point 1 (in reference configuration)
% (4) Streamwise position of point 2
% (5) Distance from vertical wall of point 2
% (6) Vertical position of point 2 (in reference configuration)
bars = [   
x1,     0,  0,  x1+h0,  0,      0, x(1)*scale, x(2)*scale
x1,     0,  0,  x1,     h1+h2,  0, x(1)*scale, x(1)*scale
x1+h0,  0,  0,  x1+h0,  h2,     0, x(2)*scale, x(2)*scale];
% Determine max vertical displacement
maxz = max(abs(panels(:,3) + panels(:,7)) + abs(panels(:,4).*panels(:,8)));

% -------------------------------------------------------------------------
% Initialize figure
% -------------------------------------------------------------------------
hFig = figure;
uicontrol('Parent',hFig,'Style','pushbutton','String','Animate','Units','pixels','Position',[10 10 60 25],'Visible','on','Callback',@anim);
hold on
% -------------------------------------------------------------------------
% Plot bars
% -------------------------------------------------------------------------
for i = 1:size(bars,1)
    pb(i) = patch(bars(i,[2,5]),bars(i,[1,4]),bars(i,[3,6])+real(bars(i,[7,8])),zeros(1,2),'FaceColor','none','EdgeColor','k','LineWidth',2);
end
% -------------------------------------------------------------------------
% Plot panels in initial configuration (undeformed)
% -------------------------------------------------------------------------
for i = 1:size(panels,1)
    xp = [panels(i,1)-panels(i,6);
          panels(i,1)-panels(i,6);
          panels(i,1)-panels(i,6)+panels(i,4);
          panels(i,1)-panels(i,6)+panels(i,4)];
    yp = [panels(i,2);
          panels(i,2)+panels(i,5);
          panels(i,2)+panels(i,5);
          panels(i,2)];
    zp = [panels(i,3);
          panels(i,3);
          panels(i,3);
          panels(i,3)];
    patch(yp,xp,real(zp),zeros(size(xp)),'FaceColor','none','EdgeColor','k','LineStyle','--');
end
% -------------------------------------------------------------------------
% Plot deformed panels
% -------------------------------------------------------------------------
for i = 1:size(panels,1)
    xp = [panels(i,1)-panels(i,6);
          panels(i,1)-panels(i,6);
          panels(i,1)-panels(i,6)+panels(i,4);
          panels(i,1)-panels(i,6)+panels(i,4)];
    yp = [panels(i,2);
          panels(i,2)+panels(i,5);
          panels(i,2)+panels(i,5);
          panels(i,2)];
    zp = [panels(i,3)+panels(i,7)+panels(i,6)*panels(i,8);
          panels(i,3)+panels(i,7)+panels(i,6)*panels(i,8);
          panels(i,3)+panels(i,7)+(panels(i,6)-panels(i,4))*panels(i,8);
          panels(i,3)+panels(i,7)+(panels(i,6)-panels(i,4))*panels(i,8)];
    pp(i) = patch(yp,xp,real(zp),zeros(size(xp)),'FaceColor',0.5*[1,1,1],'EdgeColor','k');
end
% -------------------------------------------------------------------------
% Set plot properties
% -------------------------------------------------------------------------
view(50,25)
set(gca,'color','none','xcolor','none','ycolor','none','zcolor','none','zlim',[-maxz,maxz]);
axis equal
axis tight
% -------------------------------------------------------------------------
% Animation function
% -------------------------------------------------------------------------
function anim(src,event)
    n = 2; % 2 cycles
    a = linspace(0,n*2*pi,101);
    for t = 1:length(a)
        for k = 1:size(bars,1)
            pb(k).ZData = bars(k,[3,6])+abs(bars(k,[7,8])).*cos(a(t)+angle(bars(k,[7,8])));
        end
        for k = 1:size(panels,1)
            zp = [panels(k,3)+abs(panels(k,7))*cos(a(t)+angle(panels(k,7)))+panels(k,6)*abs(panels(k,8))*cos(a(t)+angle(panels(k,8)));
                  panels(k,3)+abs(panels(k,7))*cos(a(t)+angle(panels(k,7)))+panels(k,6)*abs(panels(k,8))*cos(a(t)+angle(panels(k,8)));
                  panels(k,3)+abs(panels(k,7))*cos(a(t)+angle(panels(k,7)))+(panels(k,6)-panels(k,4))*abs(panels(k,8))*cos(a(t)+angle(panels(k,8)));
                  panels(k,3)+abs(panels(k,7))*cos(a(t)+angle(panels(k,7)))+(panels(k,6)-panels(k,4))*abs(panels(k,8))*cos(a(t)+angle(panels(k,8)))];
            pp(k).ZData = zp;
        end
        zlim([-maxz,maxz]);
        drawnow;
        pause(0.1);
    end
end
end