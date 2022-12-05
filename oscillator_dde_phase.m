clear
%==================limit cycle=================
lags = [1 2];
tspan = linspace(0,80,100);
y0 = [600,1000];
sol = dde23(@pn,lags,@(t)yhist(t,y0),tspan);

t=sol.x;
y=sol.y;
hold on
p1 = plot(y(1,60:end),y(2,60:end));
p1.LineWidth = 1;
p1.Color = '#17202A';%'#3498DB';
box on

%=================nullcline=====================
hold on
max_x = 1100;
max_y = 1200;
x_lim = [0,max_x];
y_lim = [0,max_y];
[H_null,S_null]=nullcline(x_lim,y_lim);

p2 = plot(H_null(:,1),H_null(:,2));
p2.LineWidth = 2;
p2.Color = '#2874A6';
p2.LineStyle = '-';
p3 = plot(S_null(:,1),S_null(:,2));
p3.LineWidth = 2;
p3.Color = '#CB4335';
p3.LineStyle = '-';


%===================quiver======================

%%


xList=linspace(1,max_x,11);
yList=linspace(1,max_y,11);
U=[];
V=[];
[X, Y] = meshgrid(xList, yList);
ts = 1:2;
LineLength=0.06;
for i=1:length(X)
    for j=1:length(Y)
        y00=[X(i,j),Y(i,j)];
        sol = dde23(@pn,lags,@(t)yhist(t,y00),ts);
        t=sol.x;
        y=sol.y;
        U(i,j)=y(1,end)-y(1,1);
        V(i,j)=y(2,end)-y(2,1);
        ah = annotation('arrow','headStyle','cback2','Color','#99A3A4','HeadLength',5,'HeadWidth',2.5);
        set(ah,'parent',gca);
        ah.Position=[X(i,j),Y(i,j),U(i,j)*LineLength,V(i,j)*LineLength];
    end
end
% quiver(X, Y, U, V,0.5,'Color','#99A3A4');
ylim([0,max_y]);
xlim([0,max_x]);

%% =================time traj=====================
axes('Position',[.68 .7 .2 .2])
tspan = linspace(0,80,100);
y0 = [600,1000];
sol = dde23(@pn,lags,@(t)yhist(t,y0),tspan);

t=sol.x;
y=sol.y;
p4 = plot(y(2,1:end));
p4.LineWidth = 1;
p4.Color = '#EC7063';

function dydt = pn(t,y,Z)
    ylag1 = Z(:,1);
    ylag2 = Z(:,2);
    
    as = 30.5;
    ah = 183;
    ah0 = 0.1;
    as0 = 0.1;
    beta = 4.6;
    dm = 0.3;
    dh = 3.8;
    ds = 0.2;
    Kh = 326;
    Ks = 185;
%     Khs=30;
    n1 = 3;
    n2 = 4.8;
    K = beta/dm;
    dydt = [K*ah0 + K*ah*Ks^n2/(Ks^n2+ylag2(2)^n2)-dh*y(1); %dHdt
            K*as0 + K*as*ylag1(1)^n1/(Kh^n1+ylag1(1)^n1)-ds*y(2);  %dSdt
        ];
%         dydt = [K*a0 + K*ah*Ks^n2/(Ks^n2+ylag2(2)^n2)*ylag1(1)^n1/(Khs^n1+ylag1(1)^n1)-dh*y(1); %dHdt
%             K*as*ylag1(1)^n1/(Kh^n1+ylag1(1)^n1)-ds*y(2);  %dSdt
%         ];
    
end


function s=yhist(t,y0)
    s = y0; %sir2,HAP
end

function [H_null,S_null]=nullcline(xlim,ylim)
    steps = 2000;
    x_list = linspace(xlim(1),xlim(2),steps+1);
    y_list = linspace(ylim(1),ylim(2),steps+1);
    m=1;
    n=1;
    for i=2:length(x_list)
        xitem_1 = x_list(i-1);
        xitem_2 = x_list(i);
        
        for j=2:length(y_list)
            yitem_1 = y_list(j-1);
            yitem_2 = y_list(j);
            
            [x1,y1] = quiver_(xitem_1,yitem_1);
            [x2,y2] = quiver_(xitem_2,yitem_1);
            if sign(x1)~=sign(x2)
                H_null(m,1) = xitem_1; 
                H_null(m,2) = yitem_1; 
                m=m+1;
            elseif sign(y1)~=sign(y2)
                S_null(n,1) = xitem_1; 
                S_null(n,2) = yitem_1;
                n=n+1;
            end
        end
    end
            
end

function [u,v]=quiver_(H,S)
    as = 30.5;
    ah = 183;
    ah0 = 0.1;
    as0 = 0.1;
    beta = 4.6;
    dm = 0.3;
    dh = 3.8;
    ds = 0.2;
    Kh = 326;
    Ks = 185;
%     Khs=30;
    n1 = 3;
    n2 = 4.8;
    K = beta/dm;
    u = K*ah0 + K*ah*Ks.^n2./(Ks^n2+S.^n2)-dh*H;
    v =  K*as0+K*as*H.^n1./(Kh^n1+H.^n1)-ds*S;
end


