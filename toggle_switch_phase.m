close all
clear

%% set parameters
k_hap_prod = 80;
k_hap_deg = 0.08;
k_hap_prod_base = 3; %hap4 leftmost line
KM_SonH = 328;
KM_pos = 450; %hap4 mid line

Total = 1200;
k_sir_prod = 42;
k_sir_deg = 0.015;
KM_HonS =400;
KM_pos2 = 4180; %sir2 u shape
k_sir_prod_base = 0.7; %sir2 bottom line

n1 = 4;
n3 = 4;
n4 = 2;

neg_strength1 = 0.35;
neg_strength2 = 0.95;

%%
%% get the vector for quiver plot, and save data for later plotting with python 
figure
h_2 =linspace(0,1100,10);
s_2 =linspace(0,Total,10);
[X, Y] = meshgrid(h_2,s_2);

neg_fac = (neg_strength1*X.^n4+KM_HonS^n4)./(X.^n4+KM_HonS^n4);
ds = k_sir_prod*neg_fac.*(Y.^n3./(Y.^n3+KM_pos2^n3)).*(Total-Y) + k_sir_prod_base - k_sir_deg*Y;

neg_fac2 = (neg_strength2*Y.^n4+KM_SonH^n4)./(Y.^n4+KM_SonH^n4);
dh = k_hap_prod*neg_fac2.*(X.^n3./(X.^n3+KM_pos^n3)) + k_hap_prod_base - k_hap_deg*X;
% 
%  ylim([0 1210])
%  xlim([0,1110])

box on
hq = quiver(X,Y,dh,ds,0.9,'Color','#99A3A4');
hq.AutoScaleFactor=0.5;
hq.LineWidth=1;
xlabel('hap/hem');
ylabel('sir2')


LineLength=4;
% 
% for i=1:size(X,1)
%     for j=1:size(X,2)
%         ah = annotation('arrow','headStyle','cback2','Color','#99A3A4','HeadLength',5,'HeadWidth',2.5);
%         set(ah,'parent',gca);
%         ah.Position=[X(i,j),Y(i,j),dh(i,j)*LineLength,ds(i,j)*LineLength];
%         
%     end
% end
hold on
%% find the sir2 roots at steady state
for ii = 0:110
    h(ii+1) = ii*10;
    neg_fac = (neg_strength1*h(ii+1)^n4+KM_HonS^n4)/(h(ii+1)^n4+KM_HonS^n4);
    s_roots(ii+1,:) = roots([(-k_sir_prod*neg_fac - k_sir_deg) (k_sir_prod*neg_fac*Total+k_sir_prod_base) 0 0 -k_sir_deg*KM_pos2^4  k_sir_prod_base*KM_pos2^4]);

end
%% get rid of imaginary roots
image_roots = abs(imag(s_roots))>3;
s_roots2 = s_roots;
s_roots2(image_roots) = nan;
s_line = real(s_roots2);
h_line = h;
p1=plot(h_line, s_line,  'Color',' #CB4335', 'LineWidth', 2);
% yticks([0,200,400,600,800,1000,1200])
 ylim([0 Total+100])
%  xlim([0,1100])
hold on

%% find the hap roots at steady state
for ii = 0:Total
    s(ii+1) = ii;
    neg_fac2 = (neg_strength2*s(ii+1)^n4+KM_SonH^n4)/(s(ii+1)^n4+KM_SonH^n4);
    h_roots(ii+1,:) = roots([(-k_hap_deg) (k_hap_prod_base+k_hap_prod*neg_fac2) 0 0 -k_hap_deg*KM_pos^4  k_hap_prod_base*KM_pos^4]);
end






%% get rid of imaginary roots
image_roots = abs(imag(h_roots))>80;
h_roots2 = h_roots;
h_roots2(image_roots) = nan;
h_line2 = real(h_roots2);
s_line2 = s;
plot(h_line2,s_line2, 'Color',' #2874A6', 'LineWidth', 2)

%%

% save('./data3.mat', 'X','Y', 'dh', 'ds', 's_line', 'h_line', 's_line2', 'h_line2');
% print(gcf, '-dpdf', '-r1000', 'toggle_phase')
