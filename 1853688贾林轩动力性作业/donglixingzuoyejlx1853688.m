clear all;
clc;

%%
%变量预存
n=600:0.1:4000;
m=3880;
g=9.8;
r=0.367;
yita=0.85;
f=0.013;
CDA=2.77;
i0=5.83;
If=0.218;
Iw1=1.798;
Iw2=3.598;
ig=[5.56 2.769 1.644 1 0.793];
L=3.2;
a=1.947;
b=L-a;
hg=0.9;
Ttq = -19.313+295.27*(n/1000)-165.44*(n/1000).^2+40.874*(n/1000).^3-3.8445*(n/1000).^4;

%%
%（1）驱动力与行驶阻力平衡图
ua = zeros(length(ig),length(n));
Ft = zeros(length(ig),length(n));
for i=1:1:length(ig)
    ua(i,:)=0.377*r*n/(ig(i)*i0);  %车速
    Ft(i,:)=Ttq*ig(i)*i0*yita/r;  %驱动力矩
end

Ff=m*g*f;  %行驶阻力
ua0=0:0.1:120;
Fw=CDA/21.15*ua0.^2;  %空气阻力
F=Ff+Fw;

for j = 1:1:length(ig)
    plot(ua(j,:),Ft(j,:))
    hold on
end
plot(ua0,F);
title('驱动力行驶阻力平衡图')
xlabel('ua(km/h)')
ylabel('Ft(N)')
legend('I挡','II挡','III挡','IV挡','V挡','Ff+Fw')

%%
%(2)最高车速，最大爬坡率及附着率

%最高车速
sumF=zeros(length(ig),length(n));
sumF=m*g*f+CDA/21.15*ua.^2; %将阻力曲线与ua调整为相同步长
plot(ua(5,:),sumF(5,:))
k=find(abs(sumF(5,:)-Ft(5,:))<=0.1); %找到两曲线交点
x=ua(5,k(2)); %横坐标
y=sumF(5,k(2)); %纵坐标
plot(x,y,'*') %在驱动力与行驶阻力平衡图中标出该点
legend('I挡','II挡','III挡','IV挡','V挡','Ff+Fw','Ff+Fw','uamax')
uamax=x; %求出最高车速
fprintf('最高车速为%s\n',uamax)

%最大爬坡率
arf=asin((Ft-sumF)/(m*g));
ii=tan(arf); 
figure(2);
for j = 1:1:length(ig)
    plot(ua(j,:),ii(j,:))  %画出爬坡率曲线i-ua
    hold on
end
title('爬坡度图')
xlabel('ua(km/h)')
ylabel('i（%）')

imax=max(ii(1,:));  %找到最大爬坡率
for i=1:1:length(n)
    if ii(1,i)==imax
        ki=i;
    end
end
plot(ua(1,ki),imax,'*')  %标出最大爬坡率所在点
legend('I挡','II挡','III挡','IV挡','V挡','imax')
fprintf('最大爬坡度为%s\n',imax)

%附着率
phi = imax/(a/L+hg/L*imax);
fprintf('附着率为%s\n',phi);

%%
%行驶加速度倒数曲线
Iw=Iw1+Iw2;
drt=zeros(length(ig),length(n));
for i=1:1:length(ig)
    drt(i,:)=1+1/m*Iw/(r^2)+1/m*If*ig(i).^2*i0^2*yita/(r^2);  %旋转质量换算系数
end
du_dt=1./(drt*m).*(Ft-sumF);  %加速度
a1=1./du_dt;  %加速度倒数
figure();
for j = 1:1:length(ig)
    plot(ua(j,:),du_dt(j,:))  %画出加速度曲线
    hold on
end
title('行驶加速度曲线')
xlabel('ua(km/h)')
ylabel('a(m/s^2)')
legend('I挡','II挡','III挡','IV挡','V挡')
figure();
for j = 1:1:length(ig)
    plot(ua(j,:),a1(j,:))  %画出加速度倒数曲线
    hold on
end
axis([0 80 0 6]);
title('行驶加速度倒数曲线')
xlabel('ua(km/h)')
ylabel('1/a(s^2/m)')
legend('I挡','II挡','III挡','IV挡','V挡')

%车速-时间曲线
sumt=zeros(1,length(n));
sumv=zeros(1,length(n));
sumv(1)=ua(2,1);  %在一般动力性分析计算原地起步加速时间时认为在最初时刻汽车已经具有起步档位的最低车速
s=1;ss=0;
for j=2:1:length(ig)-1
    if i==length(n)-1
       s=find(abs(ua(j,:)-max(ua(j-1,:)))<=0.01,1);
    end
    for i=s:1:length(n)-1
    ss=ss+1;
    dua=ua(j,i+1)-ua(j,i);   %表示出Δu
    dt=dua*a1(j,i);  %微元面积表示加速时间
    sumt(ss+1)=sumt(ss)+dt;  %累加微元：积出时间
    sumv(ss+1)=ua(j,i+1);
        if sumv(ss+1)>70
            break
        end
    end
    if sumv(ss+1)>70
            break
    end
end
figure();
for j=1:1:length(ig)
    plot(1/3.6*sumt,sumv);  %画出车速-时间曲线
    hold on;
end
axis([0 40 0 80]);
title('车速-时间曲线')
xlabel('t/s')
ylabel('ua/(km/h)')

