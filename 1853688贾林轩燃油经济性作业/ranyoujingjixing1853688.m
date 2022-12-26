%% info
% 1853688 贾林轩
% 汽车燃油经济性作业
%%
clear all
clc

generate_fig_one();
[uaIVnV,b] = generate_fig_two();
calculate(uaIVnV,b)
%% const 
function const = generate_const()
    const.n = 600:1:4000; %发动机转速
    const.m = 3880;%总质量
    const.g = 9.8;%重力加速度
    const.r = 0.367;%车轮半径
    const.nT = 0.85;%传动系机械效率
    const.f = 0.013;%滚动阻力系数
    const.CDA = 2.77;%空气阻力系数*迎风面积
    const.i0 = 5.83;%主减速器传动比
    const.If = 0.218;%飞轮转动惯量
    const.Iw1 = 1.798;%两前轮转动惯量
    const.Iw2 = 3.598;%四后轮转动惯量
    const.Iw = const.Iw1+const.Iw2;
    const.ig = [5.56,2.769,1.644,1.00,0.793];%变速器传动比
    const.L = 3.2;%轴距
    const.a = 1.947;%质心至前轴的距离
    const.hg = 0.9;%质心高
    const.rou = 0.74;%取汽油密度0.74g/ml
    const.Qid = 0.299;%怠速耗油量
    const.G = const.m*const.g;
    const.delta = 1+1/const.m*const.Iw/const.r^2+1/const.m*(const.If*const.ig.^2*const.i0^2*const.nT)/(const.r^2);%车辆旋转质量换算系数
    const.B = [815,1326.8,-416.46,72.379,-5.8629,0.17768;
        1207,1354.7,-303.98,36.657,-2.0553,0.043072;
        1614,1284.4,-189.75,14.524,-0.51184,0.0068164;
        2012,1122.9,-121.59,7.0035,-0.18517,0.0018555;
        2603,1141.0,-98.893,4.4763,-0.091077,0.00068906;
        3006,1051.2,-73.714,2.8593,-0.05138,0.00035032;
        3403,1233.9,-84.478,2.9788,-0.047449,0.00028230;
        3804,1129.7,-45.291,0.71113,-0.00075215,-0.000038568];%拟合系数
end



%% 汽车功率平衡图绘制
function generate_fig_one()
    config = generate_const();
    config.Ttq = -19.313+295.27*(config.n/1000)-165.44*(config.n/1000).^2+40.874*(config.n/1000).^3-3.8445*(config.n/1000).^4;%发动机转矩
    ig_length = length(config.ig);
    n_length = length(config.n);
    ua = zeros(ig_length,n_length);
    Ft = zeros(ig_length,n_length);
    for i = 1:ig_length
        ua(i,:)=0.377*config.r*config.n/(config.ig(i)*config.i0);
        Ft(i,:) = config.Ttq.*config.ig(i)*config.i0*config.nT/config.r;
    end%各档位时速及驱动力
    Pe = Ft.*ua/(config.nT*3600);%发动机输出功率
    figure()
    title('汽车功率平衡图')
    for i = 1:ig_length
        plot(ua(i,:),Pe(i,:),'-')
        hold on
    end
    ua1 = 0:0.1:120;
    Ffw1 = config.CDA*ua1.^2/21.15+config.G*config.f;%行驶阻力
    Pfw1 = Ffw1.*ua1/(3600*config.nT);%行驶阻力在发动机输出端消耗的功率
    plot(ua1,Pfw1)
    grid on
    legend('I挡','II挡','III挡','IV挡','V挡','(Pf+Pw)/nT')
    xlabel('ua/(km/h)')
    ylabel('P/(kW)')

end

%% 最高档和次高档百公里油耗
function [uaIVnV,b] = generate_fig_two()
    config = generate_const();
    [c,~] = size(config.B);
    uaIVnV = zeros(c,2);
    PeIVnV = zeros(c,2);
    b = zeros(c,2);
    uaT = (25:10:85)';
    PeT = zeros(length(uaT),2);
    for i = [length(config.ig)-1,length(config.ig)]
        uaIVnV(:,i-length(config.ig)+2)=0.377*config.r*config.B(:,1)/(config.ig(i)*config.i0);
        PeIVnV(:,i-length(config.ig)+2) = (config.CDA*uaIVnV(:,i-length(config.ig)+2).^2/21.15+config.G*config.f).*uaIVnV(:,i-length(config.ig)+2)/(config.nT*3600);
        b(:,i-length(config.ig)+2) = config.B(:,2)+config.B(:,3).*PeIVnV(:,i-length(config.ig)+2)+config.B(:,4).*PeIVnV(:,i-length(config.ig)+2).^2+config.B(:,5).*PeIVnV(:,i-length(config.ig)+2).^3+config.B(:,6).*PeIVnV(:,i-length(config.ig)+2).^4;
        PeT(:,i-length(config.ig)+2) = (config.CDA*uaT.^2/21.15+config.G*config.f).*uaT/(config.nT*3600);
    end%各转速对应时速及阻力消耗功率

    bT = zeros(length(uaT),2);
    for i = 1:length(uaT)
        for j = [1,2]
            index = find(uaIVnV(:,j)-uaT(i)<0, 1, 'last' );
            bT(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(uaT(i)-uaIVnV(index,j))+b(index,j);
        end
    end%插值得到系列车速下的次高档及最高档的燃油消耗率
    Qs = zeros(length(uaT),2);
    for i = [1,2]
        Qs(:,i) = PeT(:,i).*bT(:,i)./(1.02*uaT*config.rou*config.g);
    end%单位转换
    figure()
    title('次高档和最高档百公里油耗')
    plot(uaT,Qs(:,1),'-s')
    hold on
    plot(uaT,Qs(:,2),'->')
    legend('次高档','最高档')
    grid on
    xlabel('ua/(km/h)')
    ylabel('Qs/(L/100km)')
    
end

%% 按JB3352-83规定的六工况循环行驶计算百公里油耗
function calculate(uaIVnV,b)
    config = generate_const();
    ua_acc_1 = 25:0.001:40;%第一段加速
    Pe_acc_1 = zeros(2,length(ua_acc_1));
    acc_1 = 0.25;%第一段加速度
    for i = [length(config.ig)-1,length(config.ig)]
        Pe_acc_1(i-length(config.ig)+2,:) = (config.CDA*ua_acc_1.^2/21.15+config.G*config.f+config.delta(i)*config.m*acc_1).*ua_acc_1/(config.nT*3600);
    end%第一段加速过程中的功率
    b_acc_1 = zeros(length(ua_acc_1),2);%第一段加速过程中的耗油率
    for i = 1:length(ua_acc_1)
        for j = [1,2]
            index = find(uaIVnV(:,j)-ua_acc_1(i)<0, 1, 'last' );
            b_acc_1(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(ua_acc_1(i)-uaIVnV(index,j))+b(index,j);
        end
    end%插值得到系列车速下的次高档及最高档的燃油消耗率
    Qs = Pe_acc_1.*b_acc_1'./(367.1*config.rou*config.g);
    Q_acc_1 = sum(Qs,2)*(0.001)/(3.6*acc_1);%第一段加速时IV、V挡的耗油量

    ua_acc_2 = 40:0.001:50;%第二段加速
    Pe_acc_2 = zeros(2,length(ua_acc_2));
    acc_2 = 0.2;%第二段加速度
    for i = [length(config.ig)-1,length(config.ig)]
        Pe_acc_2(i-length(config.ig)+2,:) = (config.CDA*ua_acc_2.^2/21.15+config.G*config.f+config.delta(i)*config.m*acc_1).*ua_acc_2/(config.nT*3600);
    end%第二段加速过程中的功率
    b_acc_2 = zeros(length(ua_acc_2),2);%第二段加速过程中的耗油率
    for i = 1:length(ua_acc_2)
        for j = [1,2]
            index = find(uaIVnV(:,j)-ua_acc_2(i)<0, 1, 'last' );
            b_acc_2(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(ua_acc_2(i)-uaIVnV(index,j))+b(index,j);
        end
    end%插值得到系列车速下的次高档及最高档的燃油消耗率
    Qs = Pe_acc_2.*b_acc_2'./(367.1*config.rou*config.g);
    Q_acc_2 = sum(Qs,2)*(0.001)/(3.6*acc_2);%第二段加速时IV、V挡的耗油量

    dua_acc_3 = 25;%第三段减速
    acc_3 = 0.36;%第三段加速度
    Q_acc_3 = config.Qid*dua_acc_3/(3.6*0.36);%减速时IV、V挡的耗油量

    ua_e = [25,40,50];%匀速
    b_e = zeros(length(ua_e),2);%匀速过程中的耗油率
    for i = 1:length(ua_e)
        for j = [1,2]
            index = find(uaIVnV(:,j)-ua_e(i)<0, 1, 'last' );
            b_e(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(ua_e(i)-uaIVnV(index,j))+b(index,j);
        end
    end%插值得到系列车速下的次高档及最高档的燃油消耗率
    Pe_e = (config.CDA*ua_e.^2/21.15+config.G*config.f).*ua_e/(config.nT*3600);%匀速过程中的功率
    t_e = [7.2,22.5,18]';%匀速行驶的时间
    Q_e = sum(Pe_e'.*b_e.*t_e/(367.1*config.rou*config.g))';%匀速过程中的耗油量

    Q = Q_acc_1+Q_acc_2+Q_acc_3+Q_e;%六工况总油耗

    Qs = 100*Q/1075;
    fprintf('IV档的百公里油耗为 %f，V档的百公里油耗为 %f',Qs(1),Qs(2));

end







