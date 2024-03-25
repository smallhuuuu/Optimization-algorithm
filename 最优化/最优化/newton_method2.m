function [xk,fk,k]=newton_method2(x0,ess) 
%目标函数f，有确定最优解x*=[0 0],f(x*)=400.0
syms x1 x2; 
f=100*(2*x1^2+2)^2+x2^2;
%构造目标函数f的梯度函数 
fx=diff(f,x1); 
fy=diff(f,x2);
gf=[fx fy]';
%求Hesse阵
fxx=diff(fx,x1); 
fxy=diff(fx,x2);
fyx=diff(fy,x1); 
fyy=diff(fy,x2);
Hess=[fxx,fxy;
      fyx,fyy];
%初始点的梯度和函数值,,赋值
xk=x0'; 
fk=subs(f,[x1 x2],x0);
gk=subs(gf,[x1 x2],x0);
Hk=subs(Hess,[x1 x2],x0);
k=0;  
%进入迭代循环
while((norm(gk)>ess)&&(k<10))%迭代终止条件--满足精度条件
%确定搜索方向dk，步长为1
        dk=-Hk\gk;   
        xk=xk+dk;     
%计算下一点的函数值和梯度
        fk=subs(f,[x1 x2],xk');
        gk=subs(gf,[x1 x2],xk');
        Hk=subs(Hess,[x1 x2],xk');
%迭代次数自增
        k=k+1
end


