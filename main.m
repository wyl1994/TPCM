
%%%%%%%%%%%%%%%%%%提出的TPCM算法的部分结果
%%%%%%%%%%%%%%%%大论文图4.9的TPCM示例结果

 %%%%%%%%%%%%%%%%%%%%%%%%%%%单向耦合的情况
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%用logist模型验证收敛性，时间长度的鲁棒性,单向耦合情况

  clear
a = 0.2;
b = 0.8;

lag=0;  tau=1;   dim=5;   
rx=3.78; ry=3.77; 

ll=[200:200:4000];
 for kk=1:length(ll)
     clear X
     clear Y
     clear Z
 %X -> Y
maxL=ll(kk);
for k=1:100
X(1)=(b-a).*rand(1,1) + a; Y(1)=(b-a).*rand(1,1) + a;    Z(1)=(b-a).*rand(1,1) + a; 

tau1=1;  

for j=2:maxL
    X(j)=X(j-1)*(rx-rx*X(j-1)-0.0*Y(j-1));
    Y(j)=Y(j-1)*(ry-ry*Y(j-1)-0.08*X(j-1));
    Z(j)=Z(j-1)*(ry-ry*Z(j-1)-0.0*X(j-1));
    
end

X=X(10:end);
Y=Y(10:end);   %删除前10个序列
Z=Z(10:end);   %删除前10个序列
% [C]=CCS_egi(X',Y',0,tau,dim);                  %%%%%%%%%%%%%%%我们的ccs方法
% 
% cxy2(kk,k)=C(1);
% cyx2(kk,k)=C(2);


[C]=TPCM(X',Y',0,tau,dim);

cxy3(kk,k)=C(1);
cyx3(kk,k)=C(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%随机洗牌之后的结果

   p1=randperm(length(X)); X=X((p1(1:end)));   %%%%%%%%%%%吧原始的xx 这个打散顺序
   p1=randperm(length(Y)); Y=Y((p1(1:end)));   %%%%%%%%%%%吧原始的xx 这个打散顺序

% [C]=CCS_egi(X',Y',0,tau,dim);
% 
% cxy2_su(kk,k)=C(1);
% cyx2_su(kk,k)=C(2);


[C]=TPCM(X',Y',0,tau,dim);

cxy3_su(kk,k)=C(1);
cyx3_su(kk,k)=C(2);



        
end
 end


 
 
 %%%%%%%%%%%%%%%画图
 
 figure
 errorbar(ll,mean(cxy3'),std(cxy3'),'LineWidth',2,'marker','diamond')  
 hold on
 errorbar(ll,mean(cyx3'),std(cyx3'),'LineWidth',2,'marker','diamond')  
 hold on
  errorbar(ll,mean(cxy3_su'),std(cxy3_su'),'LineWidth',2,'marker','diamond')  
 hold on
 errorbar(ll,mean(cyx3_su'),std(cyx3_su'),'LineWidth',2,'marker','diamond')  
   set(gca,'FontSize',20);
  xlabel('Data length') 
  ylabel('TPCM') 
 legend('X->Y','Y->X','X(shuffle)->Y','Y(shuffle)->X')
 
 legend('X->Y','X(shuffle)->Y','Y->Z','Y(shuffle)->Z')
 
 
 
 
 

