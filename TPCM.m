function [C]=TPCM(x,y,lag,tau,dim)    %%%%%%%%%%%%用我们自己的相空间重构方法,用余弦距离

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%我们提出的时移排序交叉映射方法（time-shift permutation cross-mapping， TPCM）
  %%%%%先对数据进行归一化
  %x=normcdf(x);  y=normcdf(y);
%Estimates signal roughness and sets threshold
if nargin<6
    roughness=max(std(x(2:end)-x(1:end-1))/std(x),std(y(2:end)-y(1:end-1))/std(y));  
        thresh=0.05*roughness;
end
 
%Offsets time series according to the specified
%lag,这里延迟进行改变，由原来相反，这里的延迟至适用于X->Y
if lag>=0
    y=y(lag+1:end);
    x=x(1:end-lag);
else
    y=y(1:end+lag);
    x=x(1-lag:end);
end
 
%Creates phase space embedding
tau=1;  dim=2;
%Creates phase space embedding
L=length(x);
T=1+(dim-1)*tau;
Le=L-T+1;


[X]=myeig(x);
[Y]=myeig(y);
 
%Computes distance matrices，用余弦距离
dX=pdist2(X,X,'seuclidean');
dY=pdist2(Y,Y,'seuclidean');
 
 
%Removes the central diagonal with trivial neighbors
Xdiagstd=zeros(Le,1);
Ydiagstd=zeros(Le,1);
for i=0:Le-1
    Xdiagstd(i+1)=std(diag(dX,i));
    Ydiagstd(i+1)=std(diag(dY,i));
end

tmask=min(find(Xdiagstd>=nanmean(Xdiagstd),1),find(Ydiagstd>=nanmean(Ydiagstd),1));
mask=triu(ones(Le),tmask);
sm=numel(find(mask));
 
%Flattens and sorts the distance matrices
yflat=reshape(dY(mask==1),sm,1);
xflat=reshape(dX(mask==1),sm,1);


[~,idxs]=sort(xflat);
[~,idys]=sort(yflat);
rank=[1:sm].'/sm;
xr(idxs,1)=rank;
yr(idys,1)=rank;
 
%Resorts the ranks according to their order in the other reconstruction
r_xcy=xr(idys);
r_ycx=yr(idxs);
 
 
%Finds the square error of the resorted ranks 
se_xcy=(r_xcy-rank).^2;
se_ycx=(r_ycx-rank).^2;
 
 
%Normalizes the error by the expected null error 
se_null=rank.^2-rank+1/3;
diff_se_xcy=se_xcy-se_null;
diff_se_ycx=se_ycx-se_null;
 
tm=round(sm*thresh);
xtemp=[1:1:tm].'*thresh/tm;
count=[1:1:tm].';
 
rvp_xcy=cumsum(-1*diff_se_xcy(1:tm)./se_null(1:tm))./count;
rvp_ycx=cumsum(-1*diff_se_ycx(1:tm)./se_null(1:tm))./count;
 
%n=max(round(length(rvp_xcy)/200),1);    %%%%%%%%%%%%%%降采样处理
 n=1;
drvp_xcy=downsample(rvp_xcy,n);
drvp_ycx=downsample(rvp_ycx,n);
dxtemp=downsample(xtemp,n);
dcount=downsample(count,n);


%%%%%%%%%%%%%%%%%%%%%%新加的一个判断
% if (abs(drvp_xcy(1))>=abs(2*drvp_xcy(2))) 
%    drvp_xcy(1)=drvp_xcy(2);
% else
%     if (abs(2*drvp_xcy(1))<=abs(drvp_xcy(2)))
%    drvp_xcy(1)=drvp_xcy(2);     
%     end
% end
% 
% 
% if (abs(drvp_ycx(1))>=abs(2*drvp_ycx(2))) 
%    drvp_ycx(1)=drvp_ycx(2);
% else
%     if (abs(2*drvp_ycx(1))<=abs(drvp_ycx(2)))
%      drvp_ycx(1)=drvp_ycx(2);    
%     end
% end

if (abs(drvp_xcy(1))>=abs(2*drvp_xcy(2))) || (abs(2*drvp_xcy(1))<=abs(drvp_xcy(2)))
   drvp_xcy(1)=drvp_xcy(2);
end

if (abs(drvp_ycx(1))>=abs(2*drvp_ycx(2))) || (abs(2*drvp_ycx(1))<=abs(drvp_ycx(2)))
   drvp_ycx(1)=drvp_ycx(2);
end

%Fit scaled EER^2
expfit = fittype('a + b*exp(c*x)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b','c'});
 
rvp_xcyfit=fit(dxtemp,drvp_xcy,expfit,'StartPoint',[0,drvp_xcy(1),0],'weights',sqrt(dcount));
rvp_ycxfit=fit(dxtemp,drvp_ycx,expfit,'StartPoint',[0,drvp_ycx(1),0],'weights',sqrt(dcount));





%Returns the Y intercepts of the fitted curves
C=[rvp_xcyfit(0);rvp_ycxfit(0)];
C(C>1)=1;C(C<-1)=-1;
 
 
%plots
if nargout==0
    
    try
        run('Figure_Prefs.m')
    catch
    end
    
    ns=200;
    points=round(linspace(1,sm,ns+1)).';
    mse_null=zeros(ns-1,1);
    mse_xcy=zeros(ns-1,1);
    mse_ycx=zeros(ns-1,1);
 
    for i=1:ns
        mse_null(i)=mean(se_null(points(i):points(i+1)));
        mse_xcy(i)=mean(se_xcy(points(i):points(i+1)));
        mse_ycx(i)=mean(se_ycx(points(i):points(i+1)));
    end
    
    
    figure()
    tiledlayout(1,2)
    xp=points(1:ns)/(sm);
    
    nexttile()
    plot(xp,mse_xcy)
    hold on
    plot(xp,mse_ycx)
    hold on
    plot(xp,mse_null,'g')
    legend('x->y','y->x','null')
    xlabel('Rank')
    ylabel('ERR^2')
    
    
    
    nexttile()
    plot(xtemp,rvp_xcy)
    hold on
    plot(xtemp,rvp_xcyfit(xtemp))
    hold on
    plot(xtemp,rvp_ycx)
    hold on
    plot(xtemp,rvp_ycxfit(xtemp))
    legend('x->y','x->y fit','y->x','y->x fit')
    xlabel('Rank')
    ylabel('E[R^2]/E[R^2 null]')
    
    
end
end
 





function [y]=myeig(data1)
%输入为原始序列，然后创建维度为2，延迟为1的重构向量，最后输出新的重构矩阵

%C=xT*x;   然后对C特征分解，得到fi=[],每一列都是特征向量

for i=1:length(data1)-1
    x(i,1)=data1(i);
    x(i,2)=data1(i+1);
end

C=x'*x;

[fi,d]=eig(C);
y=x*fi;

end