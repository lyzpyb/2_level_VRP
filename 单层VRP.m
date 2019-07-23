clc;clear;
vNum=4;  %������
cityNum=12; %��������
vCapacity=40; %������
num_p=4;
%demands=[0,1,2,3,2,3]; %������
x=[13,72,95,119,105,96,88,97,83,109,78,100];
y=[78,55,25,83,55,66,99,118,109,100,120,4];
axis=[x' y']; %��������

rand('seed',1);
demands =[[0;0;0;0],round(rand(num_p,cityNum-1)*5)];  %�У��˿�   �У���������
rand('seed',5);
pi_p = round(rand(num_p,1)*3);  %ÿ���͵Ļ����capacity��Ӱ��


Cij = zeros(cityNum);%�������֮��ľ���
for i=1:(cityNum)
    for j=1:i
        x1 = axis(i,1);
        y1 = axis(i,2);
        x2 = axis(j,1);
        y2 = axis(j,2);
        Cij(i,j)=sqrt((x1-x2)^2+(y1-y2)^2);
        Cij(j,i)=Cij(i,j);
    end
end
%% ���߱���
Xijk=binvar(cityNum,cityNum,vNum,'full');
Yik=binvar(cityNum,vNum,'full');
Uik=sdpvar(cityNum,vNum,'full'); %Uik>demands(i),ȥ�ڵ�i��Uik��ʾǰi-1���ڵ�װ�ú󣬳��ӵ�ʣ������
%% Ŀ�꺯��
obj=0;
for i=1:cityNum
    for j=1:cityNum
        for k=1:vNum
            obj=obj+5*Cij(i,j)*Xijk(i,j,k);
        end
    end
end
f=obj;
%% Լ������
F=[];
for i=2:cityNum
    F=[F;sum(Yik(i,:))==1]; %ÿ�������i���ᱻһ��С������
end
for i=1
    F=[F;sum(Yik(i,:))==vNum];%����������ᱻ�����õ���С������
end
d1=(pi_p(1)*demands(1,:));
d2=(pi_p(2)*demands(2,:));
d3=(pi_p(3)*demands(3,:));
for k=1:vNum
    F=[F;sum(d1(:).*Yik(:,k))+sum(d2(:).*Yik(:,k))+sum(d3(:).*Yik(:,k))<=vCapacity]; %ÿ����·�ϵ�������֮��С�ڳ�������
end
for i=1:cityNum
    for j=1:cityNum
        for k=1:vNum
            if i==j
                F=[F;Xijk(i,j,k)==0];
            end
        end
    end
end
for i=1:cityNum
    for k=1:vNum
        F=[F;sum(Xijk(i,:,k))==sum(Xijk(:,i,k))];%����ƽ��
    end
end
for j=2:cityNum
    for k=1:vNum
        F=[F;sum(Xijk(:,j,k))==Yik(j,k)];%Xijk��Yik�Ĺ�ϵ
    end
end
for i=1:cityNum
    for k=1:vNum
        F=[F;sum(Xijk(i,:,k))==Yik(i,k)];%Xijk��Yik�Ĺ�ϵ
    end
end
for i=2:cityNum
    for j=2:cityNum
        for k=1:vNum
            if i~=j
                if sum(demands(:,i))+sum(demands(:,j))<=vCapacity
                    F=[F;Uik(i,k)-Uik(j,k)+vCapacity*Xijk(i,j,k)<=vCapacity-sum(demands(:,i))];%Xijk��Ui
                end
            end
        end
    end
end
for i=2:cityNum
    for k=1:vNum
        F=[F;Uik(i,k)<=vCapacity];
        F=[F;Uik(i,k)>=sum(demands(:,i))];
    end
end
%% ���
ops = sdpsettings( 'solver','cplex');
sol=solvesdp(F,f,ops);
f=value(f);
Xijk=value(Xijk);
Yik=double(Yik);
Uik=double(Uik);
%% ��ͼ
plot(axis(2:cityNum,1),axis(2:cityNum,2),'ro');hold on;
plot(axis(1,1),axis(1,2),'pm');hold on;
for i=1:cityNum
    for j=1:cityNum
        for k=1:vNum
            if Xijk(i,j,k)==1
                plot([axis(i,1),axis(j,1)],[axis(i,2),axis(j,2)],'-');
            end
        end
    end
end
for k=1:vNum    
[a,b]=find(Xijk(:,:,k));
sqe=[a,b];
sqe1=zeros(1,0);
sqe1(1)=1;
[a,b]=find(sqe(:,1)==1);
for i=2:length(sqe)+1
    [a,b]=find(sqe(:,1)==sqe1(i-1));
    sqe1(i)=sqe(a,b+1);   
end
disp(['����',num2str(k),'��·�����£�']);
disp(sqe1)
value(obj)
end
