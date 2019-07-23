clear
num_p=4;  %���ֻ���
k=3;  %3��UV

%�������
%FC�ĺ�������[0,15]��DC�ĺ�������[30,55],C�ĺ�������[80,120]
FC_cor = [13 60];
DC_cor = [37,24;46,86];
customer_cor = [72,55;95,25;119,83;105,55;96,66;88,99;97,118;83,109;109,100;78,120;100,4];


% FC_cor = [0 500];
% DC_cor = [500,1000;500,0];
% customer_cor = [1000,0;1000,100;1000,1000;1000,1100;1000,1200];
num_FC=size(FC_cor,1);
num_DC=size(DC_cor,1);
num_customer=size(customer_cor,1);

scatter(FC_cor(:,1),FC_cor(:,2),'r');
hold on;
scatter(DC_cor(:,1),DC_cor(:,2),'b');
hold on;
scatter(customer_cor(:,1),customer_cor(:,2),'g');
% hold on;
%gplot���Ի����ڽӾ���


%���߱���
%���߱���������Ӧ��Ϊһ������,����A1��A2

%u_ij
ss = 2;
u=intvar(ss,k);

%��ͨ�����
%x_ij_k  k=3
num_A1 = num_FC + num_DC;
x_1=binvar(num_A1,num_A1,'full');  %FC��DC�ľ��߱���
x_2=binvar(num_A1,num_A1,'full');  
x_3=binvar(num_A1,num_A1,'full');  

%y_ij_s  s=4
num_A2 = 1 + num_customer;
y_1=binvar(num_A2,num_A2,'full'); %DC��customer�ľ��߱���
y_2=binvar(num_A2,num_A2,'full');

%ÿ��DCֻ��1����

% y_3=binvar(num_A2,num_A2,'full');
% y_4=binvar(num_A2,num_A2,'full');

%capacity���
%w_ij_kp
% k=3 p=4 
w_11=intvar(num_A1,num_A1,'full'); w_12=intvar(num_A1,num_A1,'full');
w_13=intvar(num_A1,num_A1,'full'); w_14=intvar(num_A1,num_A1,'full');
w_21=intvar(num_A1,num_A1,'full'); w_22=intvar(num_A1,num_A1,'full');
w_23=intvar(num_A1,num_A1,'full'); w_24=intvar(num_A1,num_A1,'full');
w_31=intvar(num_A1,num_A1,'full'); w_32=intvar(num_A1,num_A1,'full');
w_33=intvar(num_A1,num_A1,'full'); w_34=intvar(num_A1,num_A1,'full');

%g_ij_s  s=4 
g_1=intvar(num_A2,num_A2,'full');
g_2=intvar(num_A2,num_A2,'full');
% g_3=intvar(num_A2,num_A2,'full');
% g_4=intvar(num_A2,num_A2,'full');

%u_sk ����k��DCs��λ�� 

%������� ���� ���뼴Ϊcost
% FCΪǰnum_FC�У����Ӧ��distanceΪǰ����
A1 = [FC_cor;DC_cor];
dis_A1 = squareform(pdist(A1));

%DCΪǮnum_DC�У����Ӧ��distanceΪǰ����
A2 = [DC_cor(1,:);customer_cor];
dis_A2 = squareform(pdist(A2));

A3 = [DC_cor(2,:);customer_cor];
dis_A3 = squareform(pdist(A3));

%���Լ��FC��DC֮�䲻�ܻ�����ͨ������Ҫʹ�ò���distanceΪ����


%Demand ���
rand('seed',1);
demand = round(rand(num_customer,num_p)*5);  %�У��˿�   �У���������
rand('seed',5);
pi_p = round(rand(num_p,1)*3);  %ÿ���͵Ļ����capacity��Ӱ��
Q1 = 60;
Q2 = 40;
%��������
hs = 10;  %

%�������
I_f1 = 40;
I_f2 = 40;
I_f3 = 40;
I_f4 = 40;

%Ŀ�꺯��
% Objective = sum(sum(x_1*dis_A1+x_2*dis_A1+x_3*dis_A1)) + sum(sum(y_1*dis_A2+y_2*dis_A2+y_3*dis_A2+y_4*dis_A2));
Objective = 1*sum(sum(x_1.*dis_A1+x_2.*dis_A1+x_3.*dis_A1)) + 3*sum(sum(y_1.*dis_A2+y_2.*dis_A3));
%ע��˴���ת�ã���variable����ת��
%��ʱ�����ǰ��˳ɱ�

CX = [];

%Լ��2
for j =1:num_A1
    CX=CX+(sum(x_1(:,j)) == sum(x_1(j,:)));
    CX=CX+(sum(x_2(:,j)) == sum(x_2(j,:)));
    CX=CX+(sum(x_3(:,j)) == sum(x_3(j,:)));
end

%Լ��3&4
for i=1:num_A1   %FC --> DC�� j is DC
   CX=CX+(sum(x_1(:,i))<=1);
   CX=CX+(sum(x_2(:,i))<=1);
   CX=CX+(sum(x_3(:,i))<=1);
end

%Լ��5 ��FC�ڵ㣬��x�ĵ�һ�У� ����ֻ��һ��FC����Լ��5��Լ��6һ��
%��Լ������С��ģ������Է��� 
CX=CX+(sum(x_1(1,:))+sum(x_2(1,:))+sum(x_3(1,:)) <= k);  %С���ܵĳ�����

%Լ��7 �Ƚϸ��ӣ�ע��ṹ
%w_ij_kp
% k=3 p=4 
for i=1:num_A1
   for j=1:num_A1
      %k=1ʱ 
      CX=CX+(pi_p(1)*w_11(i,j) + pi_p(2)*w_12(i,j) + pi_p(3)*w_13(i,j) + pi_p(4)*w_14(i,j)<=x_1(i,j)* Q1);
      %k=2ʱ
      CX=CX+(pi_p(1)*w_21(i,j) + pi_p(2)*w_22(i,j) + pi_p(3)*w_23(i,j) + pi_p(4)*w_24(i,j)<=x_2(i,j)* Q1);
      %k=3ʱ
      CX=CX+(pi_p(1)*w_31(i,j) + pi_p(2)*w_32(i,j) + pi_p(3)*w_33(i,j) + pi_p(4)*w_34(i,j)<=x_3(i,j)* Q1);
   end
end

%Լ��8  ������
% for s=2:num_A1  %�����Ǵ�2��ʼ����FD����������
    %w_ij_kp    k=3 p=4 
    %p��4��  |  ��ʽ�Ҳ࣬j����C�����ӵ�3�п�ʼ����  |  ����demand  �У��˿�   �У���������  | ��ÿ��k
%     CX=CX+( sum(w_11(:,2)+ w_12(:,2)+ w_13(:,2) + w_14(:,2)) - sum(w_11(2,:)+ w_12(2,:)+ w_13(2,:) + w_14(2,:)) == sum(sum(y_1(3:end,3:end)*demand)));
%     CX=CX+( sum(w_21(:,2)+ w_22(:,2)+ w_23(:,2) + w_24(:,2)) - sum(w_21(2,:)+ w_22(2,:)+ w_23(2,:) + w_24(2,:)) == sum(sum(y_1(3:end,3:end)*demand)));
%     CX=CX+( sum(w_11(:,3)+ w_12(:,3)+ w_13(:,3) + w_14(:,3)) - sum(w_11(3,:)+ w_12(3,:)+ w_13(3,:) + w_14(3,:)) == sum(sum(y_2(3:end,3:end)*demand)));
%     CX=CX+( sum(w_21(:,3)+ w_22(:,3)+ w_23(:,3) + w_24(:,3)) - sum(w_21(3,:)+ w_22(3,:)+ w_23(3,:) + w_24(3,:)) == sum(sum(y_2(3:end,3:end)*demand)));

    CX=CX+(sum(w_11(:,2))+sum(w_21(:,2))+sum(w_31(:,2))-(sum(w_11(2,:))+sum(w_21(2,:))+sum(w_31(2,:))) == sum(y_1*[0;demand(:,1)]));
    CX=CX+(sum(w_12(:,2))+sum(w_22(:,2))+sum(w_32(:,2))-(sum(w_12(2,:))+sum(w_22(2,:))+sum(w_32(2,:))) == sum(y_1*[0;demand(:,2)]));
    CX=CX+(sum(w_13(:,2))+sum(w_23(:,2))+sum(w_33(:,2))-(sum(w_13(2,:))+sum(w_23(2,:))+sum(w_33(2,:))) == sum(y_1*[0;demand(:,3)]));
    CX=CX+(sum(w_14(:,2))+sum(w_24(:,2))+sum(w_34(:,2))-(sum(w_14(2,:))+sum(w_24(2,:))+sum(w_34(2,:))) == sum(y_1*[0;demand(:,4)]));
   
    CX=CX+(sum(w_11(:,3))+sum(w_21(:,3))+sum(w_31(:,3))-(sum(w_11(3,:))+sum(w_21(3,:))+sum(w_31(3,:))) == sum(y_2*[0;demand(:,1)]));
    CX=CX+(sum(w_12(:,3))+sum(w_22(:,3))+sum(w_32(:,3))-(sum(w_12(3,:))+sum(w_22(3,:))+sum(w_32(3,:))) == sum(y_2*[0;demand(:,2)]));
    CX=CX+(sum(w_13(:,3))+sum(w_23(:,3))+sum(w_33(:,3))-(sum(w_13(3,:))+sum(w_23(3,:))+sum(w_33(3,:))) == sum(y_2*[0;demand(:,3)]));
    CX=CX+(sum(w_14(:,3))+sum(w_24(:,3))+sum(w_34(:,3))-(sum(w_14(3,:))+sum(w_24(3,:))+sum(w_34(3,:))) == sum(y_2*[0;demand(:,4)]));

% end
 
%Լ��9  ����FC������P
%w_ij_kp
% k=3 p=4 
CX=CX+(sum(w_11(1,:))+sum(w_21(1,:))+sum(w_31(1,:))<= I_f1);
CX=CX+(sum(w_12(1,:))+sum(w_22(1,:))+sum(w_32(1,:))<= I_f2);
CX=CX+(sum(w_13(1,:))+sum(w_23(1,:))+sum(w_33(1,:))<= I_f3);
CX=CX+(sum(w_14(1,:))+sum(w_24(1,:))+sum(w_34(1,:))<= I_f4);

%Լ��10 s ��2��ʼ
for s=2:num_A1
   CX=CX+(sum(w_11(:,s))>=(sum(w_11(s,:))));  CX=CX+(sum(w_12(:,s))>=(sum(w_12(s,:))));  CX=CX+(sum(w_13(:,s))>=(sum(w_13(s,:)))); CX=CX+(sum(w_14(:,s))>=(sum(w_14(s,:))));
   CX=CX+(sum(w_21(:,s))>=(sum(w_11(s,:))));  CX=CX+(sum(w_22(:,s))>=(sum(w_12(s,:))));  CX=CX+(sum(w_23(:,s))>=(sum(w_13(s,:)))); CX=CX+(sum(w_24(:,s))>=(sum(w_14(s,:))));
   CX=CX+(sum(w_31(:,s))>=(sum(w_11(s,:))));  CX=CX+(sum(w_32(:,s))>=(sum(w_12(s,:))));  CX=CX+(sum(w_33(:,s))>=(sum(w_13(s,:)))); CX=CX+(sum(w_34(:,s))>=(sum(w_14(s,:))));
end

%Լ��11
for i=2:num_A1
    for j=2:num_A1
        CX=CX+(u(i-1,2)+1 <= u(j-1,2) + 100000*(1-x_1(i,j)));
        CX=CX+(u(i-1,3)+1 <= u(j-1,3) + 100000*(1-x_2(i,j)));
    end
end

%Լ��12  �˴�����ڶ��㣬��s������dis_A2��ʽ�ӣ��ɴ���������ģ�ͱ����ϲ�̫һ�� ��ʵ��������һ�µ� j��3��ʼΪcustomer
for j=2:num_A2
    CX=CX+(sum(y_1(:,j))+sum(y_2(:,j))==1);
end

%Լ��13 ��ʵ���Ժ�Լ��12д��һ��
for j=2:num_A2
   CX=CX+(sum(y_1(:,j)) == sum(y_1(j,:))); 
   CX=CX+(sum(y_2(:,j)) == sum(y_2(j,:)));
end

%Լ��14 ��ʱ����

   CX=CX+(sum(y_1(1,:))<=2);
   CX=CX+(sum(y_2(1,:))<=2);


%����Լ�� ���಻���� ?�ӻ�·��
% CX=CX+(y_1(:,2) == 0); CX=CX+(y_1(2,:) == 0);
% CX=CX+(y_2(:,1) == 0); CX=CX+(y_2(1,:) == 0);


%Լ��15 �˴�sӦ��Ϊ1��2 ��Ϊ�������˱仯 
CX=CX+(sum(y_1(1,:))+sum(y_2(1,:))<=4);

%Լ��16  demand �У��˿�   �У���������
for j=2:num_A2
   CX=CX+(sum(g_1(:,j))+sum(g_2(:,j)) - (sum(g_1(j,:))+sum(g_2(j,:))) == demand(j-1,:)*pi_p);
end

%Լ��17
for i=1:num_A2
    for j=1:num_A2
       CX=CX+( g_1(i,j)<=Q2*y_1(i,j));
       CX=CX+( g_2(i,j)<=Q2*y_2(i,j));
    end
end

%��������Լ��
for i=1:num_A1
    for j=1:num_A1
       CX=CX+(w_11(i,j)>=0); CX=CX+(w_12(i,j)>=0); CX=CX+(w_13(i,j)>=0); CX=CX+(w_14(i,j)>=0);
       CX=CX+(w_21(i,j)>=0); CX=CX+(w_22(i,j)>=0); CX=CX+(w_23(i,j)>=0); CX=CX+(w_24(i,j)>=0);
       CX=CX+(w_31(i,j)>=0); CX=CX+(w_32(i,j)>=0); CX=CX+(w_33(i,j)>=0); CX=CX+(w_34(i,j)>=0);
    end
end

for i=1:num_A2
    for j =1:num_A2
        CX=CX+(g_1(i,j)>=0);
        CX=CX+(g_2(i,j)>=0);
    end
end

options=sdpsettings('verbose',1,'solver','cplex','cplex.TimeLimit',1000);
sol=optimize(CX,Objective,options);

Objective_v=value(Objective);

x_1_v = value(x_1);
x_2_v = value(x_2);
x_3_v = value(x_3);
y_1_v = value(y_1);
y_2_v = value(y_2);
g_1_v = value(g_1);
g_2_v = value(g_2);

w_11_v = value(w_11);
w_12_v = value(w_12);
w_13_v = value(w_13);
w_14_v = value(w_14);

w_21_v = value(w_21);
w_22_v = value(w_22);
w_23_v = value(w_23);
w_24_v = value(w_24);

w_31_v = value(w_31);
w_32_v = value(w_32);
w_33_v = value(w_33);
w_34_v = value(w_34);


gplot(x_1_v,A1,'b-*');
hold on
gplot(x_2_v,A1,'b-*');
hold on
gplot(x_3_v,A1,'b-*');
hold on
gplot(x_3_v+1,A1,'b-*');
hold on
gplot(y_1_v,A2,'r-o');
hold on
gplot(y_2_v,A3,'g-o');

