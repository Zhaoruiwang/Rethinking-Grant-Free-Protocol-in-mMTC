time_set= 0.01 : 0.01 : 0.1;
n = length(time_set);
m=2001;
sorted_indices_act=zeros(n,m);
sorted_indices_new=zeros(n,m);
sorted_indices_idc=zeros(n,m);
sorted_indices_sel1=zeros(n,m);
sorted_indices_sel2=zeros(n,m);
sorted_indices_sel3=zeros(n,m);

value_act = zeros(n ,1);
value_new=zeros(n,1);
value_idc=zeros(n,1);
value_sel1=zeros(n,1);
value_sel2=zeros(n,1);
value_sel3=zeros(n,1);

for i = 1 : n
    [~, sorted_indices_act(i, :)] = sort(abs(MDE_act(i, :) - FAE_act(i, :)));
    [~,sorted_indices_new(i,:)]= sort(abs(MDE_new(i,:)-FAE_new(i,:)));
    [~,sorted_indices_idc(i,:)]= sort(abs(MDE_idc(i,:)-FAE_idc(i,:)));
    [~,sorted_indices_sel1(i,:)]= sort(abs(MDE_sel1(i,:)-FAE_sel1(i,:)));
    [~,sorted_indices_sel2(i,:)]= sort(abs(MDE_sel2(i,:)-FAE_sel2(i,:)));
    [~,sorted_indices_sel3(i,:)]= sort(abs(MDE_sel3(i,:)-FAE_sel3(i,:)));
    index_act = sorted_indices_act(i, 1);
    value_act(i) = (FAE_act(i, index_act) +MDE_act(i, index_act)) / 2;
    index_new=sorted_indices_new(i,1);
    value_new(i)=(FAE_new(i,index_new)+MDE_new(i,index_new))/2;
    index_idc=sorted_indices_idc(i,1);
    value_idc(i)=(FAE_idc(i,index_idc)+MDE_idc(i,index_idc))/2;
    index_sel1=sorted_indices_sel1(i,1);
    value_sel1(i)=(FAE_sel1(i,index_sel1)+MDE_sel1(i,index_sel1))/2;
    index_sel2=sorted_indices_sel2(i,1);
    value_sel2(i)=(FAE_sel2(i,index_sel2)+MDE_sel2(i,index_sel2))/2;
    index_sel3=sorted_indices_sel3(i,1);
    value_sel3(i)=(FAE_sel3(i,index_sel3)+MDE_sel3(i,index_sel3))/2;

end

figure;
semilogy(cov_time(: ,1), value_new, '-o',LineWidth=1.5);hold on;
semilogy(cov_time(: ,2), value_idc, LineWidth=1.5);hold on;
semilogy(cov_time(: ,3), value_act, '-o',LineWidth=1.5);hold on;
semilogy(cov_time(: ,4), value_sel1, LineWidth=1.5);hold on;
semilogy(cov_time(: ,5), value_sel2, LineWidth=1.5);hold on;
semilogy(cov_time(: ,6), value_sel3, LineWidth=1.5);hold on;
legend('CD','Ideal CD', 'Active CD','D=0, alpha=0.1', 'D=0, alpha=0.01', 'D=2, alpha=0.01');
grid on;

