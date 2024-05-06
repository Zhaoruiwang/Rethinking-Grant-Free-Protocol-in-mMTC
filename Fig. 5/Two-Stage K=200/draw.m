function draw(time_set, m, MDE_new, FAE_new, MDE_sel, FAE_sel)
n=length(time_set);
sorted_indices_new=zeros(n,m);
% sorted_indices_idc=zeros(n,m);
% sorted_indices_prf=zeros(n,m);
sorted_indices_sel=zeros(n,m);

value_new=zeros(n,1);
% value_idc=zeros(n,1);
% value_prf=zeros(n,1);
value_sel=zeros(n,1);

for i=1:n
    [~,sorted_indices_new(i,:)]= sort(abs(MDE_new(i,:)-FAE_new(i,:)));
    % [~,sorted_indices_idc(i,:)]= sort(abs(MDE_idc(i,:)-FAE_idc(i,:)));
    % [~,sorted_indices_prf(i,:)]= sort(abs(MDE_prf(i,:)-FAE_prf(i,:)));
    [~,sorted_indices_sel(i,:)]= sort(abs(MDE_sel(i,:)-FAE_sel(i,:)));
    index_new=sorted_indices_new(i,1);
    value_new(i)=(FAE_new(i,index_new)+MDE_new(i,index_new))/2;
    % index_idc=sorted_indices_idc(i,1);
    % value_idc(i)=(FAE_idc(i,index_idc)+MDE_idc(i,index_idc))/2;
    % index_prf=sorted_indices_prf(i,1);
    % value_prf(i)=(FAE_prf(i,index_prf)+MDE_prf(i,index_prf))/2;
    index_sel=sorted_indices_sel(i,1);
    value_sel(i)=(FAE_sel(i,index_sel)+MDE_sel(i,index_sel))/2;

end

figure;
p1 = semilogy(time_set,value_new,'o-', LineWidth=2);hold on;
% semilogy(time_set,value_idc,'o-');hold on;
% semilogy(time_set,value_prf);hold on;
p2 = semilogy(time_set,value_sel,'o-', LineWidth=2);hold on;
legend([p1 p2], {'Random CD','Select CD with $L_{\mathrm{I}}$=5'}, 'Interpreter', 'latex');
xlabel('CPU time');
ylabel('Detection error probability');
grid on;
end