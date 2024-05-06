function createfigure(xvector1, ymatrix1)
%CREATEFIGURE(xvector1, ymatrix1)
%  XVECTOR1:  bar xvector
%  YMATRIX1:  bar 矩阵数据

%  由 MATLAB 于 15-Apr-2024 00:06:32 自动生成

% 创建 figure
figure('OuterPosition',[1032 623 576 514]);

% 创建 axes
axes1 = axes;
hold(axes1,'on');

% 使用 bar 的矩阵输入创建多行
bar1 = bar(xvector1,ymatrix1,'BarWidth',4);
set(bar1(2),'DisplayName','Number of coordinates','FaceColor',[0.690196096897125 0.690196096897125 0.690196096897125]);
set(bar1(1),'DisplayName','Number of gradients','FaceColor',[0.403921574354172 0.635294139385223 0.639215707778931]);

% 下面一行演示创建数据提示的另一种方法。
% datatip(bar1(2),2.74285714285714,4976.7644);
% 创建 datatip
datatip(bar1(2),'DataIndex',3);

% 创建 text
text('VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12,'FontName','Times New Roman','String','17637.6',...
    'Position',[0.948686635944701 17918.4833819242 0]);

% 创建 text
text('VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12,'FontName','Times New Roman','String','2824.7',...
    'Position',[1.95905529953917 2918.9944606414 0]);

% 创建 text
text('VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12,'FontName','Times New Roman','String','4976.8',...
    'Position',[2.7347465437788 5350.9778425656 0]);

% 创建 text
text('VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12,'FontName','Times New Roman','String','29262',...
    'Position',[3.26905529953917 29263 0]);

% 创建 text
text('VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12,'FontName','Times New Roman','String','1027.9',...
    'Position',[4.03509216589862 1308.7833819242 0]);

% 创建 text
text('VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12,'FontName','Times New Roman','String','6129.2',...
    'Position',[4.55592165898618 6316.7889212828 0]);

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[0.5 5]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[0 30000]);
box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'FontName','Times New Roman','FontSize',12,'XTick',[0.95 1.95 3 4.3],'XTickLabel',{'CD','Ideal CD','Active CD','FCD'},'YGrid','on');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.17738095956544 0.759126986610514 0.355357134980815 0.101190473494076],'FontSize',12);

