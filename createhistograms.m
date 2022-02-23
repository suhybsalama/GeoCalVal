function createhistograms(data,nbins,x_name,x_limi)



% Create axes

axes1 = gca;
set(gca,'FontSize',18);

xlim(axes1,x_limi);
box(axes1,'on');
hold(axes1,'all');

% Create bar

h = histfit(data,nbins,'tlocationscale') ;
set(h(1),'FaceColor',[0.5,0.5,0.5],...
    'EdgeColor',[0.5,0.5,0.5]);

set(h(2),'Color',[0 0 0 ]);
xlabel(x_name,'FontSize',18);

ylabel('Frequency','FontSize',18);
