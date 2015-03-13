function trianglePlot(LH_matrix,cellLines)
	load mycmap
	figure;
	hold on;
	plot(nan);
	xlim([0 1]);
	ylim([0 1]);
	a = [0,0,0];
	b = [1,0,0];
	c = [0.5,sqrt(3)/2,0];
	d = [0.5,sqrt(3)/6,sqrt(2)/sqrt(3)];
	m = (a+b+c+d)/4;
	fill3([a(1),b(1),c(1)],[a(2),b(2),c(2)],[a(3),b(3),c(3)],[1 1 1]);
	fill3([a(1),b(1),d(1)],[a(2),b(2),d(2)],[a(3),b(3),d(3)],[1 1 1]);
	fill3([a(1),c(1),d(1)],[a(2),c(2),d(2)],[a(3),c(3),d(3)],[1 1 1]);
	fill3([b(1),c(1),d(1)],[b(2),c(2),d(2)],[b(3),c(3),d(3)],[1 1 1]);
	axis off;
	% E->F<-A
	% E->A->F
	% E->A->F<-E
	% A->F
	text(-0.2,0,0,'E$\rightarrow$F$\leftarrow$A','interpreter','latex','fontsize',12);
	text(1.01,0,0,'E$\rightarrow$A$\rightarrow$F','interpreter','latex');
	text(0.46,0.866,0,'E$\rightarrow$F$\leftarrow$A$\leftarrow$E','interpreter','latex');
	text(0.5,sqrt(3)/6,sqrt(2)/sqrt(3),'A$\rightarrow$F','interpreter','latex');
	% plot3([a(1),m(1)],[a(2),m(2)],[a(3),m(3)],'k--');
	% plot3([b(1),m(1)],[b(2),m(2)],[b(3),m(3)],'k--');
	% plot3([c(1),m(1)],[c(2),m(2)],[c(3),m(3)],'k--');
	% plot3([d(1),m(1)],[d(2),m(2)],[d(3),m(3)],'k--');
	alpha(0)
	axis square;
	X = LH_matrix(:,[2,3,4,1]);
	for i = 1:size(LH_matrix,1)
		w = exp(X(i,:)-mean(X(i,:)));
		w = w./sum(w);
		x = w*[a;b;c;d];
		plot3(x(1),x(2),x(3),'o','markersize',10,'markeredgecolor','k','markerfacecolor',mycmap(length(mycmap)-floor(x(3)*length(mycmap)),:));
		text(x(1),x(2),x(3),cellLines{i})
	end
end

