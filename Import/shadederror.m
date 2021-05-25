function shadederror(x,d,bottombar,topbar,name)

xx=[x,fliplr(x)];
ss=[topbar ,fliplr(bottombar)];
xx=xx(~isnan(ss));
ss=ss(~isnan(ss));
hold on
h=plot(x,d,'-','DisplayName',name);
color=get(h,'Color');
patch(xx, ss, color,'linestyle', 'none','FaceAlpha','0.5');
hold off
box on;
end