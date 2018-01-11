load BANDS.RESULT
load FERMI

BANDS=BANDS-FERMI;

hold off
for jj=1:2:size(BANDS,2)
 plot(BANDS(:,jj),'LineWidth',2,'-r');
 plot(BANDS(:,jj+1),'LineWidth',2,'-b');
 hold on
end

ylim([-1.5 1.5])

set(gca,'xtick',1:49:size(BANDS,1))
set(gca,'xticklabel',{'GM','M','K', 'GM','A','L','H','A'})
xlim([1 size(BANDS,1)])

grid on
xlabel('k-vector (Spacegroup #189: P-62m')
ylabel('Energy (eV)')

print -deps -mono -solid zoom.eps
