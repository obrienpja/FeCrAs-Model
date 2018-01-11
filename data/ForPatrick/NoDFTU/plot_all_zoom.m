load BANDS.RESULT
load FERMI

BANDS=BANDS-FERMI;

hold off
for jj=1:size(BANDS,2)
 plot(BANDS(:,jj),'LineWidth',2,'-k');
 hold on
end

ylim([-0.4 0.4])

set(gca,'xtick',1:49:size(BANDS,1))
set(gca,'xticklabel',{'GM','M','K', 'GM','A','L','H','A'})
xlim([1 size(BANDS,1)])

grid on
xlabel('k-vector (Spacegroup #189: P-62m')
ylabel('Energy (eV)')

print -deps -mono -solid zoom.eps
