function cancer_parapair_commonplot(par1plotname, par2plotname, xtnew, xtlbl, ytnew, ytlbl, bp, bpstatus)
    ylabel(par1plotname,'FontSize',12,'interpreter','latex')
    xlabel(par2plotname,'FontSize',12,'interpreter','latex')
    set(gca,'YDir','normal')
    set(gca, 'XTick', xtnew, 'XTickLabel', xtlbl)
    set(gca, 'YTick', ytnew, 'YTickLabel', ytlbl)
    
if bpstatus
    plot(bp.xvalid_rotcell(bp.k_rotcell),bp.yvalid_rotcell(bp.k_rotcell),'r--','LineWidth',2);
    plot(bp.xvalid_freecell(bp.k_freecell),bp.yvalid_freecell(bp.k_freecell),'g--','LineWidth',2);
    plot(bp.xvalid_onecell(bp.k_onecell),bp.yvalid_onecell(bp.k_onecell),'k--','LineWidth',2);
end
end