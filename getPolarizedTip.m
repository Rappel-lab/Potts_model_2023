function tipmap = getPolarizedTip(nsig, ntyp, Px, Py, cmx, cmy, x, y, nbxupd, nbyupd, MeshN)
    tip_percent_thres = 0.99;
    tipmap = zeros(size(nsig),'single');
    ncell = length(Px);
    pol_ang = -1*atan2d(Py,Px); 
        
    for k = 1:ncell
        cellang = -1*(atan2d((y - cmy(k)),(x - cmx(k))));
        celldist = sqrt((y - cmy(k)).^2 + (x - cmx(k)).^2);
     
        pol_metric = celldist.*cos((pol_ang(k) - cellang)*pi()/180).*(nsig == k);
        
        tipmap = tipmap + (pol_metric >= max(pol_metric(:)*tip_percent_thres));
    end
    
    % check if tip is toucning ECM in nearest 8 neighbors
    tipidx = find(tipmap);
    touchingECM = sum(ntyp(tipidx+nbxupd*MeshN+nbyupd) == 2, [2]) > 0;
    tipmap(tipidx(~touchingECM)) = 0;

end