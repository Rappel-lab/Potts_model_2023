function dE = energies(p, xran, yran, xnb, ynb, MeshN, nsig, ntyp, nar, nl, cmx, cmy, velx, vely, Px, Py, conx, cony, nbxupd, nbyupd, S, x ,y, targetarea)
% 	[enold, ennew, epold, epnew, etold, etnew, eoldtot, enewtot] = deal(0);
    enold = 0; ennew = 0;
	nold = ntyp(yran,xran);
	nnew = ntyp(ynb,xnb);

	% First term (adhesion/cohesion/surface tension) %
%     rcsub = [yran+nbyupd;xran+nbxupd]; % row and column subscripts
%     rcsub(rcsub == (MeshN+1)) = 1;
%     rcsub(rcsub == 0) = MeshN;
%     idx = rcsub(1,:) + MeshN*(rcsub(2,:)-1);
%     enold = sum((nsig(yran,xran) ~= nsig(idx)).*(p.jma(nold+1,ntyp(idx)+1)));
%     ennew = sum((nsig(ynb,xnb) ~= nsig(idx)).*(p.jma(nnew+1,ntyp(idx)+1)));
    
    enold = 0;
	for i=1:8
        if nsig(yran,xran) ~= nsig(yran+nbyupd(i),xran+nbxupd(i))
           enold = enold + p.jma(nold+1,ntyp(yran+nbyupd(i),xran+nbxupd(i))+1);
        end
    end
    
    ennew = 0;
    for i=1:8
        if nsig(ynb,xnb) ~= nsig(yran+nbyupd(i),xran+nbxupd(i))
           ennew = ennew + p.jma(nnew+1,ntyp(yran+nbyupd(i),xran+nbxupd(i))+1);
        end
    end

	% Second term (area) %
    if (nold~=1)
		fac2 = nar(nsig(ynb,xnb)) - targetarea(nsig(ynb,xnb));
		etold = fac2*fac2;
		fac2 = nar(nsig(ynb,xnb)) + 1 - targetarea(nsig(ynb,xnb));
		etnew = fac2*fac2;
    else if (nnew~=1)
			fac1 = nar(nsig(yran,xran)) - targetarea(nsig(yran,xran));
			etold = fac1*fac1;
			fac1 = nar(nsig(yran,xran)) - 1 - targetarea(nsig(yran,xran));
			etnew = fac1*fac1;
        else
			fac1= nar(nsig(yran,xran)) - targetarea(nsig(yran,xran));
			fac2 = nar(nsig(ynb,xnb)) - targetarea(nsig(ynb,xnb));
			etold = fac1*fac1 + fac2*fac2;
			fac1 = nar(nsig(yran,xran)) - 1 - targetarea(nsig(yran,xran));
			fac2 = nar(nsig(ynb,xnb)) + 1 - targetarea(nsig(ynb,xnb));
			etnew = fac1*fac1 + fac2*fac2;
        end
    end
% 	if (nold~=1)
% 		fac2 = nar(nsig(ynb,xnb)) - p.area;
% 		etold = fac2*fac2;
% 		fac2 = nar(nsig(ynb,xnb)) + 1 - p.area;
% 		etnew = fac2*fac2;
%     else if (nnew~=1)
% 			fac1 = nar(nsig(yran,xran)) - p.area;
% 			etold = fac1*fac1;
% 			fac1 = nar(nsig(yran,xran)) - 1 - p.area;
% 			etnew = fac1*fac1;
%         else
% 			fac1= nar(nsig(yran,xran)) - p.area;
% 			fac2 = nar(nsig(ynb,xnb)) - p.area;
% 			etold = fac1*fac1 + fac2*fac2;
% 			fac1 = nar(nsig(yran,xran)) - 1 - p.area;
% 			fac2 = nar(nsig(ynb,xnb)) + 1 - p.area;
% 			etnew = fac1*fac1 + fac2*fac2;
%         end
%     end
 
	% Third term (potential, collective velocity alignment) %
    if (nsig(yran,xran) ~= 0)
        epold = conx(nsig(yran,xran))*(xran - cmx(nsig(yran,xran))) + cony(nsig(yran,xran))*(yran - cmy(nsig(yran,xran)));

        if (nsig(ynb,xnb) ~= 0)
            epnew = conx(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + cony(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
        else
            epnew = 0;
        end
    else
        epold = 0;
        epnew = conx(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + cony(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
    end
    
%     % Fourth term (self momentum by velocity) %
%     if (nsig(yran,xran) ~= 0)
%         emnew1 = -velx(nsig(yran,xran))*(xran - cmx(nsig(yran,xran))) - vely(nsig(yran,xran))*(yran - cmy(nsig(yran,xran)));
%         v1 = (velx(nsig(yran,xran))^2 + vely(nsig(yran,xran))^2)^(1/2);
%         if v1 ~= 0
%             emnew1 = emnew1/v1;
%         end
%         
%         if (nsig(ynb,xnb) ~= 0)
%             emnew2 = velx(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + vely(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
%             v2 = (velx(nsig(ynb,xnb))^2 + vely(nsig(ynb,xnb))^2)^(1/2);
%             if v2 ~= 0
%                 emnew2 = emnew2/v2;
%             end
%         else
%             emnew2 = 0;
%         end
%     else
%         emnew1 = 0;
%         emnew2 = velx(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + vely(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
%         v2 = (velx(nsig(ynb,xnb))^2 + vely(nsig(ynb,xnb))^2)^(1/2);
%         if v2 ~= 0
%             emnew2 = emnew2/v2;
%         end
%     end
    
    % Fourth term (momentum by P's) %
    if (nsig(yran,xran) ~= 0)
        emnew1 = -Px(nsig(yran,xran))*(xran - cmx(nsig(yran,xran))) - Py(nsig(yran,xran))*(yran - cmy(nsig(yran,xran)));
        if (nsig(ynb,xnb) ~= 0)
            emnew2 = Px(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + Py(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
        else
            emnew2 = 0;
        end
    else
        emnew1 = 0;
        emnew2 = Px(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + Py(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
    end

    % Fifth term (ECM) %
    eecmold = 0;
% 	for i=1:8
%         if nsig(yran,xran) ~= nsig(yran+nbyupd(i),xran+nbxupd(i))
%            eecmold = eecmold + S(yran+nbyupd(i),xran+nbxupd(i));
%         end
%     end
    if ntyp(yran,xran) == 1
        eecmold = eecmold + S(yran,xran);
    else
        eecmold = 0;
    end
    
    eecmnew = 0;
%     for i=1:8
%         if nsig(ynb,xnb) ~= nsig(yran+nbyupd(i),xran+nbxupd(i))
%            eecmnew = eecmnew + S(yran+nbyupd(i),xran+nbxupd(i));
%         end
%     end
    if ntyp(ynb,xnb) == 1
        eecmnew = eecmnew + S(yran,xran);
    else
        eecmnew = 0;
    end
    
    eperold = 0;
    epernew = 0;
%     % Sixth term (perimeter) %
%     if (nold==0)
% 		fac2 = nl(nsig(ynb,xnb)) - p.perim;
% 		eperold = fac2*fac2;
%         nlnew = nl(nsig(ynb,xnb)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(ynb,xnb))-2)*2;
% 		fac2 = nlnew - p.perim;
% 		epernew = fac2*fac2;
%     else if (nnew==0)
% 			fac1 = nl(nsig(yran,xran)) - p.perim;
% 			eperold = fac1*fac1;
%             nlnew = nl(nsig(yran,xran)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(yran,xran))-2)*(-2);
% 			fac1 = nlnew - p.perim;
% 			epernew = fac1*fac1;
%         else
% 			fac1= nl(nsig(yran,xran)) - p.perim;
% 			fac2 = nl(nsig(ynb,xnb)) - p.perim;
% 			eperold = fac1*fac1 + fac2*fac2;
%             nlnew = nl(nsig(yran,xran)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(yran,xran))-2)*(-2);
% 			fac1 = nlnew - p.perim;
%             nlnew = nl(nsig(ynb,xnb)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(ynb,xnb))-2)*2;
%             fac2 = nlnew - p.perim;
% 			epernew = fac1*fac1 + fac2*fac2;
%         end
%     end
    
    espnew1 = 0;
    espnew2 = 0;
%     % Seventh term (polarization as spring) %
%     if (nsig(yran,xran) ~= 0)
%         espnew1 = Px(nsig(yran,xran))*(xran - cmx(nsig(yran,xran)))^2 + Py(nsig(yran,xran))*(yran - cmy(nsig(yran,xran)))^2;
%         if (nsig(ynb,xnb) ~= 0)
%             espnew2 = Px(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb)))^2 + Py(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)))^2;
%         else
%             espnew2 = 0;
%         end
%     else
%         espnew1 = 0;
%         espnew2 = Px(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb)))^2 + Py(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)))^2;
%     end


    elold1 = 0;
    elold2 = 0;
    elnew1 = 0;
    elnew2 = 0;
    
    % Eighth term (preferred length) %
%     ynew = sind(theta)*xold + cosd(theta)*yold;
%     if (nsig(yran,xran) ~= 0)
%         theta = atan2d(Py(nsig(yran,xran)),Px(nsig(yran,xran)));
%         [row, col] = find(nsig == nsig(yran,xran));
%         xold = col - cmx(nsig(yran,xran));
%         yold = row - cmy(nsig(yran,xran));
%         xtrans = cosd(theta)*xold - sind(theta)*yold;
% %         lenold1 = max(xtrans) - min(xtrans);
%         lenold1 = mean(xtrans);
% 
%         rantrans = cosd(theta)*(xran - cmx(nsig(yran,xran))) - sind(theta)*(yran - cmy(nsig(yran,xran)));
%         xtrans(xtrans == rantrans) = [];
% %         lennew1 = max(xtrans) - min(xtrans);
%         lennew1 = mean(xtrans);
% 
%         elold1 = (lenold1 - p.len)^2;
%         elnew1 = (lennew1 - p.len)^2;
%     else
%         elold1 = 0;
%         elnew1 = 0;
%     end
%            
%     if (nsig(ynb,xnb) ~= 0)
%         theta = atan2d(Py(nsig(ynb,xnb)),Px(nsig(ynb,xnb)));
%         [row, col] = find(nsig == nsig(ynb,xnb));
%         xold = col - cmx(nsig(ynb,xnb));
%         yold = row - cmy(nsig(ynb,xnb));
%         xtrans = cosd(theta)*xold - sind(theta)*yold;
% %         lenold2 = max(xtrans) - min(xtrans);
%         lenold2 = mean(xtrans);
%         
%         rantrans = cosd(theta)*(xran - cmx(nsig(ynb,xnb))) - sind(theta)*(yran - cmy(nsig(ynb,xnb)));
%         xtrans = [xtrans; rantrans];
% %         lennew2 = max(xtrans) - min(xtrans);
%         lennew2 = mean(xtrans);
%         
%         elold2 = (lenold2 - p.len)^2;
%         elnew2 = (lennew2 - p.len)^2;
%     else
%         elold2 = 0;
%         elnew2 = 0;
%     end

% 	eoldtot = enold + p.xlam*etold + epold + p.ecm*eecmold + p.perimstr*eperold + p.lenstr*(elold1+elold2);
% 	enewtot = ennew + p.xlam*etnew + epnew + p.mmt*(emnew1 + emnew2) + p.ecm*eecmnew + p.perimstr*epernew + p.spring*(espnew2 - espnew1) + p.lenstr*(elnew1+elnew2);

	eoldtot = enold + p.xlam*etold + epold + p.ecm*eecmold;
	enewtot = ennew + p.xlam*etnew + epnew + p.mmt*(emnew1 + emnew2) + p.ecm*eecmnew;
    
    
	dE = enewtot - eoldtot;
end