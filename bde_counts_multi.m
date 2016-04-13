function L = bde_counts_multi(pd,nd,plevels,nlevels,S)
	L = 0;
	% Number of parents
	nParents = size(pd,1);
	if nParents == 0
		dataSamp = size(nd,2);
		priorSamp = S/nlevels;
		L = L - gln(S, dataSamp);
		for k=0:nlevels-1
			dataSamp = sum(nd==k);
			L = L + gln(priorSamp, dataSamp);
		end
	else
		% Number of configurations
		nConfig = prod(plevels);
		% Prior sample size
		priorSamp = S/(nlevels*nConfig);
		ndata = size(pd, 2);
		configs = zeros(ndata, 1);
		for i=1:ndata
			configs(i) = multivec2dec(pd(:,i), plevels);
		end
		for j=0:nConfig-1
			configDataRows = (configs==j);
			if any(configDataRows)
				dataSamp = sum(configDataRows);
				L = L - gln(priorSamp*nlevels,dataSamp);
				for k=0:nlevels-1
					dataSamp = sum(nd(configDataRows)==k);
					L = L + gln(priorSamp,dataSamp);
				end
			end
		end
	end
end

function g = gln(s,d)
	g = gammaln(s+d)-gammaln(s);
    %fprintf('G(%.2f)/G(%.2f) \n',s+d,s);
end