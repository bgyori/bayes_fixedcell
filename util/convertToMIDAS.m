cellLines = {'HS578T','184A1','BT20','SKBR3','MDA231','MCF7','HCC1806','MCF10A'};
ligands = {'IGF','EGF','BTC','HGF'};
meki = {'0.1','0'};
akti = {'1','0.25','0.0625','0'};
ts = [0, 15, 45, 90];
species = {'muERK','muAKT','muFOXO','iqrERK','iqrAKT','iqrFOXO'};

dirname = '../BEN - FOXO3a DBN/FixedCellData-TitratingMEKiAKTi/';
fname = 'FOXO_MIDAS.csv';
fh = fopen(fname,'w');
for i = 1:length(cellLines)
	fprintf(fh,'TR:%s,',cellLines{i});
end
for i=1:length(ligands)
	fprintf(fh,'TR:%s,',ligands{i});
end
%for i=1:length(akti)
	%fprintf(fh,'TR:AKTi=%s,',akti{i});
%end
fprintf(fh,'TR:MEKi,');
fprintf(fh,'TR:AKTi,');
%for i=1:length(meki)
%	fprintf(fh,'TR:MEKi=%s,',meki{i});
%end

for i=1:length(species)
	fprintf(fh,'DA:%s,',species{i});
end
for i=1:length(species)
	fprintf(fh,'DV:%s,',species{i});
end
fprintf(fh,'\n');

iqr = @(x) quantile(x,0.75)-quantile(x,0.25);

for c = 1:length(cellLines)
	cellLine = cellLines{c};
	aktSignal = load([dirname cellLine '_pAKT.mat']);
	erkSignal = load([dirname cellLine '_pERK.mat']);
	foxoSignal = load([dirname cellLine '_FOXO3a.mat']);
	
	muaktd = squeeze(cellfun(@median,aktSignal.single_pSignal(:,:,1,:,:)));
	muerkd = squeeze(cellfun(@median,erkSignal.single_pSignal(:,:,1,:,:)));
	mufoxod = squeeze(cellfun(@median,foxoSignal.single_FOXO3aRatio(:,:,1,:,:)));
	iqraktd = squeeze(cellfun(iqr,aktSignal.single_pSignal(:,:,1,:,:)));
	iqrerkd = squeeze(cellfun(iqr,erkSignal.single_pSignal(:,:,1,:,:)));
	iqrfoxod = squeeze(cellfun(iqr,foxoSignal.single_FOXO3aRatio(:,:,1,:,:)));
	
	clstr = '';
	for i = 1:length(cellLines)
		if i==c
			clstr = [clstr '1,'];
		else
			clstr = [clstr ','];
		end
	end
	
	
	% Enumerate each line
	for i=1:length(ligands)
		ligstr = '';
		for j=1:length(ligands)
			if i==j
				ligstr = [ligstr '1,'];
			else
				ligstr = [ligstr ','];
			end
		end
		
		
		for j=1:length(meki)
%			aktistr = '';
% 			for k=1:length(akti)
% 				if j==k
% 					aktistr = [aktistr '1,'];
% 				else
% 					aktistr = [aktistr ','];
% 				end
% 			end
			mekistr = [meki{j}, ','];

			
			for k=1:length(akti)
% 				mekistr = '';
% 				for l=1:length(meki)
% 					if k==l
% 						mekistr = [mekistr '1,'];
% 					else
% 						mekistr = [mekistr ','];
% 					end
% 				end
				aktistr = [akti{k}, ','];
				
				for l=1:length(ts)
					tstr = '';
					for ll=1:length(species)
						tstr = [tstr sprintf('%d,',ts(l))];
					end
					
					if l==1 && (k>1 || j>1 || i>1)
						sigstr = '';
						sigstr = [sigstr sprintf('%.3f,',muerkd(l,1,1,1))];
						sigstr = [sigstr sprintf('%.3f,',muaktd(l,1,1,1))];
						sigstr = [sigstr sprintf('%.3f,',mufoxod(l,1,1,1))];
						sigstr = [sigstr sprintf('%.3f,',iqrerkd(l,1,1,1))];
						sigstr = [sigstr sprintf('%.3f,',iqraktd(l,1,1,1))];
						sigstr = [sigstr sprintf('%.3f',iqrfoxod(l,1,1,1))];
					else
						sigstr = '';
						sigstr = [sigstr sprintf('%.3f,',muerkd(l,i,k,j))];
						sigstr = [sigstr sprintf('%.3f,',muaktd(l,i,k,j))];
						sigstr = [sigstr sprintf('%.3f,',mufoxod(l,i,k,j))];
						sigstr = [sigstr sprintf('%.3f,',iqrerkd(l,i,k,j))];
						sigstr = [sigstr sprintf('%.3f,',iqraktd(l,i,k,j))];
						sigstr = [sigstr sprintf('%.3f',iqrfoxod(l,i,k,j))];
					end
					linestr = [clstr ligstr mekistr aktistr tstr sigstr];
					fprintf(fh,'%s\n',linestr);
				end
			end
		end
	end	
end