function learnStructure(data_version,use_repeats)
	%clear variables
	addpath ..
	addpath util
	if nargin < 2
		use_repeats = 0;
		if nargin < 1
			data_version = 2;
		end
	end
	fprintf('Running data version %d\n',data_version);
	if data_version==2
		fprintf('Running repeat %d\n',use_repeats);
	end

	qm = @(x) quantile(x,0.5);
	iqr = @(x) quantile(x,0.75)-quantile(x,0.25);

	if data_version==1
		cellLines = {'HS578T','184A1','BT20','SKBR3','MDA231','MCF7','HCC1806',...
				'T47D','MCF10A10142013','MCF10A11292013','MCF10A12202013'};
	else
		cellLines = {'HS578T','184A1','BT20','SKBR3','MDA231','MCF7','HCC1806','MCF10A'};
	end

	for cl=1:length(cellLines)
		cellLine = cellLines{cl};

		switch data_version
		case 1
			dirname = 'C:\Users\SOMPONNAT\Dropbox (Somponnat workspace)\Collaborative Projects\BEN - FOXO3a DBN\FixedCellData\';
			erkSignal = load([dirname cellLine '_pERK.mat']);
			aktSignal = load([dirname cellLine '_pAKT.mat']);
			foxoSignal = load([dirname cellLine '_FOXO3a.mat']);
			% Pre-filtering
			s = size(erkSignal.single_pERK);
			if s(1)==8
				ts = [0, 15, 30, 60, 90, 120, 180, 240];
				tidx = 2:8;
			elseif s(1)==13
				ts = [0, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300, 480];
				tidx = 2:13;
			elseif s(1)==6
				ts = [0, 10 20, 30, 60, 120];
				tidx = 2:6;
			end
			
			fprintf('%s:%d\n',cellLine,s(1));

			s = size(erkSignal.single_pERK);

			% Remove infinity
			for t = 1:s(1)
				for lig = 1:s(2)
					for inh = 1:s(3)
						erkSignal.single_pERK{t,lig,inh} = erkSignal.single_pERK{t,lig,inh}(~isinf(erkSignal.single_pERK{t,lig,inh}));
						aktSignal.single_pAKT{t,lig,inh} = aktSignal.single_pAKT{t,lig,inh}(~isinf(aktSignal.single_pAKT{t,lig,inh}));
						foxoSignal.single_FOXO3a{t,lig,inh} = foxoSignal.single_FOXO3a{t,lig,inh}(~isinf(foxoSignal.single_FOXO3a{t,lig,inh}(:,1)),1);
					end
				end
			end

			% Calculate median and iqr
			for t = 1:s(1)
				for lig = 1:s(2)
					for inh = 1:s(3)
						erkSignal.qm(t,lig,inh) = qm(erkSignal.single_pERK{t,lig,inh});
						erkSignal.iqr(t,lig,inh) = iqr(erkSignal.single_pERK{t,lig,inh});
						aktSignal.qm(t,lig,inh) = qm(aktSignal.single_pAKT{t,lig,inh});
						aktSignal.iqr(t,lig,inh) = iqr(aktSignal.single_pAKT{t,lig,inh});
						foxoSignal.qm(t,lig,inh) = qm(foxoSignal.single_FOXO3a{t,lig,inh});
						foxoSignal.iqr(t,lig,inh) = iqr(foxoSignal.single_FOXO3a{t,lig,inh});
					end
				end
			end
			
			% Subtract minimum
			erkSignal.qm = normalizeMinimum(erkSignal.qm);
			erkSignal.iqr = normalizeMinimum(erkSignal.iqr);
			aktSignal.qm = normalizeMinimum(aktSignal.qm);
			aktSignal.iqr = normalizeMinimum(aktSignal.iqr);
			foxoSignal.qm = normalizeMinimum(foxoSignal.qm);
			foxoSignal.iqr = normalizeMinimum(foxoSignal.iqr);
			
			% Interpolate time points
			ts_new = [0, 15, 45, 90]';
			for lig = 1:s(2)
				for inh = 1:s(3)
					erkSignal.qm(1:length(ts_new),lig,inh) = smoothkernel(ts(:),erkSignal.qm(:,lig,inh),ts_new);
					erkSignal.iqr(1:length(ts_new),lig,inh) = smoothkernel(ts(:),erkSignal.iqr(:,lig,inh),ts_new);
					aktSignal.qm(1:length(ts_new),lig,inh) = smoothkernel(ts(:),aktSignal.qm(:,lig,inh),ts_new);
					aktSignal.iqr(1:length(ts_new),lig,inh) = smoothkernel(ts(:),aktSignal.iqr(:,lig,inh),ts_new);
					foxoSignal.qm(1:length(ts_new),lig,inh) = smoothkernel(ts(:),foxoSignal.qm(:,lig,inh),ts_new);
					foxoSignal.iqr(1:length(ts_new),lig,inh) = smoothkernel(ts(:),foxoSignal.iqr(:,lig,inh),ts_new);
				end
			end
			
			
			%HCC1806: IGF1, EGF, HRG, BTC
			%184A1: IGF1, EGF, HRG, BTC
			%MCF10A: IGF1, EGF, HRG, BTC
			%BT20: IGF1, EGF, BTC, HGF
			%MDA231: IGF1, EGF, BTC, HGF
			%HS578T: IGF1, EGF, HRG, EPR
			%MCF7: IGF1, EGF, HRG, EPR
			%SKBR3: IGF1, EGF, HRG, EPR
			ligands = {'EGF','IGF1','FGF1','HRG','HGF','EPR','BTC','NS'};
			switch cellLine
				case 'HCC1806'
					ligidx = [1,2,4,7];
				case '184A1'
					ligidx = [1,2,4,7];
				case 'MCF10A'
					ligidx = [1,2,4,7];
				case 'BT20'
					ligidx = [1,2,5,7];
				case 'MDA231'
					ligidx = [1,2,5,7];
				case 'HS578T'
					ligidx = [1,2,4,6];
				case 'MCF7'
					ligidx = [1,2,4,6];
				case 'SKBR3'
					ligidx = [1,2,4,6];
			end
				
			% Use only the chosen conditions
			tidx = 2:4;
			ts = ts(tidx); 
			%ligidx = [1,2,3,4,5,6,7];
			inhidx = 1:4;
			erkSignal.qm = erkSignal.qm(tidx,ligidx,inhidx);
			aktSignal.qm = aktSignal.qm(tidx,ligidx,inhidx);
			foxoSignal.qm = foxoSignal.qm(tidx,ligidx,inhidx);
			erkSignal.iqr = erkSignal.iqr(tidx,ligidx,inhidx);
			aktSignal.iqr = aktSignal.iqr(tidx,ligidx,inhidx);
			foxoSignal.iqr = foxoSignal.iqr(tidx,ligidx,inhidx);	
			
			s = size(erkSignal.qm);
			
			% Set
			erkSignal.val = ones(s(2:end));
			erkSignal.val(:,[3,4]) = 0;
			aktSignal.val = ones(s(2:end));
			aktSignal.val(:,[2,4]) = 0;
			foxoSignal.val = ones(s(2:end));
		case 2
			dirname = 'C:\Users\SOMPONNAT\Dropbox (Somponnat workspace)\Collaborative Projects\BEN - FOXO3a DBN\FixedCellData-TitratingMEKiAKTi\';
			erkSignal = load([dirname cellLine '_pERK.mat']);
			aktSignal = load([dirname cellLine '_pAKT.mat']);
			foxoSignal = load([dirname cellLine '_FOXO3a.mat']);
			
			
			% Pre-filtering
			s = size(erkSignal.single_pSignal);
			if s(1)==4
				ts = [0, 15, 45, 90];
				tidx = [1 2 3 4]; % 0, 15, 45, 90
			else
				ts = [0, 15, 45, 90, 135];
				tidx = [1 2 3 4]; % 0, 15, 45, 90, 135
			end

			ts = ts(tidx); 
			ligidx = [1 2 3 4];
			inhAidx = [4 3 2 1]; % 0, 0.0625, 0.25, 1
			inhMidx = [2 1]; % 0, 0.1

			% Fill in initial time point for all conditions
			for lig = 2:s(2)
				for rep = 1:s(3)
					for inhA = 1:s(4)
						for inhM = 1:s(5)
							erkSignal.single_pSignal(1,lig,rep,inhA,inhM) = ...
								erkSignal.single_pSignal(1,1,rep,inhA,inhM);
							aktSignal.single_pSignal(1,lig,rep,inhA,inhM) = ...
								aktSignal.single_pSignal(1,1,rep,inhA,inhM);
							foxoSignal.single_FOXO3aRatio(1,lig,rep,inhA,inhM) = ...
								foxoSignal.single_FOXO3aRatio(1,1,rep,inhA,inhM);
						end
					end
				end
			end
			
			
			erkSignal.qm = cellfun(qm,erkSignal.single_pSignal(tidx,ligidx,:,inhAidx,inhMidx));
			erkSignal.iqr = cellfun(iqr,erkSignal.single_pSignal(tidx,ligidx,:,inhAidx,inhMidx));
			aktSignal.qm = cellfun(qm,aktSignal.single_pSignal(tidx,ligidx,:,inhAidx,inhMidx));
			aktSignal.iqr = cellfun(iqr,aktSignal.single_pSignal(tidx,ligidx,:,inhAidx,inhMidx));
			foxoSignal.qm = cellfun(qm,foxoSignal.single_FOXO3aRatio(tidx,ligidx,:,inhAidx,inhMidx));
			foxoSignal.iqr = cellfun(iqr,foxoSignal.single_FOXO3aRatio(tidx,ligidx,:,inhAidx,inhMidx));
		
			switch use_repeats
			case 1
				erkSignal.qm = squeeze(erkSignal.qm(:,:,1,:,:));
				aktSignal.qm = squeeze(aktSignal.qm(:,:,1,:,:));
				foxoSignal.qm = squeeze(foxoSignal.qm(:,:,1,:,:));
				erkSignal.iqr = squeeze(erkSignal.iqr(:,:,1,:,:));
				aktSignal.iqr = squeeze(aktSignal.iqr(:,:,1,:,:));
				foxoSignal.iqr = squeeze(foxoSignal.iqr(:,:,1,:,:));
			case 2
				erkSignal.qm = squeeze(erkSignal.qm(:,:,2,:,:));
				aktSignal.qm = squeeze(aktSignal.qm(:,:,2,:,:));
				foxoSignal.qm = squeeze(foxoSignal.qm(:,:,2,:,:));
				erkSignal.iqr = squeeze(erkSignal.iqr(:,:,2,:,:));
				aktSignal.iqr = squeeze(aktSignal.iqr(:,:,2,:,:));
				foxoSignal.iqr = squeeze(foxoSignal.iqr(:,:,2,:,:));
			otherwise
				erkSignal.qm = squeeze(mean(erkSignal.qm,3));
				aktSignal.qm = squeeze(mean(aktSignal.qm,3));
				foxoSignal.qm = squeeze(mean(foxoSignal.qm,3));
				erkSignal.iqr = squeeze(mean(erkSignal.iqr,3));
				aktSignal.iqr = squeeze(mean(aktSignal.iqr,3));
				foxoSignal.iqr = squeeze(mean(foxoSignal.iqr,3));
			end
			s = size(erkSignal.qm);
			erkSignal.val = ones(s(2:end));
			erkSignal.val(:,:,2) = 0;
			aktSignal.val = ones(s(2:end));
			aktSignal.val(:,2:4,:) = 0;	
			foxoSignal.val = ones(s(2:end));
		end
		
		

		%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DBN learning
		%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% Temporal data matrix
		nnodes = 6;
		D = zeros(nnodes,s(1),prod(s(2:end)));
		V = zeros(nnodes,prod(s(2:end)));
		D(1,:,:) = reshape(erkSignal.qm,[s(1) prod(s(2:end))]);
		D(2,:,:) = reshape(erkSignal.iqr,[s(1) prod(s(2:end))]);
		D(3,:,:) = reshape(aktSignal.qm,[s(1) prod(s(2:end))]);
		D(4,:,:) = reshape(aktSignal.iqr,[s(1) prod(s(2:end))]);
		D(5,:,:) = reshape(foxoSignal.qm,[s(1) prod(s(2:end))]);
		D(6,:,:) = reshape(foxoSignal.iqr,[s(1) prod(s(2:end))]);
		V(1,:) = reshape(erkSignal.val,[1 prod(s(2:end))]);
		V(2,:) = V(1,:);
		V(3,:) = reshape(aktSignal.val,[1 prod(s(2:end))]);
		V(4,:) = V(3,:);
		V(5,:) = reshape(foxoSignal.val,[1 prod(s(2:end))]);
		V(6,:) = V(5,:);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% BN learning
		%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% Use only the chosen conditions
% 		tidx = 2:4;
% 		ts = ts(tidx); 
% 		ligidx = 1:7;
% 		inhidx = 1:4;
% 		erkSignal.qm = erkSignal.qm(tidx,ligidx,inhidx);
% 		aktSignal.qm = aktSignal.qm(tidx,ligidx,inhidx);
% 		foxoSignal.qm = foxoSignal.qm(tidx,ligidx,inhidx);
% 		erkSignal.iqr = erkSignal.iqr(tidx,ligidx,inhidx);
% 		aktSignal.iqr = aktSignal.iqr(tidx,ligidx,inhidx);
% 		foxoSignal.iqr = foxoSignal.iqr(tidx,ligidx,inhidx);			
		
		% Mean
		xERK = mean(reshape(erkSignal.qm,[s(1) prod(s(2:end))]));
		xERKiqr = mean(reshape(erkSignal.iqr,[s(1) prod(s(2:end))]));
		xAKT = mean(reshape(aktSignal.qm,[s(1) prod(s(2:end))]));
		xAKTiqr = mean(reshape(aktSignal.iqr,[s(1) prod(s(2:end))]));
		xFOXO3a = mean(reshape(foxoSignal.qm,[s(1) prod(s(2:end))]));
		xFOXO3aiqr = mean(reshape(foxoSignal.iqr,[s(1) prod(s(2:end))]));

		Dx(:,1) = xERK;
		Dx(:,2) = xERKiqr;
		Dx(:,3) = xAKT;
		Dx(:,4) = xAKTiqr;
		Dx(:,5) = xFOXO3a;
		Dx(:,6) = xFOXO3aiqr;

		bERK = (xERK>=otsu(xERK));
		bERKiqr = (xERKiqr>=otsu(xERKiqr));
		bAKT = (xAKT>=otsu(xAKT));
		bAKTiqr = (xAKTiqr>=otsu(xAKTiqr));
		bFOXO3a = (xFOXO3a>=otsu(xFOXO3a));
		bFOXO3aiqr = (xFOXO3aiqr>=otsu(xFOXO3aiqr));


		Db(:,1) = bERK;
		Db(:,2) = bERKiqr;
		Db(:,3) = bAKT;
		Db(:,4) = bAKTiqr;
		Db(:,5) = bFOXO3a;
		Db(:,6) = bFOXO3aiqr;

		S  = 1;
		use_iqr = 1;
		
		for af = [0 1]
			for ea = [0 1]
				for ef = [0 1]
					A = zeros(nnodes,nnodes);
					% ERK to AKT
					A(1,3) = ea;
					A(1,4) = ea & use_iqr;
					% ERK to FOXO
					A(1,5) = ef;
					A(1,6) = ef & use_iqr;
					% AKT to FOXO
					A(3,5) = af;
					A(3,6) = af & use_iqr;
					mLH = getMarginalLikelihoodDBN(A,D,V);
					mLH2 = getMarginalLikelihoodDBGe(A,D,V);
					mLH3 = getMarginalLikelihoodBN(A,S,Db,V');
					mLH4 = getMarginalLikelihoodBGe(A,Dx,V);
					LDBNx(cl,4*af+2*ea+ef+1) = mLH(3)+mLH(4)+mLH(5)+mLH(6);
					LDBN_BGex(cl,4*af+2*ea+ef+1) = mLH2(3)+mLH2(4)+mLH2(5)+mLH2(6);
					LBNx(cl,4*af+2*ea+ef+1) = mLH3(3)+mLH3(4)+mLH3(5)+mLH3(6);
					LBN_BGex(cl,4*af+2*ea+ef+1) = mLH4(3)+mLH4(4)+mLH4(5)+mLH4(6);
				end
			end
		end
		
		paf = exp(logsumexpv(LDBNx(cl,[5,6,7,8]))-logsumexpv(LDBNx(cl,:)));
		pea = exp(logsumexpv(LDBNx(cl,[3,4,7,8]))-logsumexpv(LDBNx(cl,:)));
		pef = exp(logsumexpv(LDBNx(cl,[2,4,6,8]))-logsumexpv(LDBNx(cl,:)));
		%figure;
		%title(cellLine)
		%draw_static_graph(paf,pea,pef);
		
		% Choose case where A->F is fixed
		LDBN = LDBNx(:,5:8);
		LDBN_BGe = LDBN_BGex(:,5:8);
		LBN = LBNx(:,5:8);
		LBN_BGe = LBN_BGex(:,5:8);
		
		% Normalize to sum to 1
		LDBN(cl,:) = LDBN(cl,:)-logsumexpv(LDBN(cl,:));
		LDBN_BGe(cl,:) = LDBN_BGe(cl,:)-logsumexpv(LDBN_BGe(cl,:));
		LBN(cl,:) = LBN(cl,:)-logsumexpv(LBN(cl,:));
		LBN_BGe(cl,:) = LBN_BGe(cl,:)-logsumexpv(LBN_BGe(cl,:));
	end

	% Compare AKT->FOXO3a only versus ERK also playing a role
	for cl = 1:length(cellLines)
		LDBN_AvsE(cl) = LDBN(cl,1)-logsumexpv(LDBN(cl,2:end));
		LDBN_BGe_AvsE(cl) = LDBN_BGe(cl,1)-logsumexpv(LDBN_BGe(cl,2:end));
		LBN_AvsE(cl) = LBN(cl,1)-logsumexpv(LBN(cl,2:end));
		LBN_BGe_AvsE(cl) = LBN_BGe(cl,1)-logsumexpv(LBN_BGe(cl,2:end));
	end

	figure;
    [x,idx] = sort(LDBN_AvsE,'descend');
    bar(x);
    set(gca,'xticklabel',cellLines(idx))
	lims = [min(LDBN_AvsE),max(LDBN_AvsE)];
	ylim(lims);
	ylabel('Log relative probability');
	set(gca,'ytick',lims);
	set(gca,'yticklabel',{'ERK and AKT','AKT only'});
	title('DBN learning result');

	trianglePlot(LDBN,cellLines);
	title(sprintf('DBN/Dataset %d/Repeat %d',data_version,use_repeats));
end

function signalNorm = normalizeMinimum(signal)
	s = signal(:);
	signalMin = min(s);
	signalNorm = signal-signalMin;
end

function draw_static_graph(af,ea,ef)
	r = 0.1;
	x = [0 0.5 0.25];
	y = [0.5 0.5 0];
	names = {'AKT','ERK','FOXO'};
	for i=1:length(x)
		draw_node(x(i),y(i),r,names{i});
	end
	axis off
	axis square;
	draw_edge(x([1,3]),y([1,3]),r,af);
	draw_edge(x([2,1]),y([2,1]),r,ea);
	draw_edge(x([2,3]),y([2,3]),r,ef);
end

function draw_node(x,y,r,t)
	rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'FaceColor','w');
	text(x,y,t);
end

function draw_edge(x,y,r,s)
	edge_dir = [x(2)-x(1),y(2)-y(1)]/norm([x(2)-x(1),y(2)-y(1)]);
	edge_start = [x(1),y(1)]+edge_dir*r;
	edge_end = [x(2),y(2)]-edge_dir*r;
	start_pos = d2n(edge_start);
	end_pos = d2n(edge_end);
	annotation('arrow',[start_pos(1),end_pos(1)],[start_pos(2),end_pos(2)],'LineWidth',10*s);
	edge_mid = (edge_end+edge_start)/2;
	text(edge_mid(1),edge_mid(2),sprintf('%.2f',s));
end

function xn = d2n(x)
	xl = xlim();
	yl = ylim();
	xn(1) = (x(1)-xl(1))/(xl(2)-xl(1));
	xn(2) = (x(2)-yl(1))/(yl(2)-yl(1));
end
