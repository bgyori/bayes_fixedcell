function data = read_data_table(cellLine)
	dirname = 'data/';
	erkSignal = load([dirname cellLine '_pERK.mat']);
	aktSignal = load([dirname cellLine '_pAKT.mat']);
	foxoSignal = load([dirname cellLine '_FOXO3a.mat']);
	
	% Get size of data matrix (same for ERK, AKT and FOXO)
	[num_time, num_ligand, num_inh] = size(erkSignal.single_pERK);

	% Remove infinity data points
	for t = 1:num_time
		for lig = 1:num_ligand
			for inh = 1:num_inh
				inf_idx = isinf(erkSignal.single_pERK{t,lig,inh});
				erkSignal.single_pERK{t,lig,inh} = ...
					erkSignal.single_pERK{t,lig,inh}(~inf_idx);
				inf_idx = isinf(aktSignal.single_pAKT{t,lig,inh});
				aktSignal.single_pAKT{t,lig,inh} = ...
					aktSignal.single_pAKT{t,lig,inh}(~inf_idx);
				inf_idx = isinf(foxoSignal.single_FOXO3a{t,lig,inh}(:,1));
				foxoSignal.single_FOXO3a{t,lig,inh} = ...
					foxoSignal.single_FOXO3a{t,lig,inh}(~inf_idx,1);
			end
		end
	end
	
	% Define functions to get median and iqr of a vector
	qm = @(x) quantile(x,0.5);
	iqr = @(x) quantile(x,0.75)-quantile(x,0.25);	
	
	% Calculate median and iqr of single-cell measurements, separately
    % for each ligand and inhibitor condition
    for t = 1:num_time
		for lig = 1:num_ligand
			for inh = 1:num_inh
				erkSignal.qm(t,lig,inh) = qm(erkSignal.single_pERK{t,lig,inh});
				aktSignal.qm(t,lig,inh) = qm(aktSignal.single_pAKT{t,lig,inh});
				foxoSignal.qm(t,lig,inh) = qm(foxoSignal.single_FOXO3a{t,lig,inh});
				foxoSignal.iqr(t,lig,inh) = iqr(foxoSignal.single_FOXO3a{t,lig,inh});
			end
		end
    end
	
	
	if num_time==8
		time_points = [0, 15, 30, 60, 90, 120, 180, 240];
	elseif num_time==13
		time_points = [0, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300, 480];
	end
	%TODO: check if this order is correct
	ligand_names = {'EGF', 'IGF1', 'FGF1', 'HRG', 'HGF', 'EPR', 'BTC', 'NS'};
	
	% Set values of table columns
	time_values = zeros(num_time, num_ligand, num_inh);
	ligand_values = cell(num_time, num_ligand, num_inh);
	akti_values = zeros(num_time, num_ligand, num_inh);
	meki_values = zeros(num_time, num_ligand, num_inh);
	for i=1:num_time
		time_values(i,:,:) = time_points(i);
		for j=1:num_ligand
			for k=1:num_inh
				akti_values(i,j,k) = (k==2 || k==4);
				meki_values(i,j,k) = (k==3 || k==4);
				ligand_values{i,j,k} = ligand_names{j};
			end
		end
	end
	time_col = time_values(:);
	ligand_col = ligand_values(:);
	akti_col = akti_values(:);
	meki_col = meki_values(:);
	
	variable_names = {'ligand', 'akti', 'meki', 'time', ...
		'erk_mu', 'akt_mu', 'foxo_qm', 'foxo_iqr', 'foxo_pd'};
	
    % Use normalization code for FOXO signal
	foxoSignal = normalize_foxo(foxoSignal);
	parabola_params = polyfit(foxoSignal.qm(:), foxoSignal.iqr(:), 2);
	foxoSignal.pd = get_parabola_positions(parabola_params,...
										   foxoSignal.qm(:),...
								           foxoSignal.iqr(:));
	
	data = table(ligand_col, akti_col, meki_col, time_col, ...
		erkSignal.qm(:), aktSignal.qm(:), ...
		foxoSignal.qm(:), foxoSignal.iqr(:), foxoSignal.pd(:), ...
		'VariableNames', variable_names);
	
	% Subtract first time point, inihibited condition, mean across ligands
	data.akt_mu = data.akt_mu - ...
		mean(data.akt_mu(data.time==0 & data.akti==1));
    % Normalize to maximal signal
	data.akt_mu = data.akt_mu ./ abs(max(data.akt_mu));
    % Subtract first time point, inihibited condition, mean across ligands
	data.erk_mu = data.erk_mu - ...
		mean(data.erk_mu(data.time==0 & data.meki==1));
	% Normalize to maximal signal
    data.erk_mu = data.erk_mu ./ abs(max(data.erk_mu));
end