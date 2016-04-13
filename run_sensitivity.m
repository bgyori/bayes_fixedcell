rerun = false;
ligand_type = 'none';
if rerun
	ns = 100;
	cell_lines = {'184A1', 'MCF10A', 'SKBR3', 'HCC1806', ...
				'HS578T', 'MDA231', 'BT20', 'MCF7', 'T47D'};
	for c=1:length(cell_lines)
		data{c} = read_data_table(cell_lines{c});
	end
	
	columns = {'cell_line', 'algo', 'pAE', 'pEA', 'pEF', 'pAF'};
	sens = cell2table(cell(0, 6), ...
			'VariableNames', columns);
	for c=1:length(cell_lines)
		fprintf('Sensitivities for %s...\n', cell_lines{c});
		for i=1:ns
			[LDBN_Gauss,LDBN_BGe,LDBN_BDe,LBN_BGe,LBN_BDe] = ...
				learn_structure(data{c}, ligand_type, true);
			pdbde = get_edge_probs(LDBN_BDe);
			pdbga = get_edge_probs(LDBN_Gauss);
			sens = [sens; cell2table({cell_lines{c}, ...
				'DBN_BDe', pdbde.pAE, pdbde.pEA, pdbde.pEF, pdbde.pAF}, ...
				'VariableNames', columns)];
			sens = [sens; cell2table({cell_lines{c}, ...
				'DBN_Gauss', pdbga.pAE, pdbga.pEA, pdbga.pEF, pdbga.pAF}, ...
				'VariableNames', columns)];
		end
	end
	save(['sensitivities_randn_scale0.05_',ligand_type,'.mat'], ...
			'sens', 'data', 'cell_lines')
else
	load(['sensitivities_randn_scale0.05_',ligand_type,'.mat']);
end

plot_edge_probs(sens(strcmp(sens.algo, 'DBN_BDe'),:), cell_lines)
title(['BDe ', ligand_type])
plot_edge_probs(sens(strcmp(sens.algo, 'DBN_Gauss'),:), cell_lines)
title(['Gauss ', ligand_type])

means = varfun(@mean, sens, 'inputvariables', ...
	{'pAE', 'pEA', 'pEF', 'pAF'}, ...
	'groupingvariables', {'cell_line', 'algo'});

figure; hold on;
title('BDe vs Gaussian DBN scores');
xlabel('BDe');
ylabel('Gaussian');
plot(means{strcmp(means.algo, 'DBN_BDe'), 'mean_pAE'}, ...
	means{strcmp(means.algo, 'DBN_Gauss'), 'mean_pAE'}, 'kv', ...
	'markersize',12);

plot(means{strcmp(means.algo, 'DBN_BDe'), 'mean_pEA'}, ...
	means{strcmp(means.algo, 'DBN_Gauss'), 'mean_pEA'}, 'r^', ...
	'markersize',12);

plot(means{strcmp(means.algo, 'DBN_BDe'), 'mean_pAF'}, ...
	means{strcmp(means.algo, 'DBN_Gauss'), 'mean_pAF'}, 'gd', ...
	'markersize',12);

plot(means{strcmp(means.algo, 'DBN_BDe'), 'mean_pEF'}, ...
	means{strcmp(means.algo, 'DBN_Gauss'), 'mean_pEF'}, 'bo', ...
	'markersize',12);

xlim([0, 1])
ylim([0, 1])
axis square

d_matrix = plot_inhib_effect();
figure; hold on;
for i=1:length(cell_lines)
	plot(means{strcmp(means.algo, 'DBN_BDe') & ...
			strcmp(means.cell_line, cell_lines{i}), ...
			'mean_pAE'}, abs(d_matrix(i,1)), 'kv', ...
		'markersize',12);
	x1(i) = means{strcmp(means.algo, 'DBN_BDe') & ...
			strcmp(means.cell_line, cell_lines{i}), ...
			'mean_pAE'};
	y1(i) = abs(d_matrix(i,1));

	plot(means{strcmp(means.algo, 'DBN_BDe') & ...
			strcmp(means.cell_line, cell_lines{i}), ...
			'mean_pEA'}, abs(d_matrix(i,2)), 'r^', ...
		'markersize',12);
	
	x2(i) = means{strcmp(means.algo, 'DBN_BDe') & ...
			strcmp(means.cell_line, cell_lines{i}), ...
			'mean_pEA'};
	y2(i) = abs(d_matrix(i,2));
end
xlim([0, 1]);