function d_matrix = plot_inhib_effect_new_mean()
	addpath util
    cell_lines = {'184A1','MCF10A','SKBR3','HCC1806','HS578T', ...
				'MDA231','BT20','MCF7','T47D'};
	for c = 1:length(cell_lines)
        data = read_data_table(cell_lines{c});
        ligands = unique(data.ligand);
		
        % ERK upon AKT inhibition
		akt_inhib = data((data.meki==0) & (data.akti==1),:);
		akt_noinhib = data((data.meki==0) & (data.akti==0),:);
		
		ts = akt_inhib.time(strcmp(akt_inhib.ligand, 'NS'));
		auct = @(x) auc(ts, x);
		
		akt_inhib_mean = varfun(@mean, akt_inhib, ...
			'InputVariables', 'erk_mu', ...
			'GroupingVariables', 'ligand');
        
		akt_noinhib_mean = varfun(@mean, akt_noinhib, ...
			'InputVariables', 'erk_mu', ...
			'GroupingVariables', 'ligand');
        
		akt_inhib_mean.diff = akt_inhib_mean.mean_erk_mu - akt_noinhib_mean.mean_erk_mu;

			
		% AKT upon MEK inhibition
		erk_inhib = data((data.akti==0) & (data.meki==1),:);
		erk_noinhib = data((data.akti==0) & (data.meki==0),:);
        
		erk_inhib_mean = varfun(@mean, erk_inhib, ...
			'InputVariables', 'akt_mu', ...
			'GroupingVariables', 'ligand');
        
		erk_noinhib_mean = varfun(@mean, erk_noinhib, ...
			'InputVariables', 'akt_mu', ...
			'GroupingVariables', 'ligand');
        
		erk_inhib_mean.diff = erk_inhib_mean.mean_akt_mu - erk_noinhib_mean.mean_akt_mu;

		
		% FOXO upon AKT inhibition
		foxo_akt_inhib = data((data.meki==0) & (data.akti==1),:);
		foxo_akt_noinhib = data((data.meki==0) & (data.akti==0),:);
		
        
		foxo_akt_inhib_mean = varfun(@mean, foxo_akt_inhib, ...
			'InputVariables', 'foxo_pd', ...
			'GroupingVariables', 'ligand');
        
		foxo_akt_noinhib_mean = varfun(@mean, foxo_akt_noinhib, ...
			'InputVariables', 'foxo_pd', ...
			'GroupingVariables', 'ligand');
        
		foxo_akt_inhib_mean.diff = foxo_akt_inhib_mean.mean_foxo_pd - foxo_akt_noinhib_mean.mean_foxo_pd;


		% FOXO upon ERK inhibition
		foxo_erk_inhib = data((data.meki==1) & (data.akti==0),:);
		foxo_erk_noinhib = data((data.meki==0) & (data.akti==0),:);
        
		foxo_erk_inhib_mean = varfun(@mean, foxo_erk_inhib, ...
			'InputVariables', 'foxo_pd', ...
			'GroupingVariables', 'ligand');
        
		foxo_erk_noinhib_mean = varfun(@mean, foxo_erk_noinhib, ...
			'InputVariables', 'foxo_pd', ...
			'GroupingVariables', 'ligand');
        
		foxo_erk_inhib_mean.diff = foxo_erk_inhib_mean.mean_foxo_pd - foxo_erk_noinhib_mean.mean_foxo_pd;

		
		for i = 1:length(ligands)
			dE(c,i) = akt_inhib_mean{ligands{i}, {'diff'}};
			dA(c,i) = erk_inhib_mean{ligands{i}, {'diff'}};
			dFA(c,i) = foxo_akt_inhib_mean{ligands{i}, {'diff'}};
			dFE(c,i) = foxo_erk_inhib_mean{ligands{i}, {'diff'}};
		end
	end
    
	clim =max(abs([dA(:);dE(:)]));
	
    figure;
    colormap(get_colormap());
    subplot(2,1,1);
    imagesc(dA, [-clim, clim]);
    set(gca,'xticklabelrotation', 90)
    title('MEKi driven AKT change')
    colorbar;
    set(gca, 'yticklabel', cell_lines);
    set(gca, 'xticklabel', ligands); 
    subplot(2,1,2);
    imagesc(dE, [-clim, clim]);
    set(gca,'xticklabelrotation', 90)
    title('AKTi-driven ERK change');
    colorbar;
    set(gca, 'yticklabel', cell_lines);
    set(gca, 'xticklabel', ligands);

    
	figure;
    hold on;
	
    subplot(2,1,1);
    title('MEKi driven AKT change');
    hold on;
    plot([0,length(cell_lines)+1],[0,0],'k-');
    plot(dA,'r+');
    xlim([0 length(cell_lines)+1]);
	ylim([-clim clim]);
    set(gca,'xtick',1:length(cell_lines));
    set(gca,'xticklabel',cell_lines);
    set(gca,'xticklabelrotation', 45);

    
    subplot(2,1,2);
    title('AKTi driven ERK change');
    hold on;
    plot([0,length(cell_lines)+1],[0,0],'k-');
    plot(dE,'k+');
    xlim([0 length(cell_lines)+1]);
	ylim([-clim clim]);
    set(gca,'xtick',1:length(cell_lines));
    set(gca,'xticklabel',cell_lines);
    set(gca,'xticklabelrotation', 45);
	
	clim = max(abs([dA(:);dE(:)]));

	dE_mean = mean(dE,2);
	dA_mean = mean(dA,2);
	dFA_mean = mean(dFA,2);
	dFE_mean = mean(dFE,2);
	d_matrix = [dE_mean, dA_mean, dFA_mean, dFE_mean];
	figure;
	colormap(get_colormap());
	imagesc(d_matrix', [-clim, clim]);
	set(gca, 'xticklabel', cell_lines);
	set(gca,'xticklabelrotation', 90)
	set(gca,'ytick', [1,2,3,4]);
    set(gca, 'yticklabel', {'AKTi->ERK', 'MEKi->AKT', ...
							'AKTi->FOXO3a', 'MEKi->FOXO3a'});
	colorbar;
end

function cm = get_colormap()
	g = [linspace(0.3,1,32)';flipud(linspace(0.3,1,32)')];
	r = [linspace(0.3,1,32)'; ones(32,1)];
    b = [ones(32,1); flipud(linspace(0.3,1,32)')]; 
    cm = [r, g, b];
end
