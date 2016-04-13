function p = get_edge_probs(LH)
	% AKT to ERK
	LAE = LH.E(2);
	p.pAE = exp(LAE-logsumexp(LH.E));
	% ERK to AKT
	LEA = LH.A(2);
	p.pEA = exp(LEA-logsumexp(LH.A));
	% ERK to FOXO
	LEF = logsumexp(LH.F([2,4]));
	p.pEF = exp(LEF-logsumexp(LH.F));
	% AKT to FOXO
	LAF = logsumexp(LH.F([3,4]));
	p.pAF = exp(LAF-logsumexp(LH.F));
end