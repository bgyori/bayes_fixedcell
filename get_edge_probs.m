function p = get_edge_probs(LH)
	% AKT to ERK
	LAE = LH.E(2);
	p.pAE = exp(LAE-logsumexp(LH.E));
	% ERK to AKT
	LEA = LH.A(2);
	p.pEA = exp(LEA-logsumexp(LH.A));
	% ERK to FOXO
	LF = LH.Fm + LH.Fi;
	LEF = logsumexp(LF([2,4]));
	p.pEF = exp(LEF-logsumexp(LF));
	% AKT to FOXO
	LAF = logsumexp(LF([3,4]));
	p.pAF = exp(LAF-logsumexp(LF));
end