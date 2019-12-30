.dualtable = function(dualvars) {
	tab = matrix(dualvars,  ncol = 2, byrow = TRUE)
	rownames(tab) = names(dualvars)[seq(2, length(dualvars), 2)]
	colnames(tab) = c("Upper", "Lower")
	return(tab)
}