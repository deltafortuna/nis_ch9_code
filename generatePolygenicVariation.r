makeStartingDataset <- function(n, freqvec, ofname)
{
	m <- matrix(0, ncol = length(freqvec), nrow = n);
	names <- vector();
	for (i in 1:length(freqvec)) {
		random_vecky <- sample( c( rep(1,freqvec[i]), rep(0, n-freqvec[i]) ) );
		m[,i] <- random_vecky;
		names <- c(names, paste("Locus", i, sep = ""));
	}
	d <- as.data.frame(m);
	colnames(d) <- names;
	write.table(d, file = ofname, row.names = F, col.names = T, quote= F, sep = " ");
}
