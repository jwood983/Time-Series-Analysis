# if you don't have this file, then this won't work :/
mydata <- read.table("Nile.dat", header=T)

estimate_kalman <- function(y, const, q) {
	n <- length(y);
	# create vectors
	vt <- rep(0, n);
	Ft <- rep(0, n); Ft[1] <- NA;
	at <- rep(0, n);
	Pt <- rep(0, n);
	
	Pt[2] <- (1 + q);
	at[2] <- y[1];
	for (t in 2:(n-1)) {
		vt[t] <- (y[t] - at[t]);
		Ft[t] <- (Pt[t] + 1);
		K <- (Pt[t] / Ft[t]);
		at[t+1] <- (at[t] + K*vt[t]);
		Pt[t+1] <- (Pt[t] * (1 - K) + q);
	}
	vt[n] <- (y[n] - at[n]);
	Ft[n] <- (Pt[n] + 1);
	# sigma hat
	sighat <- 0;
	Fhat <- 0;
	for (t in 2:n) {
		sighat <- (sighat + vt[t]*vt[t] / Ft[t]);
		Fhat <- (Fhat + log(Ft[t]));
	}
	sighat <- (sighat / (n-1) );
	return(const - (n-1)*log(sighat)/2 - 0.5*Fhat);
}

# estimate the sigmas
estimate_sigma <- function(y) {
	# length of the input array
	n <- length(y);
	# the constant
	C1 <- -n*log(2*pi)/2.0 - (n-1)/2.0;
	# initialize stored variables
	score <- 200;
	iter <- 0;
	psi <- 0;
	while (abs(score) >= 1e-10) {
		q <- exp(psi);
		this_log <- estimate_kalman(y, C1, q);
		score <- ((estimate_kalman(y, C1, exp(psi+0.0001)) - estimate_kalman(y, C1, exp(psi-0.0001))) / 0.0002);
		cat(iter, q, psi, score, this_log, "\n", sep="\t");
		if (iter == 0) {
			psi <- score;
		} else {
			psi <- (score + psi);
		}
		iter <- (iter + 1);
	}
	# get sigma_e^2
	vt <- rep(0, n);
	Ft <- rep(0, n); Ft[1] <- NA;
	at <- rep(0, n);
	Pt <- rep(0, n);
	
	q <- exp(psi);
	Pt[2] <- (1 + q);
	at[2] <- y[1];
	for (t in 2:(n-1)) {
		vt[t] <- (y[t] - at[t]);
		Ft[t] <- (Pt[t] + 1);
		K <- (Pt[t] / Ft[t]);
		at[t+1] <- (at[t] + K*vt[t]);
		Pt[t+1] <- (Pt[t] * (1 - K) + q);
	}
	vt[n] <- (y[n] - at[n]);
	Ft[n] <- (Pt[n] + 1);
	# sigma hat
	sigma_e2 <- 0;
	for (t in 2:n) {
		sigma_e2 <- (sigma_e2 + vt[t]*vt[t] / Ft[t]);
	}
	sigma_e2 <- (sigma_e2 / (n-1) );
	sigma_n2 <- q * sigma_e2;
	cat("sigma_e^2 = ", sigma_e2, "; sigma_n^2 = ", sigma_n2, "\n", sep="");
	return (list(v1=sigma_e2, v2=sigma_n2) );
}

# filter the data
kalman <- function(y, t, sig_e2, sig_n2, P0, a0) {
	n <- length(y);
	# create vector for Kalman filter
	at <- rep(0, n);
	vt <- rep(0, n);
	Ft <- rep(0, n);
	Kt <- rep(0, n);
	at <- rep(0, n);
	Pt <- rep(0, n);
	
	Pt[1] <- P0;
	at[1] <- a0;
	
	for (i in 1:(n-1)) {
		vt[i] <- (y[i] - at[i]);
		Ft[i] <- (Pt[i] + sig_e2);
		Kt[i] <- (Pt[i] / Ft[i]);
		at[i+1] <- (at[i] + Kt[i] * vt[i]);
		Pt[i+1] <- (Pt[i] * (1 - Kt[i]) + sig_n2);
	}
	vt[n] <- (y[n] - at[n]);
	Ft[n] <- (P[n] + sig_e2);
	Kt[n] <- (Pt[n] / Ft[n]);
	
	return(list(v1=vt,v2=Ft,v3=Kt,v4=at,v5=Pt))
}

# smooth the state
smooth_state <- function(filter) {
	# get the vectors
	vt <- filter$v1;
	Ft <- filter$v2;
	Kt <- filter$v3;
	at <- filter$v4;
	Pt <- filter$v5;
	# get lengths
	n <- length(vt);
	
	# given the filter, get the smoothed state
	rt <- rep(0, n);
	ah <- rep(0, n);
	
	for (t in n:2) {
		rt[t-1] <- (vt[t] / Ft[t] + (1 - Kt[t]) * rt[t]);
		ah[t] <- (at[t] + Pt[t] * rt[t-1]);
	}
	rt0 <- (vt[1] / Ft[1] + (1 - Kt[1]) * rt[1]);
	ah[1] <- (at[1] + Pt[1] * rt0);
	return(list(v1=rt, v2=ah))
}

# get the smoothed variance
smooth_variance <- function(filter) {
	# get the vectors
	vt <- filter$v1;
	Ft <- filter$v2;
	Kt <- filter$v3;
	at <- filter$v4;
	Pt <- filter$v5;
	# get lengths
	n <- length(vt);
	
	# given the filter, get the smoothed variances
	Nt <- rep(0, n);
	Vh <- rep(0, n);
	
	for (t in n:2) {
		Nt[t-1] <- (1/Ft[t] + (1 - Kt[t])^2 * Nt[t]);
		Vh[t] <- (Pt[t] * (1 - Pt[t] * Nt[t-1]));
	}
	Vh[1] <- (Pt[1] * (1 - Pt[1] * ((1/Ft[1] + (1 - Kt[1])^2 * Nt[1]))));
	return(list(v1=Nt,v2=Vh));
}

# predict from Kalman
predict <- function (filter, sig_n2) {
	# get the vectors
	vt <- filter$v1;
	Ft <- filter$v2;
	Kt <- filter$v3;
	at <- filter$v4;
	Pt <- filter$v5;
	# get lengths
	n <- length(vt);
	
	# copy at to anew
	anew <- rep(0,(n+30));
	Pnew <- rep(0,(n+30));
	for (i in 1:n) {
		anew[i] <- at[i];
		Pnew[i] <- Pt[i];
	}
	
	# predict by assuming n+1 -> n+30 is missing
	for (t in (n+1):(n+30)) {
		anew[t] <- anew[t-1];
		Pnew[t] <- (Pnew[n] + sig_n2);
	}
	return(list(v1=anew,v2=Pnew))
}

# given the filter, return the errors
error_moments <- function(filter) {
	# get the vectors
	vt <- filter$v1;
	Ft <- filter$v2;
	Kt <- filter$v3;
	at <- filter$v4;
	Pt <- filter$v5;
	# get length
	n <- length(vt);
	# get the error term
	e1 <- rep(0, n);
	# the moments
	m1 <- 0;
	m2 <- 0;
	m3 <- 0;
	m4 <- 0;
	for (t in 1:n) {
		e1[t] <- (vt[t] / sqrt(Ft[t]));
		m1 <- (m1 + e1[t]);
	}
	m1 <- (m1/n);
	for (t in 1:n) {
		m2 <- (m2 + (e1[t] - m1)^2);
		m3 <- (m3 + (e1[t] - m1)^3);
		m4 <- (m4 + (e1[t] - m1)^4);
	}
	m2 <- (m2/2);
	m3 <- (m3/2);
	m4 <- (m4/2);
	return(list(v0=e1, v1=m1, v2=m2, v3=m3, v4=m4));
}

# get sigma's
sigma <- estimate_sigma(mydata$flow.rate);
# get the kalman filter
kalman_filter <- kalman(mydata$flow.rate, mydata$date, sigma$v1, sigma$v2, 1e7, 0);
# get the smoothed state
smoothed_state <- smooth_state(kalman_filter);
# get the smoothed variances
smoothed_variance <- smooth_variance(kalman_filter);
# predict model
pred <- predict(kalman_filter, 1469.1);
