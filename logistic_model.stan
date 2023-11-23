data{
    int<lower = 0> N; //number of observations
    real p_m; //prior mean
    real p_sd; //prior sd
    int conc[N]; //presence/absence of concussion
    int group[N]; //grouping variable
}

parameters{
     vector[2] a; // indexed intercept of concussion for each group
}

model{
     vector[N] p;
    a ~ normal(p_m , p_sd); //variable priors

    for ( i in 1:N ) {
        p[i] = a[group[i]];
        p[i] = inv_logit(p[i]);
    }

    conc ~ binomial( 1 , p );
}

generated quantities{
    vector[N] log_lik;
     vector[N] p;
    for (i in 1:N) {
        p[i] = a[group[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:N ) log_lik[i] = binomial_lpmf( conc[i] | 1 , p[i] );
}
