functions {
  // Gaussian temperature dependence function
  real T_dep(real sbt, real Topt, real width, real exp_yn) {
    if (exp_yn == 1) {
      return exp(-0.5 * ((sbt - Topt) / width) ^ 2); // gaussian temperature-dependent function; used for T-dependent recruitment and mortality?
    } else {
      // return (-0.1 * ((sbt - Topt)/width)^2); // gaussian temperature-dependent function, not exponentiated, for temperature-dependent process that are exponentiated later (like movement)
      // I'm taking this out because I can't get it to be non-negative without the exp! 
      
      // return normal_lpdf(sbt | Topt, width);
      return log((1 / sqrt(2 * pi() * width))
      * exp(-square(sbt - Topt) / (2 * square(width))));
    }
  }
  
  vector colSums(matrix M) {
    int ncol;
    vector[cols(M)] sums;
    
    ncol = cols(M);
    
    for (i in 1 : ncol) {
      sums[i] = sum(M[ : , i]); //sums[i] = sum(col(M,i)); 
    }
    return sums;
  }
  
  // cumulative sum function
  // slightly different than the built-in Stan one -- this takes a 1D array of numbers and returns a single real number that sums elements 1:x of that array
  real csum(array[] real to_sum, int x) {
    return sum(to_sum[1 : x]);
  }
  
  // functions to calculate range quantiles
  // https://discourse.mc-stan.org/t/how-to-find-the-location-of-a-value-in-a-vector/19768/2
  
  // first: a function that calculates how many matches there are in vector x to real y (need this for counting along the abundance vector to estimate quantiles)
  int num_matches(vector x, real y) {
    int n = 0;
    for (i in 1 : rows(x)) {
      if (x[i] == y) {
        n += 1;
      }
    }
    return n;
  }
  
  // next: a function that returns the positions of those matches along vector x 
  array[] int which_equal(vector x, real y) {
    array[num_matches(x, y)] int match_positions;
    int pos = 1;
    for (i in 1 : num_elements(x)) {
      // example code used size(), got error "no matches for available argument signatures for size"; but beware of applying this function to >1D arrays 
      if (x[i] == y) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }
  
  // finally: function to calculate range quantiles
  
  real calculate_range_quantile(int np, vector patches,
  array[] real dens_by_patch, real quantile_out) {
    vector[np] csum_dens;
    vector[np] csum_to_edge;
    real quant_position;
    real cutoff;
    int cutoff_id;
    
    cutoff = sum(dens_by_patch[1 : np]) * quantile_out; // where along the vector of counts does the quantile cutoff fall? (then need to translate that into a patch position, below)  
    
    for (i in 1 : np) {
      csum_dens[i] = csum(dens_by_patch[1 : np], i); // calculate cumulative sum of density along each patch 
      if (csum_dens[i] <= cutoff) {
        csum_to_edge[i] = csum_dens[i]; // keep only the patches below the edge 
      } else {
        csum_to_edge[i] = 0;
      }
    }
    
    cutoff_id = min(which_equal(csum_to_edge, 0)); // get lowest patch with a 0 (that's the patch where the edge falls)
    
    quant_position = ((cutoff - max(csum_to_edge)) / dens_by_patch[cutoff_id])
    + cutoff_id - 1; // calculate what proportion of the edge-containing patch is "filled in" by the actual fish up to the weighted quantile (that's the decimal) and then add in the other patches 
    
    return quant_position;
  } // close quantile function
  
  /**
  * Run the population model forward in time given temperature etc. Return numbers at age by patch for each year
  * Go through and fill in this stuff would be good. 
  * @param np number of patches
  * @param ny_train number of years of training data
  * @return an array of numbers by age, patch, and year
  */
  
  array[] matrix simulate_population(int np, int ny_train, int n_ages,
  int n_lbins, int age_at_maturity,
  matrix sbt, real Topt, real width,
  real exp_yn, int T_dep_mortality,
  int T_dep_movement, matrix f_a_y,
  real m, real d, real beta_t,
  int T_dep_recruitment,
  int spawner_recruit_relationship,
  vector init_dep,
 // real beta_rec, 
  real sigma_r, vector raw,
  real r0, vector maturity_at_age,
  vector wt_at_age, real alpha, real h,
  real ssb0, vector d_at_age,
  matrix l_at_a_key,
  vector selectivity_at_bin,
  real beta_obs_int, real beta_obs,
  int number_quantiles,
  matrix init_n_at_age,
  int use_init_n_at_age,
  vector capacity) {
    //// define variables ////
    
    matrix[np, ny_train] T_adjust; // tuning parameter for sbt suitability in each patch*year
    
    array[np, n_ages, ny_train] real surv;
    
    matrix[np, np] diff_m = rep_matrix(0, np, np);
    
    matrix[np, np] adj_m; // adjacency matrix for patches 
    
    matrix[np, np] outer; // outer difference matrix 
    
    array[ny_train] matrix[np, np] T_adjust_m; // array containing temperature preference matrices
    
    array[ny_train] matrix[np, np] tax_m; // array containing taxis between patches every year 
    
    array[ny_train] matrix[np, np] mov_inst_m; // array containing instantaneous movement matrices (dispersal + temperature-dependent taxis)
    
    array[ny_train] matrix[np, np] mov_m; // Dan thinks this is the wrong syntax -- need to fix 
    
    array[ny_train] matrix[np, n_ages] n_at_age_hat;
    
    matrix[np, ny_train] ssb;
    
    vector[ny_train - 1] rec_dev; // recruitment deviates in each year (shared across all patches, for better or worse)
    
    vector[np] v_in; // vector for matrix multiplication 
    
    vector[np] v_out; // vector for matrix multiplication 
    
    array[ny_train] matrix[np, n_lbins] n_at_length_hat; // array number of years containing matrices with numbers at patch, length bin, and year 
    
    array[np, ny_train] real density_hat; // for tracking sum density 
    
    matrix[np, ny_train] theta; // Bernoulli probability of encounter  
    
    ssb = rep_matrix(0, np, ny_train);
    
    // in R this is just outer(np, np, "-") 
    for (i in 1 : np) {
      for (j in 1 : np) {
        outer[i, j] = i - j; // should just create a matrix that counts from 1 to np along rows and columns
        if (abs(outer[i, j]) == 1) {
          adj_m[i, j] = 1; // fill with 1s for adjacent patches and 0s otherwise 
        } else {
          adj_m[i, j] = 0;
        }
      }
    }
    
    //// calculate temperature based adjustments ////
    
    for (p in 1 : np) {
      for (y in 1 : cols(sbt)) {
        T_adjust[p, y] = T_dep(sbt[p, y], Topt, width, exp_yn);
      } // close years
    } // close patches
    
    //// calculate total annual mortality from instantaneous natural + fishing mortality data ////
    
    for (p in 1 : np) {
      for (a in 1 : n_ages) {
        for (y in 1 : ny_train) {
          if (T_dep_mortality == 1) {
            surv[p, a, y] = exp(-((f_a_y[a, y] + m))) * T_adjust[p, y]; // adjust survival down when off of topt
          }
          
          if (T_dep_mortality == 0) {
            surv[p, a, y] = exp(-(f_a_y[a, y] + m));
          }
        }
      }
    }
    
    //// create movement matrices ////
    
    if (T_dep_movement == 1) {
      diff_m = adj_m * d; // set up matrix as adjacency * dispersal
      
      diff_m = add_diag(diff_m, -1 * colSums(diff_m)); // rescale the columns so they will sum to 1 at the end
      
      // set up array of matrices containing 
      for (y in 1 : ny_train) {
        for (i in 1 : np) {
          for (j in 1 : np) {
            // in R this is just outer(np, np, "-") 
            
            if (exp_yn == 1) {
              T_adjust_m[y, i, j] = exp(beta_t
              * (log(T_adjust[i, y])
              - log(T_adjust[j, y])));
            } else {
              T_adjust_m[y, i, j] = fmin(500,
              exp(beta_t
              * (T_adjust[i, y]
              - T_adjust[j, y])));
            }
          }
        }
        tax_m[y] = adj_m .* T_adjust_m[y];
        
        tax_m[y] = add_diag(tax_m[y], -1 * colSums(tax_m[y])); // fill in the diagonal with within-patch "taxis" so everything sums to 1 
        
        mov_inst_m[y] = diff_m + tax_m[y]; // movement as a sum of diffusion and taxis (can cancel each other out)
        
        mov_m[y] = matrix_exp(mov_inst_m[y]); // matrix exponentiate, although see https://discourse.mc-stan.org/t/matrix-exponential-function/9595
        
        if ((sum(colSums(mov_m[y])) / np - 1) > .001) {
          print("Something has gone very wrong, movement matrix columns do not sum to 1");
          print(colSums(mov_m[y]));
          print("width is", width);
          print("Topt is", Topt);
          print(diagonal(mov_inst_m[y]));
        }
      }
    } // close T_dep_movement if 
    
    //// fill in year 1 of n_at_age_hat, initialized with mean_recruits ////
    for (p in 1 : np) {
      if (use_init_n_at_age == 0) {
        for (a in 1 : n_ages) {
          if (a == 1) {
            if (T_dep_recruitment == 1 && spawner_recruit_relationship == 0) {
              n_at_age_hat[1, p, a] = init_dep[p] * r0 * capacity[p] * T_adjust[p, 1]
              * exp(sigma_r * raw[1]
              - square(sigma_r) / 2); // initialize age 0 with mean recruitment in every patch
            }
            if (T_dep_recruitment == 0 && spawner_recruit_relationship == 0) {
              n_at_age_hat[1, p, a] = init_dep[p] * r0 * capacity[p]
              * exp(sigma_r * raw[1]
              - square(sigma_r) / 2); // initialize age 0 with mean recruitment in every patch
            }
            if (T_dep_recruitment == 0 && spawner_recruit_relationship == 1) {
              n_at_age_hat[1, p, a] = init_dep[p] * r0 * capacity[p]
              * exp(sigma_r * raw[1]
              - square(sigma_r) / 2); // scale it down a bit -- historical fishing was still occurring
            }
            if (T_dep_recruitment == 1 && spawner_recruit_relationship == 1) {
              n_at_age_hat[1, p, a] = init_dep[p] * r0 * capacity[p]
              * exp(sigma_r * raw[1]
              - square(sigma_r) / 2)
              * T_adjust[p, 1] //* beta_rec
              ;
            }
          } // close age==1 case
          else {
            n_at_age_hat[1, p, a] = n_at_age_hat[1, p, a - 1]
            * surv[p, a - 1, 1]; // initialize population with mean recruitment propogated through age classes with mortality
          }
        } // close ages
      } else {
        n_at_age_hat[1, 1 : np, 1 : n_ages] = init_n_at_age;
      }
      ssb[p, 1] = sum(to_vector(n_at_age_hat[1, p, 1 : n_ages])
      .* maturity_at_age .* wt_at_age);
    } // close patches
    
    //// run population model  ////
    
    for (y in 2 : ny_train) {
      if (y == 2) {
        rec_dev[y - 1] = sigma_r * raw[y]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
      } // close y==2 case  
      else {
        rec_dev[y - 1] = alpha * rec_dev[y - 2]
        + sqrt(1 - square(alpha)) * sigma_r * raw[y];
      } // close ifelse
      
      // describe population dynamics
      for (p in 1 : np) {
        // density-independent, temperature-dependent recruitment of age 1
        
        if (T_dep_recruitment == 1 && spawner_recruit_relationship == 0) {
          n_at_age_hat[y, p, 1] = r0 * capacity[p]
          * exp(rec_dev[y - 1] - square(sigma_r) / 2)
          * T_adjust[p, y - 1] //* beta_rec
          ;
        }
        if (T_dep_recruitment == 0 && spawner_recruit_relationship == 0) {
          n_at_age_hat[y, p, 1] = r0 * capacity[p]
          * exp(rec_dev[y - 1] - square(sigma_r) / 2);
        }
        
        if (T_dep_recruitment == 0 && spawner_recruit_relationship == 1) {
          n_at_age_hat[y, p, 1] = ((0.8 * r0 * h * sum(ssb[1:np, y - 1]))
          / (0.2 * ssb0 * (1 - h)
          + sum(ssb[1:np, y - 1]) * (h - 0.2))) * capacity[p];
          
          n_at_age_hat[y, p, 1] = n_at_age_hat[y, p, 1]
          * exp(rec_dev[y - 1] - square(sigma_r) / 2);
        }
        if (T_dep_recruitment == 1 && spawner_recruit_relationship == 1) {
          n_at_age_hat[y, p, 1] = ((0.8 * r0 * h * sum(ssb[1:np, y - 1]))
          / (0.2 * ssb0 * (1 - h)
          + sum(ssb[1:np, y - 1]) * (h - 0.2)))
          * T_adjust[p, y - 1] * capacity[p];
          
          n_at_age_hat[y, p, 1] = n_at_age_hat[y, p, 1]
          * exp(rec_dev[y - 1] - square(sigma_r) / 2);
        }
        // 
        // why estimate raw and sigma_r? we want to estimate process error
        // if we got rid of sigma_r, we would be saying that raw could be anything
        // that means raw could pick any value, and it would pick deviates to perfectly match abundance index
        // sigma_r scales the amount of process error that we say is possible
        // letting it go to infinity means the model will always fit the data perfectly 
        // raw is the realized process error 
        // exp(rec_dev[y-1] - square(sigma_r,2)/2) = random variable with mean 0 and SD sigma_r
        // allows us to have a different recruitment deviation every year even though sigma_r is constant 
      } // close recruitment by patch
      // pop dy for reproductive adults
      
      if (T_dep_movement == 0) {
        // mortality and dispersal are happening simultaneously here, between generations
        // because neither is patch-specific I don't think the order matters
        for (p in 1 : np) {
          for (a in 2 : n_ages) {
            // edge cases -- edges are reflecting
            if (p == 1) {
              n_at_age_hat[y, p, a] = n_at_age_hat[y - 1, p, a - 1]
              * surv[p, a - 1, y - 1]
              * (1 - d_at_age[a - 1])
              + n_at_age_hat[y - 1, p + 1, a - 1]
              * surv[p + 1, a - 1, y - 1]
              * d_at_age[a - 1];
            } // close patch 1 case 
            else if (p == np) {
              n_at_age_hat[y, p, a] = n_at_age_hat[y - 1, p, a - 1]
              * surv[p, a - 1, y - 1]
              * (1 - d_at_age[a - 1])
              + n_at_age_hat[y - 1, p - 1, a - 1]
              * surv[p - 1, a - 1, y - 1]
              * d_at_age[a - 1];
            } // close highest patch
            else {
              n_at_age_hat[y, p, a] = n_at_age_hat[y - 1, p, a - 1]
              * surv[p, a - 1, y - 1]
              * (1 - 2 * d_at_age[a - 1])
              + n_at_age_hat[y - 1, p - 1, a - 1]
              * surv[p - 1, a - 1, y - 1]
              * d_at_age[a - 1]
              + n_at_age_hat[y - 1, p + 1, a - 1]
              * surv[p + 1, a - 1, y - 1]
              * d_at_age[a - 1];
            } // close if/else for all other patches
          } // close ages
        } // close patches 
      } // close T-dep movement if 
      
      // this code block calculates adult population size based on survival and directional movement 
      if (T_dep_movement == 1) {
        for (p in 1 : np) {
          n_at_age_hat[y, p, 2 : n_ages] = n_at_age_hat[y - 1, p, 1 : (
            n_ages - 1)]
            .* to_row_vector(surv[p, 1 : (
              n_ages - 1), 
              y - 1]);
        }
        
        // assumption is that only mature fish move
        for (a in age_at_maturity : n_ages) {
          // some acrobatics required here, because Stan won't do matrix multiplication with an array of reals like n_at_age_hat
          // instead we do the matrix multiplication with a placeholder vector and then populate n_at_age_hat 
          
          // fill in placeholder vector with reproductive ages across patches, and do mortality   
          for (p in 1 : np) {
            v_in[p] = n_at_age_hat[y, p, a];
          }
          v_out = mov_m[y] * v_in; // redistribute each age among patches according to the movement matrix 
          
          // fill in n_at_age_hat
          for (p in 1 : np) {
            n_at_age_hat[y, p, a] = v_out[p];
          }
        } // close ages
      } // close T-dep movement if 
      
      for (p in 1 : np) {
        ssb[p, y] = sum(to_vector(n_at_age_hat[y, p, 1 : n_ages])
        .* maturity_at_age .* wt_at_age);
      }
    } // close year 2+ loop
    
    return n_at_age_hat;
  } // close popdy function
}
// close functions block

data {
  // survey data 
  
  
  int n_ages; // number of ages
  
  int np; // number of patches
  
  vector[np] patches; // easier than creating it in Stan
  
  int ny_train; // years for training
  
  int ny_proj; // number of years to forecast 
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  matrix[n_ages, n_lbins] l_at_a_key;
  
  vector[n_ages] wt_at_age;
  
  array[np, ny_train] real dens; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
  vector[np] area; // mean area swept per patch 
  
  // environmental data 
  
  matrix[np, ny_train] sbt; // temperature data for training
  
  matrix[np, ny_proj] sbt_proj;
  
  // fish data
  
  real m; // total mortality 
  
  real k;
  
  real loo;
  
  real t0;
  
  // real cv;
  
  matrix[n_ages, ny_train] f;
  
  matrix[n_ages, ny_proj + 1] f_proj;
  
  real length_50_sel_guess;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  int age_at_maturity;
  
  int<lower=0, upper=1> do_dirichlet;
  
  int<lower=0, upper=1> T_dep_recruitment;
  
  int<lower=0, upper=1> T_dep_mortality;
  
  int<lower=0, upper=1> T_dep_movement;
  
  int<lower=0, upper=1> eval_l_comps;
  
  int<lower=0, upper=1> spawner_recruit_relationship;
  
  int<lower=0, upper=1> run_forecast;
  
  int<lower=0, upper=1> exp_yn;
  
  int<lower=0, upper=1> process_error_toggle;
  
  int<lower=0, upper=1> use_poisson_link;
  
  int number_quantiles;
  
  array[number_quantiles] real quantiles_calc;
  
  array[np, n_lbins, ny_train] int n_at_length; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
  real h; // steepness
  
  real dcap; // cap on diffusion parameter
  
  // pass in all prior values for distributions 
  // all named as follows: pr (for prior) + parameter name + parameter of probability distribution as required by Stan (e.g., alpha and beta for the beta distribution, mu and sigma for the normal distribution)
  real pr_init_dep_alpha; 
  real pr_init_dep_beta; 
  real pr_beta_obs_mu;
  real pr_beta_obs_sigma;
  real pr_beta_obs_int_mu;
  real pr_beta_obs_int_sigma; 
  real pr_raw_mu;
  real pr_raw_sigma;
  real pr_sigma_r_raw_mu;
  real pr_sigma_r_raw_sigma;
  real pr_sigma_obs_mu; 
  real pr_sigma_obs_sigma; // formerly sigma_obs_cv
  real pr_log_r0_mu; 
  real pr_log_r0_sigma; 
  real pr_Topt_mu; 
  real pr_Topt_sigma; 
  real pr_width_mu; 
  real pr_width_sigma; 
  real pr_beta_t_mu; 
  real pr_beta_t_sigma; 
  real pr_beta_rec_mu; 
  real pr_beta_rec_sigma; 
  real pr_alpha_alpha; 
  real pr_alpha_beta;
  real pr_d_mu;
  real pr_d_sigma;
  real pr_p_length_50_sel_sigma; // note that there is no mu prior passed in because it's calculated in the code from length_50_sel_guess and loo 
  real pr_sel_delta_mu;
  real pr_sel_delta_sigma; 
  real pr_theta_d_mu;
  real pr_theta_d_sigma; 
  int<lower = 0, upper = 1> patch_r0s;
  
}
transformed data {
  vector[n_ages] maturity_at_age; // vector of probabilities of being mature at each age, currently binary (0/1) and taken as known
  
  matrix[np, np] adj_m; // adjacency matrix for patches 
  
  matrix[np, np] outer; // outer difference matrix 
  
  vector[np] parea = area / sum(area);
  
  //int exp_yn; 
  
  // in R this is just outer(np, np, "-") 
  for (i in 1 : np) {
    for (j in 1 : np) {
      outer[i, j] = i - j; // should just create a matrix that counts from 1 to np along rows and columns
      if (abs(outer[i, j]) == 1) {
        adj_m[i, j] = 1; // fill with 1s for adjacent patches and 0s otherwise 
      } else {
        adj_m[i, j] = 0;
      }
    }
  }
  
  // print("the outer matrix is ",outer); 
  // print("the adjacency matrix is ",adj_m); 
  
  for (a in 1 : n_ages) {
    if (a < age_at_maturity) {
      maturity_at_age[a] = 0;
    } else {
      maturity_at_age[a] = 1;
    }
  }
}
parameters {
  // real<lower = 1e-6> sigma_total;
  
  real<lower=1e-6, upper = 200> sigma_obs; // upper as a cap for edge cases where sigma_obs goes to inf
  
  real<lower=1e-6> sigma_r_raw;
  
  real<lower=0.5> width; // sensitivity to temperature variation
  
  real Topt; //  temp at which recruitment is maximized
  
  real<lower=0, upper=0.99> alpha; // autocorrelation term
  
  vector[ny_train] raw; // array of raw recruitment deviates, changed to one value per year
  
  real<upper=0.8> p_length_50_sel; // length at 50% selectivity
  
  real<lower=1e-6> sel_delta;
  
  real<lower=0> beta_obs; // controls how fast detection goes up with abundance
  
  // real<lower=0, upper=0.333> d; // dispersal fraction (0.333 = perfect admixture)
  
  real<lower=0, upper=dcap> d; // increasing bounds on this for the temperature-dependent movement model
  
  vector<lower=1e-3, upper=1>[np] init_dep;
  
  real<lower=0, upper=2> theta_d;
  
  real beta_t; // responsiveness of movement to temperature
  
//  real beta_rec; // responsivenses of mean recruits to temperature
  
  real beta_obs_int; // intercept of detection probability
  
  real log_r0;
  
  simplex[np] habitat;
}
transformed parameters {
  real length_50_sel;
  
  matrix[np, ny_train] theta; // Bernoulli probability of encounter  
  
  array[ny_train] matrix[np, n_lbins] n_at_length_hat; // array number of years containing matrices with numbers at patch, length bin, and year 
  
  array[np, ny_train] real density_hat; // for tracking sum density 
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  real ssb0;
  
  vector[n_ages] unfished;
  
  real r0;
  
  vector[ny_train] centroid;
  
  matrix[number_quantiles, ny_train] range_quantiles;
  
  // real sigma_r = sigma_total / 2;
  
  // real sigma_obs = sigma_total / 2;
  
  real sigma_r;
  
  vector[n_ages] d_at_age; // storage for diffusion at age
  
  vector[np] capacity;
  
  array[ny_train] matrix[np, n_ages] n_at_age_hat;
  
  matrix[np, n_ages] init_n_at_age;
  
  init_n_at_age = rep_matrix(0, np, n_ages);
  
  d_at_age = rep_vector(0, n_ages);
  
  unfished = rep_vector(0, n_ages);
  
  if (patch_r0s == 1){
      capacity = habitat;

  } else {
    capacity = rep_vector(1.0 / np,np);
  }
  
  for (a in 1 : n_ages) {
    if (a >= age_at_maturity) {
      d_at_age[a] = d;
    }
  }
  

  r0 = exp(log_r0);
  
  sigma_r = sigma_r_raw * process_error_toggle;
  
  ssb0 = -999;
  
  unfished[1] = r0;
  if (spawner_recruit_relationship == 1) {
    for (a in 2 : n_ages) {
      unfished[a] = unfished[a - 1] * exp(-m);
    }
    
    ssb0 = sum(unfished .* maturity_at_age .* wt_at_age);
  }
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
    // length_50_sel = loo * .01; // Dan made a note to change this sometime

  
  selectivity_at_bin = 1.0
  ./ (1
  + exp(-log(19)
  * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  n_at_age_hat = simulate_population(np, ny_train, n_ages, n_lbins,
  age_at_maturity, sbt, Topt, width,
  exp_yn, T_dep_mortality, T_dep_movement,
  f, m, d, beta_t, T_dep_recruitment,
  spawner_recruit_relationship, init_dep,
  sigma_r, raw,
  r0, maturity_at_age, wt_at_age, alpha,
  h, ssb0, d_at_age, l_at_a_key,
  selectivity_at_bin, beta_obs_int,
  beta_obs, number_quantiles,
  init_n_at_age, 0,capacity);
  
  for (y in 1 : ny_train) {
    for (p in 1 : np) {
      n_at_length_hat[y, p, 1 : n_lbins] = ((l_at_a_key'
      * to_vector(n_at_age_hat[y, p, 1 : n_ages]))
      .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      density_hat[p, y] = sum(to_vector(n_at_length_hat[y, p, 1 : n_lbins]));
      if (is_nan(density_hat[p, y])){
        density_hat[p, y] = 0;
      }
      
      if (use_poisson_link == 1){
        
        // theta[p, y] = 1 - exp(-area[p] * density_hat[p, y]);
        
        theta[p, y] = 1 - exp(-density_hat[p, y]); // not multiplied by area since implicit in population size

        
      } else {
        
        // if not using Poisson link, calculate encounter probability with lognormal dsitribution 
        
        theta[p, y] = 1
        / (1
        + exp(-(beta_obs_int
        + beta_obs * log(density_hat[p, y]))));
        
      }
    } // close patches
    
    for (q in 1 : number_quantiles) {
      // calculate every range quantile q for every year y
      range_quantiles[q, y] = calculate_range_quantile(np, patches,
      density_hat[ : , y],
      quantiles_calc[q]);
    }
    
    centroid[y] = sum(to_vector(density_hat[ : , y]) .* patches)
    / sum(to_vector(density_hat[ : , y])); // calculate center of gravity
  }
}
// close transformed parameters block

model {
  real n;
  
  vector[n_lbins] prob_hat;
  
  vector[n_lbins] prob;
  
  real dml_tmp;
  
  real test;
  
  habitat ~ dirichlet(parea);
  
  init_dep ~ beta(pr_init_dep_alpha, pr_init_dep_beta);
  
  beta_obs ~ normal(pr_beta_obs_mu, pr_beta_obs_sigma);
  
  beta_obs_int ~ normal(pr_beta_obs_int_mu, pr_beta_obs_int_sigma);
  
  raw ~ normal(pr_raw_mu, pr_raw_sigma);
  
  sigma_r_raw ~ normal(pr_sigma_r_raw_mu, pr_sigma_r_raw_sigma);
  
  sigma_obs ~ normal(pr_sigma_obs_mu, pr_sigma_obs_sigma);
  
  log_r0 ~ normal(pr_log_r0_mu, pr_log_r0_sigma);
  
  Topt ~ normal(pr_Topt_mu, pr_Topt_sigma);
  
  width ~ normal(pr_width_mu, pr_width_sigma);
  
  beta_t ~ normal(pr_beta_t_mu, pr_beta_t_sigma);
  
 // beta_rec ~ normal(pr_beta_rec_mu, pr_beta_rec_sigma);
  
  alpha ~ beta(pr_alpha_alpha, pr_alpha_beta); 
  
  d ~ normal(pr_d_mu, pr_d_sigma); 
  
  p_length_50_sel ~ normal(length_50_sel_guess / loo, pr_p_length_50_sel_sigma);
  
  sel_delta ~ normal(pr_sel_delta_mu, pr_sel_delta_sigma);
  
  theta_d ~ normal(pr_theta_d_mu, pr_theta_d_sigma);
  
  for (y in 2 : ny_train) {
    for (p in 1 : np) {
      if (dens[p, y] > 0) {
        if (eval_l_comps == 1) {
          if (sum(n_at_length[p, 1 : n_lbins, y]) > 0) {
            if (do_dirichlet == 1) {
              prob_hat = to_vector(n_at_length_hat[y, p, 1 : n_lbins])
              / sum(to_vector(n_at_length_hat[y, p, 1 : n_lbins]));
              
              prob = to_vector(n_at_length[p, 1 : n_lbins, y])
              / sum(to_vector(n_at_length[p, 1 : n_lbins, y]));
              
              n = sum(n_at_length[p, 1 : n_lbins, y]);
              
              dml_tmp = lgamma(n + 1) - sum(lgamma(n * prob + 1))
              + lgamma(theta_d * n) - lgamma(n + theta_d * n)
              + sum(lgamma(n * prob + theta_d * n * prob_hat)
              - lgamma(theta_d * n * prob_hat)); // see https://github.com/merrillrudd/LIME/blob/9dcfc7f7d5f56f280767c6900972de94dd1fea3b/src/LIME.cpp#L559 for log transformation of dirichlet-multinomial in Thorston et al. 2017
              
              target += dml_tmp;
            } else {
              n_at_length[p, 1 : n_lbins, y] ~ multinomial(to_vector(
                n_at_length_hat[y, p, 1 : n_lbins])
                / sum(to_vector(
                  n_at_length_hat[y, p, 1 : n_lbins])));
            } // close dirichlet statement
          } // close if any positive length comps
        } // close eval_length_comps
        
        if (use_poisson_link == 1){
          
          // if (theta[p,y] > 0 && density_hat[p,y] > 1e-6){
          
          log(dens[p, y]) ~ normal(log(density_hat[p, y] / (theta[p,y] + 1e-6)) - square(sigma_obs)/2, sigma_obs);
   
          // }
        } else {
          
          // print("d is ",density_hat[p, y]);
          // 
          // print("theta is ",theta[p, y]);
          // 
          // print("habitat is", capacity[p]);

            log(dens[p, y]) ~ normal(log((density_hat[p, y] + 1e-6) / (theta[p,y] + 1e-6)) - square(sigma_obs)/2, sigma_obs);
            
        }
        
        1 ~ bernoulli(theta[p, y]);
      } else {
        // only evaluate density if there are length comps to evaluate
        
        0 ~ bernoulli(theta[p, y]);
      } // close else 
    } // close patch loop
  } // close year loop
}
generated quantities {
  matrix[np, ny_train] dens_pp;
  array[ny_proj + 1] matrix[np, n_lbins] n_at_length_obs_proj;
  array[ny_proj + 1] matrix[np, n_ages] n_at_age_proj; 
  array[ny_proj + 1] matrix[np, n_lbins] n_at_length_proj;
  array[np, ny_proj + 1] real density_obs_proj;
  array[np, ny_proj + 1] real density_proj;
  vector[ny_proj + 1] total_density_proj; 
  vector[ny_proj + 1] raw_proj;
  // array[ny_proj] matrix[np, n_lbins] proj_n_at_length_hat;
  vector[ny_proj + 1] centroid_proj;
  matrix[number_quantiles, ny_proj + 1] range_quantiles_proj;
  matrix[np, ny_proj + 1] theta_proj; // Bernoulli probability of encounter  
  
  //// generate posterior predictive distributions for training data ////
  
  for (y in 1 : ny_train) {
    for (p in 1 : np) {
      // ignoring error around length sampling for now
      // add in poisson toggle here
      
      if (use_poisson_link == 1){
        
        if (theta[p,y] > 0){
          
          dens_pp[p, y] = bernoulli_rng(theta[p, y])
          * exp(normal_rng(log(density_hat[p, y]) - square(sigma_obs)/2,
          sigma_obs));
        } else {
          dens_pp[p, y] = 0;
        }
      } else {
        
        dens_pp[p, y] = bernoulli_rng(theta[p, y])
        * exp(normal_rng(log(density_hat[p, y]) - square(sigma_obs) / 2,
        sigma_obs));
      }
      
    }
  }
  
  if (run_forecast == 1) {
    for (y in 1 : (ny_proj + 1)) {
      raw_proj[y] = normal_rng(0, 1); // draw a raw value 
    }
    
    n_at_age_proj = simulate_population(np, ny_proj + 1, n_ages, n_lbins,
    age_at_maturity, sbt_proj, Topt,
    width, exp_yn, T_dep_mortality,
    T_dep_movement, f_proj, m, d, beta_t,
    T_dep_recruitment,
    spawner_recruit_relationship,
    init_dep,// beta_rec,
    sigma_r, raw_proj, r0,
    maturity_at_age, wt_at_age, alpha, h,
    ssb0, d_at_age, l_at_a_key,
    selectivity_at_bin, beta_obs_int,
    beta_obs, number_quantiles,
    n_at_age_hat[ny_train,  : ,  : ], 1,capacity);
    
    for (y in 1 : (ny_proj + 1)) {
      
      for (p in 1 : np) {
        n_at_length_proj[y, p, 1 : n_lbins] = ((l_at_a_key'
        * to_vector(n_at_age_proj[y, p, 1 : n_ages]))
        .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
        // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
        
        // NOTE, ignoring length sampling process at this point. In theory, need multinomial_rng or a custom dirichlet-multinomial rng here to generate simulated length comps
        
        n_at_length_obs_proj[y, p, 1 : n_lbins] = n_at_length_proj[y, p, 1 : n_lbins]; // this is where length sampling process would go
        
        density_proj[p, y] = sum(to_vector(n_at_length_proj[y, p, 1 : n_lbins])); // true population
        
        
      if (use_poisson_link == 1){
        
        // theta_proj[p, y] = 1 - exp(-area[p] * density_proj[p, y]);
        
        theta_proj[p, y] = 1 - exp(-density_proj[p, y]);

      } else {
        theta_proj[p, y] = 1
        / (1
        + exp(-(beta_obs_int
        + beta_obs
        * log(density_proj[p, y]))));
      }
        
        density_obs_proj[p, y] = bernoulli_rng(theta_proj[p, y])
        * exp(normal_rng(log((density_proj[p, y] / theta_proj[p, y])) - square(sigma_obs) / 2,
        sigma_obs));
        // the observed densitieis as opposed to the true densities
      } // close patches 
      
      for (q in 1 : number_quantiles) {
        // calculate every range quantile q for every year y
        range_quantiles_proj[q, y] = calculate_range_quantile(np, patches,
        density_obs_proj[ : , y],
        quantiles_calc[q]);
      }
      
      centroid_proj[y] = sum(to_vector(density_obs_proj[ : , y]) .* patches)
      / sum(to_vector(density_obs_proj[ : , y])); // calculate center of gravity
      
      total_density_proj[y] = sum(density_obs_proj[ : , y]); // get total abundance for summary statistics 
      
    } // close projection years loop 
  } // close run_forecast
}
// close generated quantities run_forecast

