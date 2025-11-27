# Robust Site Selection: single-file runnable script
# Save this as e.g. run_dispersion.R and source() or paste into R console.

# ---- 1) Write C++ source and compile with Rcpp ----
cpp_code <- '
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
using namespace Rcpp;

// ---------- Helper non-exported functions ----------
static double calc_objective_sum_c(const NumericMatrix &dist, const IntegerVector &selected) {
  int n = selected.size();
  double total = 0.0;
  for(int i = 0; i < n - 1; ++i) {
    for(int j = i + 1; j < n; ++j) {
      int ii = selected[i] - 1;
      int jj = selected[j] - 1;
      total += dist(ii, jj);
    }
  }
  return total;
}

static double calc_objective_maxmin_c(const NumericMatrix &dist, const IntegerVector &selected) {
  int n = selected.size();
  double min_dist = R_PosInf;
  for(int i = 0; i < n - 1; ++i) {
    for(int j = i + 1; j < n; ++j) {
      int ii = selected[i] - 1;
      int jj = selected[j] - 1;
      double d = dist(ii, jj);
      if(d < min_dist) min_dist = d;
    }
  }
  return min_dist;
}

static double calc_delta_sum_c(const NumericMatrix &dist, const IntegerVector &selected, int idx_sel, int i_new) {
  int i_old = selected[idx_sel];
  int n_sel = selected.size();
  double delta = 0.0;

  int old_0idx = i_old - 1;
  int new_0idx = i_new - 1;

  for(int k = 0; k < n_sel; ++k) {
    if(k == idx_sel) continue;
    int other_0idx = selected[k] - 1;
    delta -= dist(old_0idx, other_0idx);
    delta += dist(new_0idx, other_0idx);
  }

  return delta;
}

// ---------- Exported objective wrappers ----------
 // [[Rcpp::export]]
double calc_objective_sum(NumericMatrix dist, IntegerVector selected) {
  return calc_objective_sum_c(dist, selected);
}

 // [[Rcpp::export]]
double calc_objective_maxmin(NumericMatrix dist, IntegerVector selected) {
  return calc_objective_maxmin_c(dist, selected);
}

// ---------- Exported local search swap ----------
 // [[Rcpp::export]]
List local_search_swap(NumericMatrix dist, IntegerVector selected, IntegerVector candidates, std::string objective, int max_iter) {
  int n_sel = selected.size();
  int n_cand = candidates.size();
  double current_obj = -R_PosInf;

  if(objective == "sum") {
    current_obj = calc_objective_sum_c(dist, selected);
  } else {
    current_obj = calc_objective_maxmin_c(dist, selected);
  }

  bool improved = true;
  int iter = 0;

  while(improved && iter < max_iter) {
    improved = false;

    for(int i = 0; i < n_sel; ++i) {
      for(int j = 0; j < n_cand; ++j) {
        double new_obj;

        if(objective == "sum") {
          double delta = calc_delta_sum_c(dist, selected, i, candidates[j]);
          new_obj = current_obj + delta;
        } else {
          IntegerVector temp_selected = clone(selected);
          temp_selected[i] = candidates[j];
          new_obj = calc_objective_maxmin_c(dist, temp_selected);
        }

        if(new_obj > current_obj) {
          selected[i] = candidates[j];
          current_obj = new_obj;
          improved = true;
          break; // restart scanning selected positions
        }
      }
      if(improved) break;
    }

    ++iter;
  }

  return List::create(
    Named("selected") = selected,
    Named("objective") = current_obj,
    Named("iterations") = iter
  );
}
'

# write to temp file and compile
tmp_cpp <- tempfile(fileext = ".cpp")
writeLines(cpp_code, tmp_cpp)
message("Compiling C++ (Rcpp)...")
if(!requireNamespace("Rcpp", quietly = TRUE)) stop("Please install package 'Rcpp' first.")
Rcpp::sourceCpp(tmp_cpp)
message("C++ compiled and functions exported: calc_objective_sum, calc_objective_maxmin, local_search_swap")

# ---- 2) R helper functions 

coerce_to_3d_array <- function(distances, n_sites) {
  if (is.null(distances)) stop("'distances' is NULL")
  # vector
  if (is.atomic(distances) && is.vector(distances) && !is.list(distances) && !is.matrix(distances) && !is.array(distances)) {
    if (length(distances) == n_sites * n_sites) {
      mat <- matrix(as.numeric(distances), nrow = n_sites, ncol = n_sites)
      return(array(mat, dim = c(n_sites, n_sites, 1)))
    } else {
      stop(sprintf("distance vector length %d doesn't match %d sites (expected %d)", length(distances), n_sites, n_sites*n_sites))
    }
  }
  # matrix
  if (is.matrix(distances)) {
    if (nrow(distances) != n_sites || ncol(distances) != n_sites) stop("distance matrix dims != n_sites")
    return(array(as.numeric(distances), dim = c(n_sites, n_sites, 1)))
  }
  # array
  if (is.array(distances)) {
    dims <- dim(distances)
    if (length(dims) == 2) {
      if (dims[1] != n_sites || dims[2] != n_sites) stop("2D array dims mismatch")
      return(array(as.numeric(distances), dim = c(n_sites, n_sites, 1)))
    } else if (length(dims) == 3) {
      if (dims[1] != n_sites || dims[2] != n_sites) stop("3D array dims mismatch")
      return(distances)
    } else {
      stop("distances array must be 2D or 3D")
    }
  }
  stop("Unsupported 'distances' type")
}

great_circle_distance <- function(lat1, lon1, lat2, lon2) {
  # Haversine in km
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad; lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad; lon2 <- lon2 * to_rad
  dlat <- lat2 - lat1; dlon <- lon2 - lon1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * asin(pmin(1, sqrt(a)))
  r <- 6371 # radius of earth. 
  return(r * c)
}

update_distances_jitter <- function(distances, sites_df, uncertain_idx) {
  n <- nrow(sites_df)
  dist_arr <- coerce_to_3d_array(distances, n)
  K <- dim(dist_arr)[3]
  dist_boot <- array(dist_arr, dim = c(n, n, K))
  coords_boot <- data.frame(
    lat = sites_df$lat, 
    lon = sites_df$lon
  )

  if (length(uncertain_idx) > 0) {
    for (i in uncertain_idx) {
      jitter_sd <- sites_df$coord_uncertainty[i]
      if (is.na(jitter_sd) || jitter_sd <= 0) next
      coords_boot$lat[i] <- coords_boot$lat[i] + rnorm(1, 0, jitter_sd)
      coords_boot$lon[i] <- coords_boot$lon[i] + rnorm(1, 0, jitter_sd)
    }

    # recompute primary distance matrix
    for (i in uncertain_idx) {
      for (j in seq_len(n)) {
        d <- great_circle_distance(
          coords_boot$lat[i], coords_boot$lon[i],
          coords_boot$lat[j], coords_boot$lon[j]
        )
        dist_boot[i, j, 1] <- d
        dist_boot[j, i, 1] <- d
      }
    }
    # perturb other matrices proportionally (if present)
    if (K >= 2) {
      for (k in 2:K) {
        for (i in uncertain_idx) {
          jitter_factor <- rnorm(1, 1, 0.05)
          dist_boot[i, , k] <- dist_boot[i, , k] * jitter_factor
          dist_boot[, i, k] <- dist_boot[, i, k] * jitter_factor
        }
      }
    }
  }
  dist_boot
}

# Greedy initialization
greedy_initialize <- function(dist_matrix, n_sites, seeds = integer(0), objective = "sum") {
  n_total <- nrow(dist_matrix)
  selected <- as.integer(seeds)
  if (length(selected) > n_sites) selected <- selected[1:n_sites]
  candidates <- setdiff(seq_len(n_total), selected)

  while (length(selected) < n_sites && length(candidates) > 0) {
    best_site <- NULL
    best_obj <- -Inf
    for (site in candidates) {
      test_selected <- c(selected, site)
      if (objective == "sum") {
        obj <- calc_objective_sum(dist_matrix, as.integer(test_selected))
      } else {
        obj <- calc_objective_maxmin(dist_matrix, as.integer(test_selected))
      }
      if (obj > best_obj) { best_obj <- obj; best_site <- site }
    }
    if (is.null(best_site)) break
    selected <- c(selected, best_site)
    candidates <- setdiff(candidates, best_site)
  }
  as.integer(selected)
}

# Optimiation function
maximize_dispersion_robust <- function(
  input_data,
  n_sites = 5,
  weight_1 = 1.0,
  weight_2 = 0.0,
  n_bootstrap = 99,
  dropout_prob = 0.15,
  objective = c("sum", "maxmin"),
  n_local_search_iter = 100,
  n_restarts = 3,
  track_top_n = 10,
  seed = NULL,
  verbose = TRUE
) {

  objective <- match.arg(objective)
  if(!is.null(seed)) set.seed(seed)

  cooccur <- matrix(0, nrow(sites_df), nrow(sites_df))
  distances <- input_data$distances
  sites_df <- input_data$sites
  n_total <- nrow(sites_df)
  if(n_total <= 0) stop("No sites provided")

  # coerce to 3D
  distances <- coerce_to_3d_array(distances, n_total)

  # combine first two matrices by weights (safe)
  dist_combined <- distances[,,1] * weight_1
  if (dim(distances)[3] >= 2 && weight_2 > 0) {
    dist_combined <- dist_combined + distances[,,2] * weight_2
  }

  # check seeds and uncertain
  if(!"required" %in% names(sites_df)) sites_df$required <- FALSE
  sites_df$required <- as.logical(sites_df$required)
  seeds <- which(sites_df$required)
  uncertain_idx <- which(!is.na(sites_df$coord_uncertainty) & sites_df$coord_uncertainty > 0)

  if(verbose) cat(sprintf("Sites: %d | Seeds: %d | Uncertain: %d\n", n_total, length(seeds), length(uncertain_idx)))

  selection_counts <- integer(n_total)
  all_solutions <- list()
  solution_counter <- 1

  for (b in seq_len(n_bootstrap)) {
    available_sites <- seq_len(n_total)

    # dropout (non-seed only)
    if (dropout_prob > 0) {
      droppable <- setdiff(available_sites, seeds)
      n_drop <- floor(length(droppable) * dropout_prob)
      if (n_drop > 0 && length(droppable) >= n_drop) {
        dropped <- sample(droppable, n_drop)
        available_sites <- setdiff(available_sites, dropped)
      }
    }

    # jitter distances if needed
    if (length(uncertain_idx) > 0) {
      dist_boot_array <- update_distances_jitter(distances, sites_df, uncertain_idx)
      dist_boot <- dist_boot_array[,,1] * weight_1
      if(dim(dist_boot_array)[3] >= 2 && weight_2 > 0) dist_boot <- dist_boot + dist_boot_array[,,2] * weight_2
    } else {
      dist_boot <- dist_combined
    }

    best_solution <- NULL
    best_objective <- -Inf

    for (restart in seq_len(n_restarts)) {
      current_solution <- greedy_initialize(dist_boot, n_sites, seeds, objective)
      # ensure selection length
      if (length(current_solution) < n_sites) {
        # try to pad with available_sites
        extra <- setdiff(available_sites, current_solution)
        if (length(extra) > 0) {
          needed <- n_sites - length(current_solution)
          current_solution <- c(current_solution, head(extra, needed))
        }
      }
      candidates <- setdiff(available_sites, current_solution)
      if (length(candidates) > 0 && length(current_solution) == n_sites) {
        res <- local_search_swap(dist_boot, as.integer(current_solution), as.integer(candidates), objective, n_local_search_iter)
        sol <- as.integer(res$selected)
        objv <- as.numeric(res$objective)
        if (objv > best_objective) {
          best_objective <- objv
          best_solution <- sol
        }
      } else {
        # compute objective for current_solution
        if (length(current_solution) > 0) {
          if (objective == "sum") {
            curv <- calc_objective_sum(dist_boot, as.integer(current_solution))
          } else {
            curv <- calc_objective_maxmin(dist_boot, as.integer(current_solution))
          }
          if (curv > best_objective) {
            best_objective <- curv
            best_solution <- as.integer(current_solution)
          }
        }
      }
    } # restarts

    if (!is.null(best_solution)) {
      idx <- best_solution
      cooccur[idx, idx] <- cooccur[idx, idx] + 1
  }
    
    if(verbose && (b %% max(1, floor(n_bootstrap/10)) == 0)) cat(sprintf("Bootstrap %d/%d\n", b, n_bootstrap))
  } # bootstrap

  # final greedy + local search on combined (non-boot) distances
  final_solution <- greedy_initialize(dist_combined, n_sites, seeds, objective)
  candidates_final <- setdiff(seq_len(n_total), final_solution)
  if (length(final_solution) == n_sites && length(candidates_final) > 0) {

    final_res <- local_search_swap(
      dist_combined,
      as.integer(final_solution),
      as.integer(candidates_final), 
      objective,
      n_local_search_iter * 2)
    
    final_solution <- as.integer(final_res$selected)
    final_objective <- as.numeric(final_res$objective)
  } else {
    if (objective == "sum")
      final_objective <- calc_objective_sum(dist_combined, as.integer(final_solution)) else 
      final_objective <- calc_objective_maxmin(dist_combined, as.integer(final_solution))
  }

  cooccurrence_strength <- rowSums(cooccur)

  stability <- data.frame(
    site_id = sites_df$site_id,
    cooccur_strength = cooccurrence_strength,
    is_seed = sites_df$required
  )

  # Highest "teamwork" sites first
  stability <- stability[order(-stability$cooccur_strength), ]

  # top solutions aggregation
  top_solutions <- data.frame(
      rank = integer(), 
      sites = I(list()), 
      objective = numeric(), 
      frequency = numeric()
    )

  if (length(all_solutions) == 0) {
    top_solutions 
  } else {
    solution_strings <- sapply(all_solutions, function(x) paste(sort(x$sites), collapse = "-"))
    solution_objectives <- sapply(all_solutions, function(x) x$objective)
    unique_solutions <- unique(solution_strings)
    for (i in seq_along(unique_solutions)) {
      idx <- which(solution_strings == unique_solutions[i])
      sites_vector <- all_solutions[[idx[1]]]$sites
      top_solutions <- rbind(
        top_solutions, 
        data.frame(
          rank = i, 
          sites = I(list(sites_vector)), 
          objective = mean(solution_objectives[idx]), 
          frequency = length(idx)/max(1,n_bootstrap)
        )
      )
    }
    top_solutions <- top_solutions[order(-top_solutions$frequency, -top_solutions$objective), ]
    if (nrow(top_solutions) > track_top_n) top_solutions <- top_solutions[1:track_top_n, ]
  }

  # return these objects back to the user. 
  list(
    selected_sites = final_solution,
    objective_value = final_objective,
    stability = stability,
    top_solutions = top_solutions,
    settings = list(
      n_sites = n_sites,
      n_bootstrap = n_bootstrap, 
      objective = objective, 
      dropout_prob = dropout_prob, 
      n_uncertain = length(uncertain_idx)
    )
  )
}



# Examples 

# function to create dummy data. 
create_example_data <- function(n_sites = 30, seed = 101) {
  set.seed(seed)
  df <- data.frame(
    site_id = seq_len(n_sites),
    lat = runif(n_sites, 25, 50),
    lon = runif(n_sites, -125, -65),
    required = FALSE,
    coord_uncertainty = 0
  )

  ## identify the center of all points.
  # the algorithm requires the center of gravity
  # otherwise it will just select a 'ring' of sites around the
  # exterior of occurrences.  
  dists2c <- great_circle_distance(
    mean(df$lat), 
    mean(df$lon), 
    df$lat, 
    df$lon
  )

  ### any of these can be the seeds. 
  df[order(dists2c)[1:3],'required'] <- TRUE

  # a few required seeds
  uncertain_sites <- sample(setdiff(seq_len(n_sites), which(df$required)), size = min(6, n_sites-3))
  df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 0.001, 0.01)
  # distance matrix

  n <- nrow(df)
  dist_mat <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- great_circle_distance(df$lat[i], df$lon[i], df$lat[j], df$lon[j])
      dist_mat[i,j] <- d; dist_mat[j,i] <- d
    }
  }
  list(distances = dist_mat, sites = df)
}

# small quick run (fast) to verify everything
test_data <- create_example_data(n_sites = 50, seed = 23)
res <- maximize_dispersion_robust(
  input_data = test_data,
  n_sites = 25,
  n_bootstrap = 99,
  dropout_prob = 0.10,
  objective = "sum",
  n_local_search_iter = 50,
  n_restarts = 2,
  seed = 42,
  verbose = TRUE
)
print(res$selected_sites)
print(res$objective_value)
print(head(res$stability, 10))

## plot results
library(ggplot2)

plottable = merge(test_data$sites, res$stability, on = 'site_id', how = 'left')
plottable = merge(
  plottable,
  data.frame(
    site_id = res$selected_sites, 
    selected = T
  ), on = 'site_id',all.x = TRUE
)

ggplot(data = plottable, 
  aes(
    x = lon, 
    y = lat, 
    shape = required, 
    size = selection_freq,
    color = selected
    )
  ) +
  geom_point()
