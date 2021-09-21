
#' IUB codes used by CutFree
#' @export
IUB_CODES <- list(
  A = c("A"),
  C = c("C"),
  G = c("G"),
  T = c("T"),
  R = c("A", "G"),
  Y = c("C", "T"),
  M = c("A", "C"),
  K = c("G", "T"),
  S = c("C", "G"),
  W = c("A", "T"),
  H = c("A", "C", "T"),
  B = c("C", "G", "T"),
  V = c("A", "C", "G"),
  D = c("A", "G", "T"),
  N = c("A", "C", "G", "T")
)

printf <- function(fmt, ...) {
  cat(sprintf(fmt, ...))
}

# Find codes that do not contain a code.
# Example: if code is A, then blocking codes are
#   C,G,T,Y,K,S,B (all codes without A)
# If names==TRUE, return the codes (character vector of bases);
# otherwise return the code names.
get_blocking_codes <- function(code, codes=IUB_CODES, names=FALSE) {
  to_block <- codes[code]
  blocking <- Filter(Negate(function(bases) any(to_block %in% bases)), codes)
  if (names) {
    return(names(blocking))
  } else {
    return(blocking)
  }
}

# Turn a string into a vector of single characters.
str_to_vector <- function(str) {
  if (length(str) == 1) {
    strsplit(str, "")[[1]]
  } else {
    str
  }
}

#' Text summary of degenerate oligos
#'
#' Create textual representations of an oligo describing the base content at
#' each position.
#'
#' \code{make_oligo_block} creates a text block showing which bases are allowed
#' by the oligo. It returns four strings denoting the presence
#' of each base.
#'
#' \code{print_oligo_block} prints the output of
#' \code{make_oligo_block} plus the oligo.
#'
#' @examples
#' make_oligo_block("ACGTRYMKSWHBVDN")
#'
#' @export
make_oligo_block <- function(oligo, codes=IUB_CODES) {
  oligo <- str_to_vector(oligo)

  bases <- c("A", "C", "G", "T")
  strs <- c(A="", C="", G="", T="")
  for (base in bases) {
    strs[base] <- paste(sapply(oligo, function(code) if (base %in% codes[[code]]) base else "-"),
                        collapse="")
  }
  return(strs)
}

#' @rdname make_oligo_block
#'
#' @examples
#' print_oligo_block("ACGTRYMKSWHBVDN")
#'
#' @export
print_oligo_block <- function(oligo, indent="") {
  cat(paste0(indent, oligo, collapse=""), "\n")
  cat(paste0(indent, make_oligo_block(oligo), collapse="\n"))
}

`%<in>%` <- function(x,y) all(x %in% y) & all(y %in% x)
which_name <- function(tf) names(tf)[which(tf)]

subcodes <- function(code, codes=IUB_CODES) {
  which_name(sapply(codes, function(x) all(x %in% codes[[code]])))
}

#' Reversing and complementing degenerate oligos
#'
#' @describeIn complement Complement degenerate oligos
#'
#' @export
complement <- function(strings, codes=IUB_CODES) {
  base_comps <- c(A="T", C="G", G="C", T="A")
  complemented <- lapply(codes, function(x) unname(base_comps[x]))
  comp_codes <- sapply(codes, function(x) "")
  for (code in names(comp_codes)) {
    comp_codes[code] <- which_name(sapply(codes, `%<in>%`, complemented[[code]]))
  }
  # comp_codes now contains the complement code for each code in codes
  # For IUB_CODES:
  #   A   C   G   T   R   Y   M   K   S   W   H   B   V   D   N
  #  "T" "G" "C" "A" "Y" "R" "K" "M" "S" "W" "D" "V" "B" "H" "N"

  chartr(paste(names(comp_codes), collapse=""),
         paste(comp_codes, collapse=""),
         strings)
}

#' @describeIn complement Reverse a degenerate oligo
#'
#' @export
reverse <- function(strings) {
  sapply(lapply(strsplit(strings, NULL), rev), paste, collapse="")
}

#' @describeIn complement Reverse complement of a degenerate oligo
#'
#' @export
reverse_complement <- function(strings, codes=IUB_CODES) reverse(complement(strings, codes=codes))

expand_asymmetric <- function(sites) {
  unique(c(sites, reverse_complement(sites)))
}

str2char <- function(x) strsplit(x, NULL)[[1]]
char2str <- function(x) paste(x, collapse="")

#' Calculate number of sequences in a degerate oligo
#'
#' @export
degeneracy <- function(sequence, codes=IUB_CODES) {
  dups <- vapply(codes, length, 0)
  sapply(sequence, function(x) prod(dups[str2char(x)]))
}

#' Design degenerate oligos free of restriction sites
#'
#' Design pools of DNA oligos that are guaranteed to not
#' contain specified restriction enzyme cut sites. All cut
#' sites are blocked while ensuring the optimal number of oligos
#' remain in the pool.
#'
#' @param len Length of the random oligo. Must be provided if
#' \code{starting_oligo} is not given.
#' @param sites Character vector of restriction sites to block.
#' @param starting_oligo Starting oligo from which sites will be
#' removed. If not given, defaults to a string of N's.
#' @param min_blocks Minimum number of blocks at each site, i.e.
#' the minimum number of changes that need to be made for a cut
#' site to appear anywhere in the oligo.
#' @param obj_weights Objective function weights for each code's
#' degeneracy.
#' @param re_randomize If \code{TRUE} (default), re-run the optimization
#' to randomize the codes while maintaining the same number of oligos.
#' @param obj_frac Fraction of the objective function that must be
#' maintained during re-randomization.
#' @param seed Seed for the random number generator.
#' @param quiet Run silently with no output from Gurobi.
#' @param maxtime Maximum time (in seconds) allowed for the solver.
#'
#' @return A list containing \code{code} -- the randomize oligo, and a
#' fields describing the MILP and the solution statistics.
#'
#' @examples
#' result <- cutfree(m=20, sites=c("GGTCTC", "GATATC"), min_blocks=2)
#' print_degenerate_oligo(result$code)
#'
#' @export
cutfree <- function(len=20, sites=c(), codes=IUB_CODES,
                    starting_oligo=strrep("N",len),
                    min_blocks=1,
                    obj_weights=log(1:4),
                    re_randomize=TRUE,
                    seed=NULL,
                    obj_frac=1.0,
                    quiet=FALSE,
                    maxtime=30) {
  m <- len
  sites <- expand_asymmetric(sites)  # include both strands if asymmetric
  ks <- sapply(sites, nchar)  # length of each site
  nB <- length(codes)  # number of variables for position in randomer
  starting_oligo <- str2char(starting_oligo)
  m <- length(starting_oligo)  # number of positions in randomer

  ncons <- m + sum(m - ks + 1)
  nvars <- nB * m
  idx <- function(b,i) nB*(i-1) + which(names(codes) == b)

  A <- matrix(0, nrow=ncons, ncol=nvars,
              dimnames=list(NULL,
                            paste0(names(codes), rep(1:m, each=nB))))
  rhs = rep(0, ncons)
  sense <- rep(">", ncons)
  sense[1:m] <- "="
  lb = numeric(nvars)
  ub = numeric(nvars)

  # open up only codes allowed by the starting oligo
  for (i in 1:m) {
    subs <- subcodes(starting_oligo[i])
    for (b in subs) {
      ub[idx(b,i)] <- 1
    }
  }

  # allow only a single code per location
  for (i in 1:m) {
    A[i,(nB*(i-1)+1):(nB*i)] <- 1
    rhs[i] <- 1
  }

  # add constraints to remove the RS
  curr_row <- m + 1
  for (s in seq_along(sites)) {
    rs <- str2char(sites[s])
    for (j in 1:(m-ks[s]+1)) {
      for (i in 1:ks[s]) {
        blocked <- get_blocking_codes(rs[i], names=TRUE)
        for (b in blocked) {
          A[curr_row,idx(b,i+j-1)] <- 1
        }
      }
      rhs[curr_row] <- min_blocks
      curr_row <- curr_row + 1
    }
  }

  # define the objective
  codelens <- obj_weights[vapply(codes, length, 0)]
  obj <- rep(codelens, m)

  model <- list(
    A = A,
    obj = obj,
    modelsense = "max",
    rhs = rhs,
    lb = lb,
    ub = ub,
    sense = sense,
    vtype = "B"
  )
  params <- list(
    OutputFlag=as.integer(!quiet),
    TimeLimit=maxtime
  )

  times <- system.time(
    orig_result <- gurobi::gurobi(model, params)
  )

  if (re_randomize && orig_result$status %in% c("OPTIMAL", "TIME_LIMIT")) {
    model2 <- add_objective_constraint(model, orig_result, frac=obj_frac)
    set.seed(seed)
    model2$obj <- runif(length(model2$obj))
    result <- gurobi::gurobi(model2, params)
  } else {
    result <- orig_result
  }

  if (result$status %in% c("OPTIMAL", "TIME_LIMIT"))
    code <- paste(rep(names(codes),m)[result$x[1:nvars] == 1], collapse="")
  else
    code <- NA

  return(list(code = code,
              model = model,
              time = times,
              sol_orig = orig_result,
              sol_final = result))
}

create_barcodes <- function(randomer, n=100, min_distance=5,
                            file=NULL, codes=IUB_CODES,
                            seed=NULL) {
  randomer <- str2char(randomer)
  m <- length(randomer)

  nvars <- 4*m
  ncons <- m + n
  nuc_names <- c("A", "C", "G", "T")
  idx <- function(b,i) 4*(i-1) + which(nuc_names == b)

  model = list(
    modelsense = "max",
    vtype = "B"
  )
  model$A <- matrix(0, nrow=ncons, ncol=nvars,
                    dimnames=list(NULL,
                                  paste0(nuc_names, rep(1:m, each=4))))
  model$rhs = rep(0, ncons)
  model$rhs[1:m] <- 1
  model$sense <- rep("=", ncons)
  model$lb = numeric(nvars)
  model$ub = numeric(nvars)

  # open up only codes allowed by the starting oligo
  for (i in 1:m) {
    allowed <- codes[randomer[i]]
    for (b in allowed) {
      model$ub[idx(b,i)] <- 1
    }
  }

  # allow only a single code per location
  for (i in 1:m) {
    model$A[i,(4*(i-1)+1):(4*i)] <- 1
  }

  barcodes <- as.character(rep(NA, n))

  params <- list(OutputFlag=0)
  set.seed(seed)
  for (i in 1:n) {
    model$obj <- runif(nvars) - 0.5
    elapsed <- system.time(
      result <- gurobi::gurobi(model, params),
      gcFirst = FALSE
    )
    if (result$status != "OPTIMAL") {
      printf("   No additional barcodes possible.\n")
      return(barcodes)
    }
    barcodes[i] <- paste(rep(nuc_names,m)[result$x == 1], collapse="")

    printf("   Found barcode %i (%s) in %f seconds.\n", i, barcodes[i], elapsed[3])

    if (i < n) {
      # add existing barcode as constraint
      model$A[m+i, ] <- result$x
      model$sense[m+i] <- "<"
      model$rhs[m+i] <- m - min_distance
    }
  }
  return(barcodes)
}

#' @export
find_match <- function(query, subject, codes=IUB_CODES) {
  query_bases <- codes[str_to_vector(query)]
  subject_bases <- codes[str_to_vector(subject)]
  nq <- length(query_bases)
  ns <- length(subject_bases)
  nm <- ns - nq + 1  # number of match positions to try
  is_match <- !logical(nm)
  for (i in 1:nm) {
    for (j in 1:nq) {
      if (!any(query_bases[[j]] %in% subject_bases[[i+j-1]])) {
        is_match[i] <- FALSE
        break
      }
    }
  }
  return(which(is_match))
}


 result <- cutfree(m=20, sites=c("GGTCTC", "GATATC"), min_blocks=2)
 #bcs <- create_barcodes("NNNNNNN", min_distance=3, n=100)

 