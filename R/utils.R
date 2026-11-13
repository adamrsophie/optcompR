# Internal helpers shared by the modeling layer.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Sum values v into n bins given by integer index i (duplicates add).
acc_sum <- function(i, v, n) {
  out <- numeric(n)
  if (length(i) == 0L) return(out)
  s <- rowsum(v, i, reorder = TRUE)
  out[sort(unique(i))] <- s[, 1L]
  out
}

# Gather triplet blocks: for each new position r, copy all triplets whose
# element index equals src[r]. Returns the positions into the triplet arrays
# and the new element index of each copied triplet.
gather_index <- function(li, len, src) {
  if (length(li) == 0L || length(src) == 0L) {
    return(list(pos = integer(), i = integer()))
  }
  if (length(li) == len && li[1L] == 1L && li[len] == len &&
      !is.unsorted(li, strictly = TRUE)) {
    # one triplet per element in order (plain variable blocks)
    return(list(pos = src, i = seq_along(src)))
  }
  counts <- tabulate(li, nbins = len)
  ord <- order(li)
  start <- cumsum(c(0L, counts))[seq_len(len)]
  cnt <- counts[src]
  total <- sum(cnt)
  if (total == 0L) return(list(pos = integer(), i = integer()))
  r <- rep.int(seq_along(src), cnt)
  off <- sequence(cnt) - 1L
  list(pos = ord[start[src[r]] + off + 1L], i = r)
}

check_scalar_string <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a single non-empty string.", what), call. = FALSE)
  }
  invisible(x)
}

check_numeric <- function(x, what) {
  if (!is.numeric(x) || anyNA(x)) {
    stop(sprintf("`%s` must be numeric without missing values.", what),
         call. = FALSE)
  }
  invisible(x)
}

# Recycle a bound or value vector to length n; only length 1 or n allowed.
recycle_to <- function(x, n, what) {
  if (length(x) == 1L) return(rep.int(x, n))
  if (length(x) == n) return(as.vector(x))
  stop(sprintf("`%s` has length %d; it must have length 1 or %d.",
               what, length(x), n), call. = FALSE)
}

# Element labels for a block with optional dim and dimnames.
element_labels <- function(name, len, dim = NULL, dimnames = NULL) {
  if (is.null(dim) && len == 1L) return(name)
  if (is.null(dim)) dim <- len
  if (length(dim) == 1L) {
    lab <- if (!is.null(dimnames) && !is.null(dimnames[[1L]])) {
      dimnames[[1L]]
    } else {
      as.character(seq_len(dim))
    }
    return(sprintf("%s[%s]", name, lab))
  }
  idx <- lapply(seq_along(dim), function(k) {
    if (!is.null(dimnames) && !is.null(dimnames[[k]])) dimnames[[k]]
    else as.character(seq_len(dim[k]))
  })
  grid <- expand.grid(idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  sprintf("%s[%s]", name, do.call(paste, c(grid, sep = ",")))
}

fmt_num <- function(x) sprintf("%.6g", x)
