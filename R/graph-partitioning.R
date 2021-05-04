.prep_edge_list <- function(from, to, weight_data){
  stopifnot(length(from) == length(to))
  id <- unique(c(from, to))
  n <- length(id)
  from <- match(from, id, nomatch = n + 1L) - 1L
  to   <- match(to  , id, nomatch = n + 1L) - 1L

  if(is.null(weight_data)){
    weights_ids = integer()
    weights = numeric()
  } else {
    # checks
    stopifnot(is.list(weight_data),
              c("id", "weight") %in% names(weight_data),
              length(weight_data$id) == length(weight_data$weight),
              is.vector(weight_data$id), is.finite(weight_data$id),
              is.vector(weight_data$weight), is.finite(weight_data$weight))

    # get the elements
    weights_ids <- match(weight_data$id, id, nomatch = n  + 1L) - 1L
    weights <- weight_data$weight
  }

  list(from = from, to = to, id = id, weights_ids = weights_ids,
       weights = weights)
}

#' Finds the Biconnected Components Given
#' @export
get_biconnected_components <- function(from, to){
  dat <- .prep_edge_list(from = from, to = to, weight_data = NULL)
  id <- dat$id

  out <- with(dat, .get_biconnected_components(from, to, weights_ids,
                                               weights))
  lapply(out, function(x)
    structure(sort(id[x]),
              cut_verices = sort(id[attr(x, "cut_verices")])))
}

.pedigree_to_from_to <- function(id, father.id, mother.id, id_weight = NULL){
  n <- length(id)
  stopifnot(length(father.id) == n, length(mother.id) == n)
  to <- c(id, id)
  from <- c(father.id, mother.id)
  keep <- !is.na(from) & from %in% id

  if(is.null(id_weight)){
    weights_ids = integer()
    weights = numeric()
  } else {
    stopifnot(length(id_weight) == n)
    weights_ids <- id
    weights <- id_weight
  }

  list(to = to[keep], from = from[keep],
       weight_data = list(id = weights_ids, weight = weights))
}

#' @rdname get_biconnected_components
#' @export
get_biconnected_components_pedigree <- function(id, father.id, mother.id){
  dat <- .pedigree_to_from_to(id = id, father.id = father.id,
                              mother.id = mother.id)

  with(dat, get_biconnected_components(from = from, to = to))
}

#' Creates a Block Cut Tree Like Object
#' @export
get_block_cut_tree <- function(from, to){
  dat <- .prep_edge_list(from = from, to = to, weight_data = NULL)
  id <- dat$id

  out <- with(dat, .get_block_cut_tree(from, to, weights_ids,
                                       weights))
  . <- function(x)
    list(vertices     = sort(id[x$vertices]),
         cut_vertices = sort(id[x$cut_vertices]),
         leafs        = lapply(x$leafs, .))
  .(out)
}

#' @rdname get_block_cut_tree
#' @export
get_block_cut_tree_pedigree <- function(id, father.id, mother.id){
  dat <- .pedigree_to_from_to(id = id, father.id = father.id,
                              mother.id = mother.id)

  with(dat, get_block_cut_tree(from = from, to = to))
}

#' Finds an Approximately Balanced Connected Partition
#' @description
#' Uses the method suggested by Chlebíková (1996) to construct an approximate
#' maximally balanced connected partition. A furhter refinement step can be made
#' to reduce the cost of the cut edges.
#' @export
get_max_balanced_partition <- function(from, to, weight_data = NULL,
                                       slack = 0., max_kl_it_inner = 50L,
                                       max_kl_it = 10000L, trace = 0L){
  dat <- .prep_edge_list(from = from, to = to, weight_data = weight_data)
  id <- dat$id

  out <- with(dat, .get_max_balanced_partition(
    from, to, weights_ids,  weights, slack = slack,
    max_kl_it_inner = max_kl_it_inner, max_kl_it = max_kl_it, trace = trace))
  out$removed_edges[] <- id[out$removed_edges]
  out$set_1 <- sort(id[out$set_1])
  out$set_2 <- sort(id[out$set_2])
  out
}

#' @rdname get_max_balanced_partition
#' @export
get_max_balanced_partition_pedigree <- function(id, father.id, mother.id,
                                                id_weight = NULL, slack = 0.,
                                                max_kl_it_inner = 50L,
                                                max_kl_it = 10000L,
                                                trace = 0L){
  dat <- .pedigree_to_from_to(id = id, father.id = father.id,
                              mother.id = mother.id, id_weight = id_weight)

  with(dat, get_max_balanced_partition(
    from = from, to = to, weight_data = weight_data, slack = slack,
    max_kl_it_inner = max_kl_it_inner, max_kl_it = max_kl_it, trace = trace))
}
