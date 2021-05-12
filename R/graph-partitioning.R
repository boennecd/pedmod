.prep_edge_list <- function(from, to, weight_data, edge_weights,
                            init = integer()){
  if(is.null(edge_weights))
    edge_weights <- rep(1, length(from))

  stopifnot(length(from) == length(to), length(from) == length(edge_weights),
            all(is.finite(from)), all(is.finite(to)),
            all(is.finite(edge_weights)))
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

  init <- match(init, id, nomatch = NA_integer_)
  init <- init[!is.na(init)] - 1L

  list(from = from, to = to, id = id, weights_ids = weights_ids,
       weights = weights, edge_weights = edge_weights, init = init)
}

#' Finds the Biconnected Components
#'
#' @description
#' Finds the biconnected components and the cut vertices (articulation points)
#' using the methods suggested by Hopcroft et al. (1973).
#'
#' @param from integer vector with one of the vertex ids.
#' @param to integer vector with one of the vertex ids.
#'
#' @seealso
#' \code{\link{get_block_cut_tree}} and
#' \code{\link{get_max_balanced_partition}}.
#'
#' @references
#' Hopcroft, J., & Tarjan, R. (1973).
#' \emph{Algorithm 447: efficient algorithms for graph manipulation}.
#' Communications of the ACM, 16(6), 372-378.
#'
#' @export
get_biconnected_components <- function(from, to){
  dat <- .prep_edge_list(from = from, to = to, weight_data = NULL,
                         edge_weights = NULL)
  id <- dat$id

  out <- with(dat, .get_biconnected_components(from, to, weights_ids,
                                               weights, edge_weights))
  lapply(out, function(x)
    structure(sort(id[x]),
              cut_verices = sort(id[attr(x, "cut_verices")])))
}

.pedigree_to_from_to <- function(id, father.id, mother.id, id_weight = NULL,
                                 father_weight = NULL, mother_weight = NULL){
  n <- length(id)
  stopifnot(length(father.id) == n, length(mother.id) == n)
  to <- c(id, id)
  from <- c(father.id, mother.id)
  keep <- !is.na(from) & from %in% id

  if(is.null(id_weight)){
    weights_ids = integer()
    weights = numeric()
  } else {
    stopifnot(length(id_weight) == n, all(is.finite(id_weight)))
    weights_ids <- id
    weights <- id_weight
  }

  if(is.null(father_weight))
    father_weight <- rep(1, n)
  else
    stopifnot(length(father_weight) == n, all(is.finite(father_weight)))
  if(is.null(mother_weight))
    mother_weight <- rep(1, n)
  else
    stopifnot(length(mother_weight) == n, all(is.finite(mother_weight)))
  edge_weights <- c(father_weight, mother_weight)

  list(to = to[keep], from = from[keep], edge_weights = edge_weights[keep],
       weight_data = list(id = weights_ids, weight = weights))
}

#' @rdname get_biconnected_components
#'
#' @param id integer vector with the child id.
#' @param father.id integer vector with the father id. May be \code{NA} if
#' it is missing.
#' @param mother.id integer vector with the mother id. May be \code{NA} if
#' it is missing.
#'
#' @export
get_biconnected_components_pedigree <- function(id, father.id, mother.id){
  dat <- .pedigree_to_from_to(id = id, father.id = father.id,
                              mother.id = mother.id)

  with(dat, get_biconnected_components(from = from, to = to))
}

#' Creates a Block-cut Tree Like Object
#'
#' @description
#' Creates a block-cut tree like structure computed using the method suggested
#' by Hopcroft et al. (1973).
#'
#' @inheritParams get_biconnected_components
#'
#' @seealso
#' \code{\link{get_biconnected_components}} and
#' \code{\link{get_max_balanced_partition}}.
#'
#' @references
#' Hopcroft, J., & Tarjan, R. (1973).
#' \emph{Algorithm 447: efficient algorithms for graph manipulation}.
#' Communications of the ACM, 16(6), 372-378.
#'
#' @export
get_block_cut_tree <- function(from, to){
  dat <- .prep_edge_list(from = from, to = to, weight_data = NULL,
                         edge_weights = NULL)
  id <- dat$id

  out <- with(dat, .get_block_cut_tree(from, to, weights_ids,
                                       weights, edge_weights))
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

.sort_partition_result <- function(out, id){
  removed_edges <- out$removed_edges
  removed_edges[] <- id[removed_edges]
  removed_edges[] <- t(apply(removed_edges, 1L, sort))
  removed_edges <- removed_edges[
    order(removed_edges[, 1], removed_edges[, 2]), , drop = FALSE]
  out$removed_edges <- removed_edges

  out$set_1 <- sort(id[out$set_1])
  out$set_2 <- sort(id[out$set_2])

  out
}

#' Finds an Approximately Balanced Connected Partition
#' @description
#' Uses the method suggested by Chlebíková (1996) to construct an approximate
#' maximally balanced connected partition. A further refinement step can be made
#' to reduce the cost of the cut edges. See
#' \code{vignette("pedigree_partitioning", package = "pedmod")} for further
#' details.
#'
#' @inheritParams get_biconnected_components
#' @param weight_data list with two elements called \code{"id"} for the id and
#' \code{"weight"} for the vertex weight. All vertices that are not in this list
#' have a weight of one. Use \code{NULL} if all vertices have a weight of one.
#' @param slack fraction between zero and 0.5 for the allowed amount of
#' deviation from the balance criterion that is allowed to reduce the cost of
#' the cut edges.
#' @param max_kl_it_inner maximum number of moves to consider in each
#' iteration when \code{slack > 0}.
#' @param max_kl_it maximum number of iterations to use when reducing the
#' cost of the cut edges. Typically the method converges quickly and this
#' argument is not needed.
#' @param trace integer where larger values yields more information printed to
#' the console during the procedure.
#' @param edge_weights numeric vector with weights for each edge. Needs to have
#' the same length as \code{from} and \code{to}. Use \code{NULL} if all edges
#' should have a weight of one.
#' @param check_weights logical for whether to check the weights in each
#' biconnected component. This may fail if the graph is not connected in which
#' case the results will likely be wrong. It may also fail for large graphs
#' because of floating-point arithmetic. The latter is not an error and the
#' reason for this argument.
#' @param do_reorder logical for whether the implementation should reorder the
#' vertices. This may reduce the computation time for some data sets.
#'
#' @seealso
#' \code{\link{get_biconnected_components}} and
#' \code{\link{get_block_cut_tree}}.
#'
#' @references
#' Chlebíková, J. (1996).
#' \emph{Approximating the maximally balanced connected partition problem in
#' graphs}. Information Processing Letters, 60(5), 225-230.
#'
#' Hopcroft, J., & Tarjan, R. (1973).
#' \emph{Algorithm 447: efficient algorithms for graph manipulation}.
#' Communications of the ACM, 16(6), 372-378.
#'
#' @export
get_max_balanced_partition <- function(from, to, weight_data = NULL,
                                       edge_weights = NULL,
                                       slack = 0., max_kl_it_inner = 50L,
                                       max_kl_it = 10000L, trace = 0L,
                                       check_weights = TRUE,
                                       do_reorder = FALSE){
  dat <- .prep_edge_list(from = from, to = to, weight_data = weight_data,
                         edge_weights = edge_weights)
  id <- dat$id

  out <- with(dat, .get_max_balanced_partition(
    from, to, weights_ids,  weights, slack = slack, edge_weights = edge_weights,
    max_kl_it_inner = max_kl_it_inner, max_kl_it = max_kl_it, trace = trace,
    check_weights = check_weights, do_reorder = do_reorder))

  # set and sort the removed edges
  .sort_partition_result(out, id)
}

#' @rdname get_max_balanced_partition
#'
#' @param id_weight numeric vector with the weight to use for each vertex
#' (individual). \code{NULL} yields a weight of one for all.
#' @param father_weight weights of the edges created between the fathers
#' and the children. Use \code{NULL} if all should have a weight of one.
#' @param mother_weight weights of the edges created between the mothers
#' and the children. Use \code{NULL} if all should have a weight of one.
#'
#' @export
get_max_balanced_partition_pedigree <- function(
  id, father.id, mother.id, id_weight = NULL, father_weight = NULL,
  mother_weight = NULL, slack = 0., max_kl_it_inner = 50L, max_kl_it = 10000L,
  trace = 0L, check_weights = TRUE, do_reorder = FALSE){
  dat <- .pedigree_to_from_to(id = id, father.id = father.id,
                              mother.id = mother.id, id_weight = id_weight,
                              father_weight = father_weight,
                              mother_weight = mother_weight)

  with(dat, get_max_balanced_partition(
    from = from, to = to, weight_data = weight_data, slack = slack,
    max_kl_it_inner = max_kl_it_inner, max_kl_it = max_kl_it, trace = trace,
    edge_weights = edge_weights, do_reorder = do_reorder,
    check_weights = check_weights))
}

#' Finds an Approximately Balanced Partition
#' @inheritParams get_max_balanced_partition
#' @export
get_unconnected_partition <- function(from, to, weight_data = NULL,
                                      edge_weights = NULL,
                                      slack = 0., max_kl_it_inner = 50L,
                                      max_kl_it = 10000L, trace = 0L,
                                      init = integer()){
  dat <- .prep_edge_list(from = from, to = to, weight_data = weight_data,
                         edge_weights = edge_weights, init = init)
  id <- dat$id

  out <- with(dat, .get_unconnected_partition(
    from, to, weights_ids,  weights, slack = slack, edge_weights = edge_weights,
    max_kl_it_inner = max_kl_it_inner, max_kl_it = max_kl_it, trace = trace,
    init = init))

  # set and sort the removed edges
  .sort_partition_result(out, id)
}

#' @rdname get_unconnected_partition
#' @export
get_unconnected_partition_pedigree <- function(
  id, father.id, mother.id, id_weight = NULL, father_weight = NULL,
  mother_weight = NULL, slack = 0., max_kl_it_inner = 50L, max_kl_it = 10000L,
  trace = 0L){
  dat <- .pedigree_to_from_to(id = id, father.id = father.id,
                              mother.id = mother.id, id_weight = id_weight,
                              father_weight = father_weight,
                              mother_weight = mother_weight)

  with(dat, get_unconnected_partition(
    from = from, to = to, weight_data = weight_data, slack = slack,
    max_kl_it_inner = max_kl_it_inner, max_kl_it = max_kl_it, trace = trace,
    edge_weights = edge_weights))
}
