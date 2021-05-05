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
#'
#' @param id_weight numeric vector with the weight to use for each vertex
#' (individual). \code{NULL} yields a weight of one for all.
#'
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
