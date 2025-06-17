#' all topological orderings of a DAG
#'
#' This function computes all topological orderings of a graph
#' using the recursive algorithm described in Knuth and Szwarcfiter (1974).
#'
#' @param n number of nodes in the DAG
#' @param adj_list edges given as an adjacency list
#'
#' @return Returns a list of topological orderings.
#'
#' @examples
#' n = 4
#' alltopsorts_recursion(n, list(c(1,3),c(2,4)))
#'
#'
#' @references Knuth, D. E. and J. L. Szwarcfiter (1974).
#' A structured program to generate all topological sorting arrangements.
#' Information Processing Letters 2(6), 153–157.
#'
#' @export
alltopsorts_recursion = function(n, adj_list){

  if((n%%1!=0)|(n<=1)){
    return('Error: n has to be an integer larger than 1')
  }

  if(!is.list(adj_list)){
    return("Error: adj_list has to be an adjacency list")
  }

  if(length(adj_list) == 0){
    Dlist = 1:n
  } else{
    countmat = sapply(1:n, function(s)
      sapply(adj_list, function(a) a[2]==s))
    if(is.matrix(countmat)){
      counts = colSums(countmat)
    } else{
      counts = countmat + 0
    }
    Dlist = which(counts==0)
  }
  PERMUTATIONRESULTS <- list(numeric(n))

  skipbeginning = FALSE
  decrease_rec_depth = FALSE
  alltopsorts_list = list(list(n=n,
                               adj_list=adj_list,
                               Dlist=Dlist,
                               deleted_vertices=c()))
  rec_depth = 1
  while(rec_depth >= 1){
    if(!skipbeginning){
      alltopsorts_list[[rec_depth]]$base =
        alltopsorts_list[[rec_depth]]$
        Dlist[length(alltopsorts_list[[rec_depth]]$Dlist)]
      alltopsorts_list[[rec_depth]]$repeatuntil = FALSE
      skipbeginning = TRUE
    }
    if(alltopsorts_list[[rec_depth]]$repeatuntil){
      decrease_rec_depth = TRUE
    }

    while(!alltopsorts_list[[rec_depth]]$repeatuntil){
      alltopsorts_list[[rec_depth]]$q =
        alltopsorts_list[[rec_depth]]$
        Dlist[length(alltopsorts_list[[rec_depth]]$Dlist)]
      alltopsorts_list[[rec_depth]]$deleted_vertices_new =
        c(alltopsorts_list[[rec_depth]]$
            deleted_vertices,alltopsorts_list[[rec_depth]]$q)
      if(length(alltopsorts_list[[rec_depth]]$adj_list)==0){
        alltopsorts_list[[rec_depth]]$counts =
          rep(0,alltopsorts_list[[rec_depth]]$n)
        alltopsorts_list[[rec_depth]]$adj_list_new =
          alltopsorts_list[[rec_depth]]$adj_list
      } else{
        alltopsorts_list[[rec_depth]]$adj_list_new =
          alltopsorts_list[[rec_depth]]$
          adj_list[!sapply(alltopsorts_list[[rec_depth]]$
                             adj_list,function(a) a[1]==
                             alltopsorts_list[[rec_depth]]$q)]
        if(length(alltopsorts_list[[rec_depth]]$adj_list_new)==0){
          alltopsorts_list[[rec_depth]]$counts = rep(0,n)
        } else{
          alltopsorts_list[[rec_depth]]$countmat =
            sapply(1:alltopsorts_list[[rec_depth]]$n, function(s)
              sapply(alltopsorts_list[[rec_depth]]$adj_list_new,
                     function(a) a[2]==s))
          if(is.matrix(alltopsorts_list[[rec_depth]]$countmat)){
            alltopsorts_list[[rec_depth]]$counts =
              colSums(alltopsorts_list[[rec_depth]]$countmat)
          } else{
            alltopsorts_list[[rec_depth]]$counts =
              alltopsorts_list[[rec_depth]]$countmat + 0
          }
        }
      }
      alltopsorts_list[[rec_depth]]$Dlist_new =
        ((1:alltopsorts_list[[rec_depth]]$n)*
           (alltopsorts_list[[rec_depth]]$
              counts==0))[-alltopsorts_list[[rec_depth]]$deleted_vertices_new]
      alltopsorts_list[[rec_depth]]$Dlist_new =
        alltopsorts_list[[rec_depth]]$
        Dlist_new[alltopsorts_list[[rec_depth]]$Dlist_new!=0]
      alltopsorts_list[[rec_depth]]$rescopy = PERMUTATIONRESULTS
      alltopsorts_list[[rec_depth]]$
        rescopy[[length(alltopsorts_list[[rec_depth]]$
                          rescopy)]][length(alltopsorts_list[[rec_depth]]$
                                              deleted_vertices_new)] <-
        alltopsorts_list[[rec_depth]]$q
      PERMUTATIONRESULTS <- alltopsorts_list[[rec_depth]]$rescopy
      if(length(alltopsorts_list[[rec_depth]]$deleted_vertices_new) ==
         alltopsorts_list[[rec_depth]]$n-1){
        alltopsorts_list[[rec_depth]]$rescopy = PERMUTATIONRESULTS
        alltopsorts_list[[rec_depth]]$
          rescopy[[length(alltopsorts_list[[rec_depth]]$
                            rescopy)]][length(alltopsorts_list[[rec_depth]]$
                                                deleted_vertices_new)] <-
          alltopsorts_list[[rec_depth]]$q
        PERMUTATIONRESULTS <- alltopsorts_list[[rec_depth]]$rescopy
        alltopsorts_list[[rec_depth]]$recopy = PERMUTATIONRESULTS
        alltopsorts_list[[rec_depth]]$
          rescopy[[length(alltopsorts_list[[rec_depth]]$
                            rescopy)]][alltopsorts_list[[rec_depth]]$n] <-
          (1:alltopsorts_list[[rec_depth]]$n)[-alltopsorts_list[[rec_depth]]$
                                                deleted_vertices_new]
        alltopsorts_list[[rec_depth]]$
          rescopy[[length(alltopsorts_list[[rec_depth]]$rescopy)+1]] =
          numeric(alltopsorts_list[[rec_depth]]$n)
        PERMUTATIONRESULTS <- alltopsorts_list[[rec_depth]]$rescopy
        recursion = FALSE
      } else{
        recursion = TRUE
        skipbeginning = FALSE
        alltopsorts_list[[rec_depth + 1]] =
          list(n = alltopsorts_list[[rec_depth]]$n,
               adj_list = alltopsorts_list[[rec_depth]]$adj_list_new,
               Dlist = alltopsorts_list[[rec_depth]]$Dlist_new,
               deleted_vertices = alltopsorts_list[[rec_depth]]$
                 deleted_vertices_new)
      }
      alltopsorts_list[[rec_depth]]$Dlist =
        c(alltopsorts_list[[rec_depth]]$q,
          alltopsorts_list[[rec_depth]]$
            Dlist[-length(alltopsorts_list[[rec_depth]]$Dlist)])
      alltopsorts_list[[rec_depth]]$repeatuntil =
        (alltopsorts_list[[rec_depth]]$base ==
           alltopsorts_list[[rec_depth]]$
           Dlist[length(alltopsorts_list[[rec_depth]]$Dlist)])
      if(recursion){
        rec_depth = rec_depth + 1
        alltopsorts_list[[rec_depth]]$repeatuntil = TRUE
      } else{
        if(alltopsorts_list[[rec_depth]]$repeatuntil){
          decrease_rec_depth = TRUE
        }
      }
    }
    if(decrease_rec_depth){
      rec_depth = rec_depth-1
      decrease_rec_depth = FALSE
    }
  }

  PERMUTATIONRESULTS = PERMUTATIONRESULTS[-length(PERMUTATIONRESULTS)]

  for(s in seq_along(PERMUTATIONRESULTS)[-1]){
    PERMUTATIONRESULTS[[s]][PERMUTATIONRESULTS[[s]]==0] =
      PERMUTATIONRESULTS[[s-1]][PERMUTATIONRESULTS[[s]]==0]
  }

  return(PERMUTATIONRESULTS)
}
