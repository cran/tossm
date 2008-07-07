`calc.effort` <-
function( catch.locs,yr){
  #effort is average distance from shore of harvested animals
  if (!is.null(catch.locs[[yr]])) {mean(catch.locs[[yr]][,2])  
  } else 0/0
}

