"calc.effort" <-
function( catch.locs,yr){
  #effort is average distance from shore of harvested animals
  if (!is.null(catch.locs[[yr]])) {mean(catch.locs[[yr]][,3])  
  } else 0/0
}


