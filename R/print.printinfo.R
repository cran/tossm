`print.printinfo` <-
function( x, ...){
  x <- attr( x, 'printinfo')
  NextMethod()
}

