#  James Rekow

firstDeriv = function(y, deltax = 1){
  
  #  ARGS:  y - numeric vector of y-values
  #         deltax - mesh width of corresponding x-values
  #
  #  RETURNS:  a numeric vector of the same length as y whose elements are the
  #            centered first derivative of y at each interior point, and one-sided
  #            derivatives of y at the endpoints
  
  N = length(y)
  return((0.5 / deltax) * c(2 * (y[2] - y[1]), (y[3:N] - y[1:(N - 2)]), 2 * (y[N] - y[N - 1])))
  
} #  end firstDeriv function
