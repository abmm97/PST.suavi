#' @name estimacionBrownC
#' @title Quadratic Brown Estimation
#'
#' @param datos Array of data with columns t, Zt
#' @param posicion_estimada estimation period
#' @param alfa parameter value <1
#'
#' @return returns the value of the smoothed Zt of the requested period
#' @export estimacionBrownC
estimacionBrownC<-function(datos,posicion_estimada,alfa){
  b1t=suavi_brownC(datos,alfa)$Valores_B1t[nrow(datos)]
  b2t=suavi_brownC(datos,alfa)$Valores_B2t[nrow(datos)]
  b3t=suavi_brownC(datos,alfa)$Valores_B3t[nrow(datos)]
  h=posicion_estimada-nrow(datos)
  Estimado=b1t+b2t*h+b3t*h^2
  list(Estimado=Estimado)
}
