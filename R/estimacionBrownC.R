function(datos,col_zt,anio,mes,A,C,D,metodo=c("aditivo","multiplicativo"),prediccion=NULL){
  Ft_as=rep(0,nrow(datos))
  Zt_as=rep(0,nrow(datos))
  Tt_as=rep(0,nrow(datos))
  Zt_suav=rep(0,nrow(datos))
  at=rep(0,nrow(datos))
  media_anio1=mean(datos[1:12,col_zt])
  Zt_as[12]=media_anio1

  if (metodo=="multiplicativo"){
    for (i in 1:12) {
      Ft_as[i]=datos[i,col_zt]/media_anio1
    }
    for (j in 13:nrow(datos)) {
      Zt_as[j]=A*(datos[j,col_zt]/Ft_as[(j-12)])+(1-A)*(Zt_as[(j-1)]+Tt_as[(j-1)])
      Tt_as[j]=C*(Zt_as[j]-Zt_as[(j-1)])+(1-C)*Tt_as[(j-1)]
      Ft_as[j]=D*(datos[j,col_zt]/Zt_as[j])+(1-D)*Ft_as[(j-12)]
      Zt_suav[j]=(Zt_as[(j-1)]+1*Tt_as[(j-1)])*Ft_as[(j-12)]
      at[j]=datos[j,col_zt]-Zt_suav[j]
    }
  } else {
    for (i in 1:12) {
      Ft_as[i]=datos[i,col_zt]-media_anio1
    }
    for (j in 13:nrow(datos)) {
      Zt_as[j]=A*(datos[j,col_zt]-Ft_as[(j-12)])+(1-A)*(Zt_as[(j-1)]+Tt_as[(j-1)])
      Tt_as[j]=C*(Zt_as[j]-Zt_as[(j-1)])+(1-C)*Tt_as[(j-1)]
      Ft_as[j]=D*(datos[j,col_zt]-Zt_as[j])+(1-D)*Ft_as[(j-12)]
      Zt_suav[j]=(Zt_as[(j-1)]+1*Tt_as[(j-1)])+Ft_as[(j-12)]
      at[j]=datos[j,col_zt]-Zt_suav[j]
    }
  }

  Zt_as_p=Zt_as[nrow(datos)]
  Tt_as_p=Tt_as[nrow(datos)]
  Ft_as_p=tail(Ft_as,12)

  ###para prediccion
  if(!is.null(prediccion)){
    h=prediccion-nrow(datos)
    mes=prediccion %% 12
    for (i in 1:12) {
      if (mes==i){
        Ft_as_pr=Ft_as_p[i]
      }
    }
    if(metodo=="multiplicativo"){
      Estimado=(Zt_as_p + h*Tt_as_p)*Ft_as_pr
    } else {
      Estimado=(Zt_as_p + h*Tt_as_p)+Ft_as_pr
    }
    SC=sum(at^2)
    resul=data.frame(cbind(Ft_as,Zt_as,Tt_as,Zt_suav,at))
    list(HoltWinters_Multiplicativo=resul,
         Suma_Cuadrados_at=SC,
         Prediccion=Estimado)
  } else {
    SC=sum(at^2)
    resul=data.frame(cbind(Ft_as,Zt_as,Tt_as,Zt_suav,at))
    list(HoltWinters_Multiplicativo=resul,
         Suma_Cuadrados_at=SC)
  }
}
