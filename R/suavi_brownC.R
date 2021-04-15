function(datos,alfa){
  Zt_a=matrix(datos[,2],nrow(datos),length(alfa))
  Zt_aa=matrix(datos[,2],nrow(datos),length(alfa))
  Zt_aaa=matrix(datos[,2],nrow(datos),length(alfa))
  b1t=matrix(datos[,2],nrow(datos),length(alfa))
  b2t=matrix(0,nrow(datos),length(alfa))
  b3t=matrix(0,nrow(datos),length(alfa))
  Zt_sua=matrix(0,nrow(datos),length(alfa))
  at=matrix(0,nrow(datos),length(alfa))
  for (i in 2:nrow(datos)) {
    for (j in 1:length(alfa)) {
      Zt_a[i,j]=datos[i,2]*alfa[j] + Zt_a[(i-1),j]*(1-alfa[j])
      Zt_aa[i,j]=Zt_a[i,j]*alfa[j] + Zt_aa[(i-1),j]*(1-alfa[j])
      Zt_aaa[i,j]=Zt_aa[i,j]*alfa[j] + Zt_aaa[(i-1),j]*(1-alfa[j])
      b1t[i,j]=3*Zt_a[i,j]-3*Zt_aa[i,j]+Zt_aaa[i,j]
      b2t[i,j]=(alfa[j]/(2*(1-alfa[j])^2))*((6-5*alfa[j])*Zt_a[i,j]-
                                              2*(5-4*alfa[j])*Zt_aa[i,j]+
                                              (4-3*alfa[j])*Zt_aaa[i,j])
      b3t[i,j]=((alfa[j]/(1-alfa[j]))^2)*(Zt_a[i,j]-
                                            2*Zt_aa[i,j]+
                                            Zt_aaa[i,j])
      Zt_sua[i,j]=b1t[(i-1),j]+1*b2t[(i-1),j]+(1^2)*b3t[(i-1),j]
      at[i,j]=datos[i,2]-Zt_sua[i,j]
    }
  }
  SC=rep(0,length(alfa))
  for (k in 1:length(alfa)) {
    SC[k]=sum(at[,k]^2)
  }
  list(ValoresSuavizados_Zt=Zt_sua,
       Valores_B1t=b1t,
       Valores_B2t=b2t,
       Valores_B3t=b3t,
       Valores_at=at,
       Suma_Cuadrados=SC)
}

