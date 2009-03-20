	
dyn.load("./libtwovarcomp.so.0.01")

FastMixedModel<-function(Response, Explan, Kin, Covariates=NULL, nu_naught=10^(-5), gamma_naught=10^5)
{
 Response=as.matrix(Response)
 Explan=as.matrix(Explan)
 n=length(Response)
 if(is.null(Covariates))
    Covariates=matrix(1,n,1)
 
 Covariates=as.matrix(Covariates)

 ret=.Call("rint_flmm", as.double(Explan), as.double(Response), as.integer(dim(Explan)[1]), as.integer(dim(Explan)[2]), as.double(t(Covariates)), as.integer(dim(Covariates)[2]), as.double(t(Kin)), as.double(nu_naught), as.double(gamma_naught))
  ret
}

Result=FastMixedModel(Response, Explan, Kin)

Result$chi.sq


