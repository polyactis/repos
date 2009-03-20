#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "twovarcomp.h"

#ifdef _OPENMP
#include <omp.h>
#define OMP_GET_MAX_THREADS omp_get_max_threads()
#else
#define OMP_GET_MAX_THREADS 1
#endif // _OPENMP

// for debug
#include "moorepenrose.h"
extern "C" {

#include "main.h"


SEXP rint_flmm(SEXP pexplan_sexp, SEXP presp_sexp, SEXP pn_sexp, SEXP pp_sexp, SEXP pcovar_sexp, SEXP pp_covar_sexp, SEXP pVar2_sexp, SEXP nu_naught_sexp, SEXP gamma_naught_sexp)
{

  double *pexplan, *presp,  *pnu_naught, *pgamma_naught, *pcovar, *pVar2;
  double* pchisq;
  double* pherit;
  unsigned int *pn, *pp_covar, *pp;
 

  char* pret_names[2]={"chi.sq", "herit"};

  SEXP preturn_list_SEXP, preturn_names_SEXP, paname_SEXP;
    
  SEXP pchisq_SEXP;
  SEXP pherit_SEXP;

  gsl_matrix* pvar1_mat, *pvar2_mat;
  gsl_matrix* pcovar_mat;
 
  gsl_vector* presponse_vec;
      
  // really must check all gsl returns
  //gsl_set_error_handler_off();

  // C side pointers to R objects
  pexplan=(double*) REAL(pexplan_sexp);
  presp=(double*) REAL(presp_sexp);
  pn=(unsigned int*) INTEGER(pn_sexp);
  pp=(unsigned int*) INTEGER(pp_sexp);
  pcovar=(double*) REAL(pcovar_sexp);
  pp_covar=(unsigned int*) INTEGER(pp_covar_sexp);
  pVar2=(double*) REAL(pVar2_sexp);
  pnu_naught=(double*)REAL(nu_naught_sexp);
  pgamma_naught=(double*)REAL(gamma_naught_sexp);
  

  presponse_vec=&(gsl_vector_view_array(presp, *pn).vector);
  pvar2_mat=&(gsl_matrix_view_array(pVar2, *pn, *pn).matrix);
  pvar1_mat=gsl_matrix_alloc(*pn, *pn);
  pcovar_mat=&(gsl_matrix_view_array(pcovar, *pn, *pp_covar).matrix);
  
  gsl_matrix* pincid1_mat, *pincid2_mat;
 
  pincid1_mat=gsl_matrix_alloc(*pn,*pn);
  pincid2_mat=gsl_matrix_alloc(*pn,*pn);

  gsl_matrix_set_identity(pvar1_mat);
  // gsl_matrix_set_identity(pvar2_mat); 
  gsl_matrix_set_identity(pincid1_mat);
  gsl_matrix_set_identity(pincid2_mat); 

   
  PROTECT(pchisq_SEXP=NEW_NUMERIC(*pp));
  pchisq=NUMERIC_POINTER(pchisq_SEXP);
  PROTECT(pherit_SEXP=NEW_NUMERIC(*pp));
  pherit=NUMERIC_POINTER(pherit_SEXP);

    
  TwoVarCompModel DaddyTwoVarCompModel(presponse_vec, pcovar_mat, pvar1_mat, pvar2_mat, pincid1_mat, pincid2_mat);  
  double nullminimand=0.5;
  double altminimand;
  double nulldev=DaddyTwoVarCompModel.MinimiseNullDeviance(&nullminimand);
  //std::cout << "null sigmasq=" << nullminimand<<std::endl<<std::endl;
  //DaddyTwoVarCompModel.NullDeviance(0.5);
  
#pragma omp parallel shared(pexplan, pp, pn, pchisq, pherit, nulldev, nullminimand) private(altminimand)
    {
      TwoVarCompModel ChildTwoVarCompModel(DaddyTwoVarCompModel);
      
#pragma omp for
      for(int it=0;it<*pp;it++)
	{
#pragma omp critical
	  { 
	    std::cout<<".";
	  }
	  TwoVarCompModel ATwoVarCompModel(DaddyTwoVarCompModel);
	  ChildTwoVarCompModel.SetExplan(&(gsl_vector_view_array(pexplan+(*pn)*it, *pn).vector));
	  altminimand=nullminimand;
	  pchisq[it]=nulldev-ChildTwoVarCompModel.MinimiseDeviance(&altminimand);
	
	}    
    }
  
  
  PROTECT(preturn_list_SEXP=allocVector(VECSXP,1));
  SET_VECTOR_ELT(preturn_list_SEXP, 0,pchisq_SEXP);
  //SET_VECTOR_ELT(preturn_list_SEXP, 1,pherit_SEXP);
 
 
  PROTECT(preturn_names_SEXP=allocVector(STRSXP,1));

  
  for(int it=0;it<1;it++)
    {
      std::cout<<pret_names[it];
      //paname_SEXP=Rf_mkChar(pret_names[it]);
      SET_STRING_ELT(preturn_names_SEXP,it,Rf_mkChar(pret_names[it]));
    }
  setAttrib(preturn_list_SEXP, R_NamesSymbol,preturn_names_SEXP);
  
  UNPROTECT(4);

  return preturn_list_SEXP;
}

}

