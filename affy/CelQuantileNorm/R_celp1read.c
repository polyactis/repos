/**
 * Copyright (c) 2006 Hin-Tak Leung
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 **/

/* 
   This code is endian specific. Please don't expect it to work on
   anything other than i686/x86_64 unices.
*/

#include <R.h>
#include <Rinternals.h>

#define CEL_INTENSITY_ARRAY_SIZE 6553600

SEXP R_celp1read(SEXP filename)
{
  FILE *fhandle = NULL;
  SEXP ans = R_NilValue;
  PROTECT(filename = coerceVector(filename, STRSXP));
  PROTECT(ans = allocVector(REALSXP, CEL_INTENSITY_ARRAY_SIZE));
  
  fhandle = fopen(CHAR(STRING_ELT(filename, 0)), "rb");
  size_t readsize = fread(REAL(ans), sizeof(double), CEL_INTENSITY_ARRAY_SIZE, fhandle);
  Rprintf("Read %d of %d\n", readsize, CEL_INTENSITY_ARRAY_SIZE);
  fclose(fhandle);
  
  UNPROTECT(2);
  return ans;
}
