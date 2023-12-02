#include <cmath>       /* tgamma */
#include <map>

//#define _USE_MATH_DEFINES

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb);

const __int128_t one128 = 1;

/******************************************************************************/
/**************************   MODEL COMPLEXITY   ******************************/
/******************************************************************************/

/********* of a Sub-Complete Model based on m basis operators     *************/
/******************************************************************************/
// for 1 <= m <= n . Rem C(m=1) = log(pi)
double GeomComplexity_ICC(unsigned int m)     // Geometric complexity
{
  double pow = (double) ( one128 << (m-1) );
  return (log(M_PI) * pow - lgamma(pow)  );   // lgamma(x) = log(gamma(x))
}

double ParamComplexity_ICC(unsigned int m, unsigned int N)  // Parameter Complexity
{
  __uint32_t K = (one128 << m) - 1;  // number of interactions
  return K * log(((double) N)/2./M_PI) / 2.;
}

/********* of the Minimally Complex Model (MCM) defined by "Partition"   ******/
/******************************************************************************/
// Compute separately: -- the first order complexity    --> stored in C_param
//                     -- and the geometric complexity  --> stored in C_geom

double Complexity_MCM(std::map<unsigned int, __int128_t> Partition, unsigned int N, double *C_param, double *C_geom)
{
  *C_param = 0;   *C_geom = 0;
  unsigned int m_i = 0;  // number of elements in the i-th part "Ai"

  for (std::map<unsigned int, __int128_t>::iterator Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    m_i = Bitset_count((*Part).second);
    (*C_param) += ParamComplexity_ICC(m_i, N);
    (*C_geom) += GeomComplexity_ICC(m_i);
  }  

  return (*C_param) + (*C_geom);
}

