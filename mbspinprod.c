#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mcbsp.h"
#include "mcbsp-affinity.h"
#include "mbspdebug.h"
#include "mbsptime.h"

#include "mbsputil.h"
#include "mbsp-discover.h"

#define SZDBL (sizeof(double))
//N0=1; 
int N=4000; 

unsigned long long clvl1_comm[10];
unsigned long long clvl2_comm[10];
unsigned long long clvl1_sync[10];
unsigned long long clvl2_sync[10];


// start: initizing the counters
void clean_counters() {
  int i;
  for (i=0; i<10; i++) {
    clvl1_comm[i] = 0;
    clvl2_comm[i] = 0;
    clvl1_sync[i] = 0;
    clvl1_sync[i] = 0;
  }
}

void print_counters(int vsize ) {
  int val1_comm, val1_sync;
  int val2_comm, val2_sync, i;
  
  val2_comm = 0;
  val2_sync = 0;
  for (i=0; i<10; i++) {
    if (clvl2_comm[i] != 0) val2_comm += clvl2_comm[i];
    if (clvl2_sync[i] != 0) val2_sync += clvl2_sync[i];

  }

  val1_comm = 0;
  val1_sync = 0;

  for (i=0; i<10; i++) {
    if (clvl1_comm[i] != 0) val1_comm += clvl1_comm[i];
    if (clvl1_sync[i] != 0) val1_sync += clvl1_sync[i];

  }
  printf ("%d %llu %llu %llu %llu #counters \n", vsize, val1_comm, val1_sync, val2_comm, val2_sync);

}
// end: initizing the counters

// root from multibsp tree
multibsp_tree_node_t tree_root = NULL;


/* 
 * Compute number of local components of procesor s for vector of length n
 * distributed cyclically over p processors.
 */
int nloc (int p, int s, int n) 
{
  return (n+p-s-1)/p;
}


/*
 * Compute inner product of vectors x and y of length n>=0
 * << Secuencial code >>
 */
double bspip(double* x, double* y, int n) 
{
  /* Compute inner product of vectors x and y of length n>=0 */

  int nloc(int p, int s, int n);
  double inprod, *Inprod, alpha;
  int i, t;

  inprod= 0.0;
  for (i=0; i<n; i++) inprod += x[i]*y[i];

  return inprod;
}


/*
 * ask for memory with the new size 
 * and make a copy of values 
 */
double* get_own_memory(double* v, int p, int s, int n) 
{
  int nl, iglob, i;
  double* x;

  nl= nloc(p,s,n);
  x= vecallocd(nl);
  for (i=0; i<nl; i++){
    iglob= i*p+s;
    x[i]= v[iglob];
  }
  return x;
}


void multibsp_level1_inprod()
{
  double r1, inprod, alpha;
  double* vn2; 
  double* vown1; 
  double* r1_per_threads;
  double* Inprod;
  int i, t;
  int nl2, nl1;

  int N1=0;
  size_t pin1[100];
  multibsp_tree_node_t tnode     = NULL;

  tnode = tree_root;
  tnode = tnode->sons[bsp_pid()];
  N1 = tnode->length;
  for(i=0;i<N1; i++) pin1[i] = tnode->sons[i]->index;

  //printf("N1: %d, pid: %d\n", N1, bsp_pid());
  mcbsp_set_pinning( pin1, N1 );
  bsp_begin(N1);
  // getting n2 from up level
  
  bsp_push_reg(&nl2, SZINT);
  bsp_sync();  
  clvl1_sync[bsp_pid()] ++;
  

  bsp_sync();  
  clvl1_sync[bsp_pid()] ++;
  // getting the vector from up level
  vn2 = vecallocd(SZDBL*nl2);
  bsp_push_reg(vn2, SZDBL*nl2);
  bsp_sync();  
  clvl1_sync[bsp_pid()] ++;
  
  
  bsp_sync();  
  clvl1_sync[bsp_pid()] ++;

  nl1   = nloc(bsp_nprocs(), bsp_pid(), nl2);
  //log_info("nl1 = %d\n", nl1);
  vown1 = get_own_memory(vn2, bsp_nprocs(), bsp_pid(), nl2);
  r1 = bspip(vown1, vown1, nl1);


  bsp_sync(); 
  clvl1_sync[bsp_pid()] ++;


  r1_per_threads = vecallocd(bsp_nprocs()); 
  r1_per_threads[bsp_pid()] = r1;
  bsp_push_reg(r1_per_threads, SZDBL*bsp_nprocs());


  bsp_sync(); 
  clvl1_sync[bsp_pid()] ++;

  bsp_put(0, &r1, r1_per_threads, bsp_pid()*SZDBL, SZDBL);
  clvl1_comm[bsp_pid()] += SZDBL;
  bsp_sync(); 
  clvl1_sync[bsp_pid()] ++;

  bsp_end();  

}


double multibsp_level2_inprod(double* v, int n)
{
  double result;
  double* vown2; 
  double* level1_results;
  int i, nl2, total;

  int N2=0;
  size_t pin2[100];
  int N1=0;
  size_t pin1[100];


  multibsp_tree_node_t tnode     = NULL;


  // start: Multibsp level 2
  tnode = tree_root;
  N2 = tnode->length;
  //printf("pin2: [");
  for(i=0;i<N2; i++)  {  
    pin2[i] = tnode->sons[i]->index;
//  printf(" %d ", pin2[i]); 
  }
 // printf("]\n");

  mcbsp_set_pinning( pin2, N2 );
  bsp_begin(N2);
  nl2   = nloc(bsp_nprocs(), bsp_pid(), n);
  vown2 = get_own_memory(v, bsp_nprocs(), bsp_pid(), n);
  //for(i=0; i<nl2; i++) log_info("%d, vown2[%d] = %f\n",bsp_pid(), i, vown2[i]);


  // start: segundo nivel
  bsp_init(&multibsp_level1_inprod, 0, NULL);

  // start: Multibsp level 1
  tnode = tree_root;

  tnode = tnode->sons[bsp_pid()];
  N1 = tnode->length;
  //printf("pid: %d, n2: %d, n1: %d\n", bsp_pid(), N2, N1);

  //printf("pid: %d pin1: [", bsp_pid());
  for(i=0;i<N1; i++) {
     pin1[i] = tnode->sons[i]->index;
     //printf(" %d ", pin1[i]); 

  }
  //printf("]\n");
  mcbsp_set_pinning( pin1, N1 );

  bsp_begin(N1);

  // input params for the next level
  bsp_push_reg(&nl2 , SZINT     );
  bsp_sync(); 
  clvl1_sync[bsp_pid()]++;

  for (i=0; i<bsp_nprocs(); i++) bsp_put(i, &nl2, &nl2, 0, SZINT);
  bsp_sync(); 
  clvl1_sync[bsp_pid()]++;
  clvl1_comm[bsp_pid()]+= bsp_nprocs() * SZINT;
  
  
  bsp_push_reg(vown2, SZDBL*nl2 );
  bsp_sync();  
  clvl1_sync[bsp_pid()]++;
 
  for (i=0; i<bsp_nprocs(); i++) bsp_put(i, vown2, vown2, 0, SZDBL*nl2);
  bsp_sync(); 
  clvl1_sync[bsp_pid()]++;
  clvl1_comm[bsp_pid()]+= bsp_nprocs()*SZDBL*nl2;
  // waiting for read the params in the next level
  int nl1 = nloc(bsp_nprocs(), bsp_pid(), nl2); 
  double* vown1 = get_own_memory(vown2, bsp_nprocs(), bsp_pid(), nl2);
  
  double r1 = bspip(vown1, vown1, nl1);

  //log_info("%d: r1 = %f\n", bsp_pid(), r1);
  bsp_sync(); 
  clvl1_sync[bsp_pid()]++;
  
  // getting results from level 1
  level1_results = vecallocd(bsp_nprocs());
  level1_results[bsp_pid()] = r1;
  bsp_push_reg(level1_results, SZDBL*bsp_nprocs());
  
  
  bsp_sync(); 
  clvl1_sync[bsp_pid()]++;
  // waiting for results in level 1
 
  bsp_sync();
  clvl1_sync[bsp_pid()]++;

  result = 0;
  if (bsp_pid() == 0) { // master
    for (i=0; i< bsp_nprocs(); i++) {
      result += level1_results[i];
    }
  }

  bsp_end();
  // end: Multibsp level 1


  double* level2_results = vecallocd(bsp_nprocs());
  bsp_push_reg(level2_results, SZDBL*bsp_nprocs());
  bsp_sync();
  clvl2_sync[bsp_pid()]++;
  bsp_put(0, &result, level2_results, bsp_pid()*SZDBL, SZDBL);
  bsp_sync();
  clvl2_sync[bsp_pid()]++;
  clvl2_comm[bsp_pid()]+= SZDBL;
  total = 0;
  if(bsp_pid() == 0 ) {
    for(i=0; i<bsp_nprocs(); i++) {
      total += level2_results[i];
    }
  }
  bsp_end();
  // end: Multibsp level 2

  return total;

} /* end bspinprod */


int main(int argc, char **argv)
{

  double* v;// es cuadratico por eso uso solo 1
  double r;
  double time;
  CLOCK_TYPE chrono; /* Total time of execution */

  int amount;

  // start: define tree with multibsp-discover 
  tree_root = NULL;
  //multibsp_tree_node_t tnode     = NULL;
  int i;
  int j;

  
  tree_root = multibsp_discover_new();
  //printf("N2: %d\n", N2);
  //printf("N1: %d\n", N1);

  // end: define tree with multibsp-discover 
  /*
  switch(argc) {
    case 1: amount = 1024; break;
    case 2: amount = atoi(argv[1]); break;
    default: 
            printf("Usage: ./prodin <<size of vector>>. If you don't put size, the default is 1024\n"); break;
  };
  */
  for (amount = 1024; amount < ((1024*1000000)+1); amount*=2) {
    clean_counters();
    //printf("x%d\n", amount);
    v = vecallocd(amount);
    for (i=0; i<amount; i++) v[i] = 1;
    CLOCK_Start(chrono);
    r = multibsp_level2_inprod(v, amount);

    CLOCK_End(chrono, time);
    //printf("\nResult = %f\n", r);
    printf("%d %.6lf\n", i, time/10);
    print_counters(amount);
    free(v);

  }
  
  exit(0);

} /* end main */


