#include "globals.h"

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "dyft.input" , "r" ) ;

  fscanf(inp , "%d", &Nsp);
  fgets( tt, 80, inp);
  fscanf( inp , "%d %d" , &Nda , &Ndb ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHA , &Nha ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHB , &Nhb ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%lf %d" , &phiHC , &Nhc ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%lf" , &CG_ratio ) ;
  fgets( tt , 80 , inp ) ;

  
  // Blank line //
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &phiP ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d %d" , &sigma, &Ng, &uni_sig ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Rp ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Xi ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &A_partics ) ;
  fgets( tt , 80 , inp ) ;
  
  
  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf" , &C ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %lf %lf" , &chiAB, &chiAC, &chiBC ) ;
  fgets( tt , 80 , inp ) ;
  cout<<"chiAB "<<chiAB<<" chiAC "<<chiAC<<" chiBC "<<chiBC<<endl;
  fscanf( inp , "%lf %lf" , &kappa ,&kappa_p) ;
  fgets( tt , 80 , inp ) ;
   cout<<"kappa "<<kappa<<" "<<kappa_p<<endl;
  fscanf( inp , " %lf %lf %lf %lf %lf" , &Diff[0] , &Diff[1], &Diff[3] ,&Diff[2],&Diff_rot ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp , "%lf" , &L[i] ) ; 
  fgets( tt , 80 , inp ) ;

  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp , "%d" , &Nx[i] ) ; 
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf", &delt ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &pmeorder ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%d" , &nsteps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &print_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &pre_equil_steps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d %d" , &sample_wait , &sample_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &stress_freq ) ;
  fgets( tt , 80 , inp ) ;
  printf("stress_freq: %d\n" , stress_freq ) ;
  fscanf( inp , "%d" , &frame_freq ) ;
  fgets( tt , 80 , inp ) ;

  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d %d" , &optm_L,&L_dim ) ;
//  if(optm_L>0){ 
    printf("doing L optm on dim %d\n", L_dim) ;
    fgets( tt , 80 , inp ) ;
    fscanf( inp , "%d " , &L_fren ) ;
    fgets( tt , 80 , inp ) ;
    fscanf( inp , "%d " , &L_aver ) ;
    printf("Lfren %d and Laver %d ",L_fren, L_aver);
    fgets( tt , 80 , inp ) ;
    fscanf( inp , "%lf " , &lam_L ) ;
    printf("lam_L %f\n",lam_L);
 

	
/*  }
  else{
    fgets( tt , 80 , inp ) 	;	
    fgets( tt , 80 , inp )      ;
    fgets( tt , 80 , inp )      ;
  }*/
  fgets( tt , 80 , inp )      ;
  fgets( tt , 80 , inp )  ;

 
  fscanf( inp , "%d %lf %lf %lf %lf " , &wall_para,&wall_lamb[0],&wall_lamb[1],&wall_lamb[2],&wall_lamb[3]) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%lf " , &wall_thick ) ;
  fgets( tt , 80 , inp ) ;

  fgets( tt , 80 , inp ) ;

  fscanf( inp ,"%d %d ", &graft_surface , &graft_surf_wait); 
  fgets( tt , 80 , inp ) ;
   fscanf( inp ,"%lf ",&surf_bond_k); 
   fgets( tt , 80 , inp ) ;
  fclose( inp ) ;
}
