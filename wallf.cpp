#include "globals.h"


void find_new_surface_bonds( void ) ;
void surface_bond_forces( void ) ;




void wallf(){

int i, j =Dim -1 ;

#pragma omp parallel for private(i)
   for( i = 0; i<nstot ; i++){
	if ( x[i][j] < wall_thick )
	   f[i][j] += wall_lamb[tp[i]] ; 
	if ( x[i][j] > (L[j] - wall_thick) )
	   f[i][j] -= wall_lamb[tp[i]]  ;
   }

  
  if ( graft_surface && step > graft_surf_wait ) {
    surface_bond_forces() ;
    find_new_surface_bonds() ;
  }

}

void write_surface_bonds() {
  int i , j ;
  FILE *otp ;
  otp = fopen( "surface_bonds.out", "w" ) ;

  for ( i=0 ; i<n_surf_bonds ; i++ ) {
    fprintf( otp , "%d %d  " , i , surf_bonds[i] ) ;
    fprintf( otp , "%d  " , tp[ surf_bonds[i] ] ) ;
    for ( j=0 ; j<Dim ; j++ ) 
      fprintf( otp , "%lf " , surf_bond_req[i][j] ) ;
    fprintf(otp, "\n" ) ;
  }

  fclose( otp ) ;

}

void surface_bond_forces() {

  int i , ind, j ;
  double mdr2, dr[Dim] ;
  
  surf_bond_energy = 0.0 ;

#pragma omp parallel for \
  private(i, ind, dr, mdr2) \
  reduction(+:surf_bond_energy)
  for ( i=0 ; i<n_surf_bonds ; i++ ) {
    ind = surf_bonds[i] ;

    mdr2 = pbc_mdr2( x[ind], surf_bond_req[i], dr ) ;

    surf_bond_energy += mdr2 * 0.5 * surf_bond_k ;

    for ( j=0 ; j<Dim ; j++ ) 
      f[ind][j] -= surf_bond_k * dr[j] ;
  }


}

void find_new_surface_bonds( ) {

  int i, j, ind ;
  ind = 0 ;

  // Check diblock ends for new bonds //
  for ( i=0 ; i<nD ; i++ ) {
    if ( surf_bond_flags[ind] == 0 && x[ind][Dim-1] < wall_thick ) {
      surf_bonds[ n_surf_bonds ] = ind ;
      surf_bond_flags[ind] = 1 ;
      xc[ind] = "C" ;
      for ( j=0 ; j<Dim ; j++ ){ 
       if(j != Dim-1) surf_bond_req[ n_surf_bonds ][j] = x[ind][j] ;
       else   surf_bond_req[ n_surf_bonds ][j] = wall_thick;
      }

      n_surf_bonds++ ;
    }

    ind += Nda + Ndb ;
  }


  // Check homopolymer ends for new bonds //
  for ( i=0 ; i<nA ; i++ ) {
    if ( surf_bond_flags[ind] == 0 && x[ind][Dim-1] < wall_thick ) {
      surf_bonds[ n_surf_bonds ] = ind ;
      surf_bond_flags[ind] = 1 ;
      xc[ind] = "C" ;
      for ( j=0 ; j<Dim ; j++ ) 
        surf_bond_req[ n_surf_bonds ][j] = x[ind][j] ;

      n_surf_bonds++ ;
    }

    ind += Nha ;
  }

  // Check homopolymer ends for new bonds //
  for ( i=0 ; i<nB ; i++ ) {
    if ( surf_bond_flags[ind] == 0 && x[ind][Dim-1] < wall_thick ) {
      surf_bonds[ n_surf_bonds ] = ind ;
      surf_bond_flags[ind] = 1 ;
      xc[ind] = "C" ;
      for ( j=0 ; j<Dim ; j++ ) 
        surf_bond_req[ n_surf_bonds ][j] = x[ind][j] ;

      surf_bond_flags[ ind ] = 1;

      n_surf_bonds++ ;
    }

    ind += Nhb ;
  }
}


void init_surface_bonds() {
  int i ;
  int max_bonds = 2*(nA+nB+nD) ;

  n_surf_bonds = 0 ;

  surf_bonds = ( int* ) calloc( max_bonds,  sizeof( int ) ) ;
  surf_bond_flags = ( int* ) calloc( nstot , sizeof( int ) ) ;
  surf_bond_req = ( double** ) calloc( max_bonds , sizeof( double* ) ) ;

  for ( i=0 ; i<max_bonds ; i++ ) 
    surf_bond_req[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;

  for ( i=0 ; i<nstot ; i++ ) 
    surf_bond_flags[i] = 0 ;

}

