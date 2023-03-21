function two_foci_demo ( )

%*****************************************************************************80
%
%% TWO_FOCI_DEMO demonstrates MESH2D on a box with two foci.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    13 May 2017
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TWO_FOCI_DEMO:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Demonstrate MESH2D on a the two foci problem.\n' );

  warning off
%
%  1) Do the region.
%
  figure ( )
  v = [ 0.0, 0.0; 1.0, 0.0; 1.0, 1.0; 0.0, 1.0 ];
  [ p, t ] = mesh2d ( v );
  print ( '-dpng', 'two_foci_1.png' );
%
%  2) Do the region with a maximum H.
%
  figure ( )
  v = [ 0.0, 0.0; 1.0, 0.0; 1.0, 1.0; 0.0, 1.0 ];
  hdata = [];
  hdata.hmax = 0.1;
  [ p, t ] = mesh2d ( v, [], hdata );
  print ( '-dpng', 'two_foci_2.png' );
%
%  3) Do the region, with H specifying two special foci.
%
  figure ( )
  v = [ 0.0, 0.0; 1.0, 0.0; 1.0, 1.0; 0.0, 1.0 ];
  hdata = [];
  hdata.fun = @h_function;
  [ p, t ] = mesh2d ( v, [], hdata );
  print ( '-dpng', 'two_foci_3.png' );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TWO_FOCI_INTERIOR_DEMO:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  return
end
function h = h_function ( x, y )

%*****************************************************************************80
%
%% H_FUNCTION is a size function for the two foci problem.
%
%  Discussion:
%
%    The mesh size is 0.01 near (0.25,1.0) and near (0.75,1.0).
%
  h1 = 0.01 + 0.1 * sqrt ( ( x - 0.25 ).^2  + ( y - 1.0 ).^2 );
  h2 = 0.01 + 0.1 * sqrt ( ( x - 0.75 ).^2  + ( y - 1.0 ).^2 );
  h = min ( h1, h2 );

  return
end