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