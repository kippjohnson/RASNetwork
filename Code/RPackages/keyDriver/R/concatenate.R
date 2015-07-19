concatenate <- function( myvect , mysep="" )
{
  if ( is.null( myvect ) )
  {
    return ( "" )
  }
  else if ( isTRUE( all.equal( length( myvect ) , 1 ) ) )
  {
    return ( as.character( myvect ) )
  }
  return( paste( as.character( myvect ), sep = "" , collapse = mysep ) )
}

