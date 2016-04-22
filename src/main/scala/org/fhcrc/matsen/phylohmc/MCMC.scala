package org.fhcrc.matsen.phylohmc

object MCMC {

  def apply[T](start: T, length: Int)(f: T => T): TraversableOnce[T] = (0 to length).view.scanLeft(start)((t, _) => f(t))

}
