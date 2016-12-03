package group.matsen.phylohmc

object MCMC {

  def apply[T](start: T)(f: T => T): Stream[T] = Stream.from(0).scanLeft(start)((t, _) => f(t))

}
