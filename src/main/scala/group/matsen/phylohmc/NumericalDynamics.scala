package group.matsen.phylohmc

trait NumericalDynamics[R, N, D <: Int with Singleton] extends PhyloHMC[R, N, D] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N]

  def simulateDynamics(z: Z[R, N]): Z[R, N] = (0 until L).foldLeft(z)((z, _) => leapprog(eps)(z))

}
