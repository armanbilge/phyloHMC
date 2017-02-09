package group.matsen.phylohmc

trait RandomizedLength[R, N, D <: Int with Singleton] extends NumericalDynamics[R, N, D] {

  override def simulateDynamics(z: Z[R, N, G]): Z[R, N, G] = (0 until rng.nextInt(1, L)).foldLeft(z)((z, _) => leapprog(eps)(z))

}
