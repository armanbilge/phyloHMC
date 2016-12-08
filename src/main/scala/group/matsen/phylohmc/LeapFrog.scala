package group.matsen.phylohmc

import spire.std.seq._
import spire.syntax.vectorSpace._

trait LeapFrog[R, N, D <: Int with Singleton] extends PhyloHMC[R, N, D] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val (_, dK) = K(pp)
    val qp = z.q.modifyLengths(_ + eps *: dK)
    val Up = U(qp)
    val ppp = pp - halfEps *: Up._2
    Z(qp, ppp)(Up, K(ppp))
  }

}
