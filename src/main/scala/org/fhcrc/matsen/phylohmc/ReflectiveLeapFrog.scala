package org.fhcrc.matsen.phylohmc

import spire.algebra.Signed
import spire.std.seq._
import spire.syntax.vectorSpace._
import spire.syntax.order._

trait ReflectiveLeapFrog[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val (_, dK) = K(pp)
    val qp = z.q.modifyLengths(_ + eps *: dK)
    val qpp = qp.modifyLengths(_.map(Signed[R].abs))
    val Up = U(qpp)
    val ppp = (qp.lengths, pp).zipped.map((qi, pi) => if (qi < 0) -pi else pi) - halfEps *: Up._2
    Z(qpp, ppp)(Up, K(ppp))
  }

}
