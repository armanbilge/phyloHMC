package org.fhcrc.matsen.phylohmc

import spire.algebra.Signed
import spire.syntax.field._
import spire.syntax.order._

trait ReflectiveLeapFrog[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val zp = leapfrog(eps)(z)
    val qp = zp.q.modifyLengths(_.map(Signed[R].abs(_)))
    val pp = (zp.q.lengths, zp.p).zipped.map((qi, pi) => if (qi < 0) -pi else pi)
    zp.copy(q = qp, p = pp)(U(qp), K(pp))
  }

}
