package org.fhcrc.matsen.phylohmc

import spire.algebra.Signed
import spire.syntax.field._
import spire.syntax.order._

trait ReflectiveLeapFrog[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val zp = leapfrog(eps)(z)
    val qp = zp.q.mapLengths(Signed[R].abs(_))
    val pp = zp.p.mapLengths((b, l) => if (zp.q(b) < 0) -l else l)
    zp.copy(q = qp, p = pp)(U(qp), K(zp.Minv)(pp))
  }

}
