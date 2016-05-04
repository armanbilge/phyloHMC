package org.fhcrc.matsen.phylohmc

import spire.algebra.{Field, Signed}
import spire.algebra.Order._
import spire.std.seq._
import spire.syntax.vectorSpace._
import spire.syntax.order._

trait ReflectiveLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val Kp = K(pp)
    val qp = z.q.modifyLengths(l => (l + eps *: Kp._2).map(Signed[R].abs))
    val zp = (z.q.lengths, Kp._2).zipped.map(- _ / _).zipWithIndex.filter(Function.tupled((e, _) => Field[R].zero <= e && e <= eps)).sortBy(_._1).view.map(_._2).foldLeft(z.copy(q = qp, p = pp)(U(qp), Kp)) { (z, i) =>
      val q = (if (z.q.isInternal(i)) rng.nextInt(1) else 0) match {
        case 0 => z.q
        case 1 => z.q.nni(i, false)
        case 2 => z.q.nni(i, true)
      }
      val p = z.p.updated(i, -z.p(i))
      z.copy(q, p)(U(q), K(p))
    }
    val ppp = zp.p - halfEps *: zp.dU
    zp.copy(p = ppp)(_K = K(ppp))
  }

}
