package org.fhcrc.matsen.phylohmc

import spire.algebra.Order._
import spire.algebra.Signed
import spire.std.seq._
import spire.syntax.vectorSpace._
import spire.syntax.order._

trait ReflectiveLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val (_, dK) = K(pp)
    val qp = z.q.modifyLengths(l => (l + eps *: dK).map(Signed[R].abs))
    val zp = qp.lengths.indices.view.map(solveForEps(z)).zipWithIndex.filter(_._1.isDefined).map(ei => (ei._1.get, ei._2)).filter(_._1 <= eps).sortBy(_._1).view.map(_._2).foldLeft(z.copy(q = qp)(_U = U(qp))) { (z, i) =>
      val q = (if (z.q.isInternal(i)) rng.nextInt(3) else 0) match {
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
