package org.fhcrc.matsen.phylohmc

import spire.algebra.Order._
import spire.algebra.Signed
import spire.syntax.field._
import spire.syntax.order._

trait ReflectiveLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val zp = leapfrog(eps)(z)
    val qp = zp.q.modifyLengths(_.map(Signed[R].abs))
    qp.lengths.indices.view.map(solveForEps(z)).zipWithIndex.filter(_._1.isDefined).map(ei => (ei._1.get, ei._2)).filter(_._1 <= eps).sortBy(_._1).view.map(_._2).foldLeft(zp.copy(q = qp)(_U = U(qp))) { (z, i) =>
      val q = (if (z.q.isInternal(i)) rng.nextInt(3) else 0) match {
        case 0 => z.q
        case 1 => z.q.nni(i, false)
        case 2 => z.q.nni(i, true)
      }
      val p = z.p.updated(i, -z.p(i))
      z.copy(q = q, p = p)(U(q), K(p))
    }
  }

}
