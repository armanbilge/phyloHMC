package org.fhcrc.matsen.phylohmc

import spire.algebra.Signed
import spire.syntax.field._
import spire.syntax.order._

trait ReflectiveLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val zp = leapfrog(eps)(z)
    val zpp = zp.copy(q = zp.q.mapLengths(Signed[R].abs(_)))()
    zpp.q.branches.map(b => (b, solveForEps(z)(b))).filter(x => x._2 <= eps).toList.sortWith(_._2 < _._2).view.map(_._1).foldLeft(zpp) { (z, b) =>
      val zp = z.copy(p = z.p.updated(b, z.p(b)))(_K = (z.k, -z.dK))
      (if (zp.q.isInternal(b)) rng.nextInt(3) else 0) match {
        case 0 => zp
        case 1 => zp.nni(b, false, U, K)
        case 2 => zp.nni(b, true, U, K)
      }
    }
  }

}
