package org.fhcrc.matsen.phylohmc

import spire.algebra.Signed
import spire.syntax.order._

trait ReflectiveLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val zp = leapfrog(eps)(z)
    val zpp = zp.copy(q = zp.q.mapLengths(Signed[R].abs(_)))()
    z.q.branches.map(b => (b, solveForEps(z)(b))).filter(x => x._2 <= eps).toList.sortWith(_._2 < _._2).view.map(_._1).foldLeft(zpp) { (z, b) =>
      rng.nextInt(3) match {
        case 0 => z
        case 1 => z.nni(b, false, U, K)
        case 2 => z.nni(b, true, U, K)
      }
    }
  }

}
