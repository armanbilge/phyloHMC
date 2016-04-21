package org.fhcrc.matsen.phylohmc

import spire.algebra.Field
import spire.syntax.field._
import spire.syntax.order._

trait VuPhyloHMC[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val (zp, e) = z.q.branches.map(b => (b, solveForEps(z)(b))).filter(x => x._2 <= eps).toList.sortWith(_._2 < _._2).foldLeft((z, Field[R].zero)) { (ze, bt) =>
      val (z, e) = ze
      val (b, t) = bt
      val zp = leapfrog(t - e)(z)
      (rng.nextInt(3) match {
        case 0 => zp
        case 1 => zp.nni(b, false, U, K)
        case 2 => zp.nni(b, true, U, K)
      }, t)
    }
    leapfrog(eps - e)(zp)
  }

}
