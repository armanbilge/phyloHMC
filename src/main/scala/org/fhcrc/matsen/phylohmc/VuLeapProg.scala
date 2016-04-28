package org.fhcrc.matsen.phylohmc

import spire.algebra.Field
import spire.syntax.field._
import spire.syntax.order._

trait VuLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val (zp, e) = z.q.lengths.indices.view.map(solveForEps(z)).filter(_ <= eps).zipWithIndex.sortWith(_._1 < _._1).foldLeft((z, Field[R].zero)) { (ze, ti) =>
      val (z, e) = ze
      val (t, i) = ti
      val zp = leapfrog(t - e)(z)
      val q = (if (zp.q.isInternal(i)) rng.nextInt(3) else 0) match {
        case 0 => zp.q
        case 1 => zp.q.nni(i, false)
        case 2 => zp.q.nni(i, true)
      }
      val p = z.p.updated(i, -z.p(i))
      (zp.copy(q = q, p = p)(U(q), K(p)), t)
    }
    leapfrog(eps - e)(zp)
  }

}
