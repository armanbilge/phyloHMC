package org.fhcrc.matsen.phylohmc

import spire.algebra.{Field, Order}
import spire.syntax.field._
import spire.syntax.order._

import scala.annotation.tailrec

trait VuLeapProg[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {

    @tailrec
    def recurse(z: Z[R, N], eps: R): Z[R, N] = {
      z.q.lengths.indices.view.map(solveForEps(z)).zipWithIndex.filter(_._1.isDefined).map(ei => (ei._1.get, ei._2)).filter(_._1 <= eps).reduceOption(Order.by[(R, Int), R](_._1).min) match {
        case Some((e, i)) =>
          val zp = leapfrog(e)(z)
          val q = ((if (zp.q.isInternal(i)) rng.nextInt(3) else 0) match {
            case 0 => zp.q
            case 1 => zp.q.nni(i, false)
            case 2 => zp.q.nni(i, true)
          }).modifyLengths(_.updated(i, Field[R].zero))
          val p = zp.p.updated(i, -zp.p(i))
          recurse(zp.copy(q = q, p = p)(U(q), K(p)), eps - e)
        case None => leapfrog(eps)(z)
      }
    }

    recurse(z, eps)

  }

}
