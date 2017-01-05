package group.matsen.phylohmc

import spire.algebra.{Field, NRoot, Order, Signed}
import spire.std.seq._
import spire.syntax.innerProductSpace._
import spire.syntax.order._

import scala.annotation.tailrec

trait VuLeapProg[R, N, D <: Int with Singleton] extends NumericalDynamics[R, N, D] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {

    @tailrec
    def recurse(z: Z[R, N], eps: R): Z[R, N] = {

      def solveForEps(z: Z[R, N])(i: Int): Option[R] = {
        val twoa = - invM.rows(i) dot z.dU
        val b = invM.rows(i) dot z.p
        val c = z.q.lengths(i)
        val mbdiv2a = - b / twoa
        val discriminantdiv2a = NRoot[R].sqrt(b * b - 2 * twoa * c) / twoa
        val eps1 = mbdiv2a + discriminantdiv2a
        val eps2 = mbdiv2a - discriminantdiv2a
        (eps1 > 0, eps2 > 0) match {
          case (true, true) => Some(eps1 min eps2)
          case (false, false) => None
          case _ => Some(eps1 max eps2)
        }
      }

      def leapfrog(eps: R, i: Option[Int] = None)(z: Z[R, N]): Z[R, N] = {
        val halfEps = eps / 2
        val pp = z.p - halfEps *: z.dU
        val (_, dK) = K(pp)
        val qp = z.q.modifyLengths(_ + eps *: dK).modifyLengths(l => i.fold(l)(l.updated(_, Field[R].zero)))
        val Up = U(qp)
        val ppp = pp - halfEps *: Up._2
        Z(qp, ppp)(Up, K(ppp))
      }

      z.q.lengths.indices.view.map(solveForEps(z)).zipWithIndex.filter(_._1.isDefined).map(ei => (ei._1.get, ei._2)).filter(_._1 <= eps).reduceOption(Order.by[(R, Int), R](_._1).min) match {
        case Some((e, i)) =>
          val zp = leapfrog(e, Some(i))(z)
          val q = (if (zp.q.isInternal(i)) rng.nextInt(3) else 0) match {
            case 0 => zp.q
            case 1 => zp.q.nni(i, false)
            case 2 => zp.q.nni(i, true)
          }
          // TODO What is the best way to handle this?
          val p = zp.p.updated(i, Signed[R].abs(zp.p(i)))
          recurse(zp.copy(q = q, p = p)(U(q), K(p)), eps - e)
        case None => leapfrog(eps)(z)
      }
    }

    recurse(z, eps)

  }

}
