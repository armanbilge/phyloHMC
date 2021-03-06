package group.matsen.phylohmc

import spire.algebra.Order._
import spire.algebra.{Field, Signed}
import spire.std.seq._
import spire.syntax.order._
import spire.syntax.vectorSpace._

trait ReflectiveLeapProg[R, N, D <: Int with Singleton] extends NumericalDynamics[R, N, D] {

  val barrier = Field[R].fromDouble(0.0)

  def leapprog(eps: R)(z: ZZ): ZZ = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val Kp = K(pp)
    val qp = z.q.modifyLengths(l => l + eps *: Kp._2)
    val qpp = qp.modifyLengths(l => l.zipWithIndex.map(Function.tupled((li, i) => if (qp.isInternal(i)) Signed[R].abs(li) else Signed[R].abs(li - barrier) + barrier)))
    val ppp = (qp.lengths, pp, pp.indices).zipped.map((qi, pi, i) => if (qp.isInternal(i)) pi else if (qi < barrier) -pi else pi)
    val zp = (z.q.lengths, Kp._2).zipped.map(- _ / _).zipWithIndex.filter(Function.tupled((_, i) => z.q.isInternal(i))).filter(Function.tupled((e, _) => Field[R].zero <= e && e <= eps)).sortBy(_._1).view.map(_._2).foldLeft(z.copy(q = qpp, p = ppp)(U(qpp), K(ppp))) { (z, i) =>
      val q = (if (z.q.isInternal(i)) rng.nextInt(3) else 0) match {
        case 0 => z.q
        case 1 => z.q.nni(i, false)
        case 2 => z.q.nni(i, true)
      }
      val p = z.p.updated(i, -z.p(i))
      z.copy(q, p)(U(q), K(p))
    }
    val pppp = zp.p - halfEps *: zp.dU
    zp.copy(p = pppp)(_K = K(pppp))
  }

}
