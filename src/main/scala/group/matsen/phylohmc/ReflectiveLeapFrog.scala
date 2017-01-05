package group.matsen.phylohmc

import spire.algebra.{Field, Signed}
import spire.std.seq._
import spire.syntax.order._
import spire.syntax.vectorSpace._

trait ReflectiveLeapFrog[R, N, D <: Int with Singleton] extends NumericalDynamics[R, N, D] {

  val barrier = Field[R].fromDouble(0.0)

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val (_, dK) = K(pp)
    val qp = z.q.modifyLengths(_ + eps *: dK)
    val qpp = qp.modifyLengths(_.zipWithIndex.map(Function.tupled((li, i) => if (qp.isInternal(i)) Signed[R].abs(li) else Signed[R].abs(li - barrier) + barrier)))
    val Up = U(qpp)
    val ppp = (qp.lengths, pp, pp.indices).zipped.map((qi, pi, i) => if (qi < (if (qp.isInternal(i)) Field[R].zero else barrier)) -pi else pi) - halfEps *: Up._2
    Z(qpp, ppp)(Up, K(ppp))
  }

}
