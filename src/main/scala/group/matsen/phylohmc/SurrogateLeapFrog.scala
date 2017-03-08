package group.matsen.phylohmc

import spire.algebra.{Field, Signed}
import spire.std.seq._
import spire.syntax.order._
import spire.syntax.vectorSpace._

trait SurrogateLeapFrog[R, N, D <: Int with Singleton] extends NumericalDynamics[R, N, D] {

  val barrier = Field[R].fromDouble(0.0)

  val delta: R

  def g(t: R): R = if (t >= delta)
    t
  else
    1 / (2 * delta) * (t * t + delta * delta)

  def `g'`(t: R): R = if (t >= delta)
    Field[R].one
  else
    t / delta

  override def U(q: Tree[R, N]) = (
    super.U(q)._1,
    (super.U(q.modifyLengths(_.map(g)))._2, q.lengths.view.map(`g'`)).zipped.map(Field[R].times)
  )

  def surrogateU(q: Tree[R, N]) = super.U(q.modifyLengths(_.map(g)))._1

  def leapprog(eps: R)(z: ZZ): ZZ = {
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
