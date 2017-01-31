package group.matsen.phylohmc

import spire.algebra.Order._
import spire.algebra.{Field, NRoot}
import spire.std.seq._
import spire.syntax.order._
import spire.syntax.vectorSpace._

trait SurrogateLeapProg[R, N, D <: Int with Singleton] extends NumericalDynamics[R, N, D] {

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
    val Kp = K(pp)
    val zp = (z.q.lengths, Kp._2).zipped.map(- _ / _).zipWithIndex.sortBy(_._1).view.foldLeft(z.copy(p = pp)(_K = Kp)) { (z, ei) =>
      val (e, i) = ei
      val pi = z.p(i)
      if (Field[R].zero <= e && e <= eps) {
        val q = z.q.modifyLengths(_.updated(i, Field[R].zero))
        val qp = (if (q.isInternal(i)) rng.nextInt(3) else 0) match {
          case 0 => q
          case 1 => q.nni(i, false)
          case 2 => q.nni(i, true)
        }
        val Up = U(q)
        val deltaU = surrogateU(qp) - surrogateU(q)
        if (pi * pi > 2 * deltaU) {
          val p = z.p.updated(i, NRoot[R].sqrt(pi * pi - 2 * deltaU))
          val qpp = qp.modifyLengths(_.updated(i, (eps - e) * p(i)))
          z.copy(qpp, p)(U(qpp), K(p))
        } else {
          val p = z.p.updated(i, -z.p(i))
          val qpp = q.modifyLengths(_.updated(i, (eps - e) * p(i)))
          z.copy(qpp, p)(U(qpp), K(p))
        }
      } else {
        val q = z.q.modifyLengths(_.updated(i, z.q.lengths(i) + eps * z.dK(i)))
        z.copy(q = q)(U(q))
      }
    }
    val ppp = zp.p - halfEps *: zp.dU
    zp.copy(p = ppp)(_K = K(ppp))
  }

}
