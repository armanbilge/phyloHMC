package group.matsen.phylohmc

import shapeless.Witness
import spire.algebra.{Field, Ring, VectorSpace}
import spire.std.seq._
import spire.syntax.innerProductSpace._

class Matrix[N <: Int with Singleton : Witness.Aux, R : Field](val values: Vector[R]) {

  val size = implicitly[Witness.Aux[N]].value
  require(values.length == size * size)

  @inline final def index(i: Int, j: Int): Int = size * i + j

  def apply(i: Int, j: Int): R = values(index(i, j))

  def updated(i: Int, j: Int, r: R): Matrix[N, R] = new Matrix(values.updated(index(i, j), r))

  lazy val indices: Traversable[(Int, Int)] = (for (i <- 0 until size; j <- 0 until size) yield (i, j)).toVector

  lazy val rows: IndexedSeq[IndexedSeq[R]] = values.sliding(size, size).toVector

  lazy val columns: IndexedSeq[IndexedSeq[R]] = (for (j <- 0 until size) yield (for (i <- 0 until size) yield apply(i, j)).toVector).toVector

  def *(x: IndexedSeq[R]): IndexedSeq[R] = rows.map(_ dot x)

}

object Matrix {

  def apply[N <: Int with Singleton : Witness.Aux, R : Field](values: R*): Matrix[N, R] = new Matrix(values.toVector)

  def apply[N <: Int with Singleton : Witness.Aux, R : Field](f: (Int, Int) => R): Matrix[N, R] = {
    val n = implicitly[Witness.Aux[N]].value
    new Matrix(Vector.tabulate(n * n)(k => f(k / n, k % n)))
  }

  implicit def matrixIsRing[N <: Int with Singleton : Witness.Aux, R : Field]: Ring[Matrix[N, R]] = new Ring[Matrix[N, R]] {

    override def negate(x: Matrix[N, R]): Matrix[N, R] = Matrix((i, j) => -x(i, j))

    override def zero: Matrix[N, R] = Matrix((i, j) => Field[R].zero)

    override def plus(x: Matrix[N, R], y: Matrix[N, R]): Matrix[N, R] = Matrix((i, j) => x(i, j) + y(i, j))

    override def minus(x: Matrix[N, R], y: Matrix[N, R]): Matrix[N, R] = Matrix((i, j) => x(i, j) - y(i, j))

    override def one: Matrix[N, R] = Matrix((i, j) => if (i == j) Field[R].one else Field[R].zero)

    override def times(x: Matrix[N, R], y: Matrix[N, R]): Matrix[N, R] = Matrix((i, j) => x.rows(i) dot y.columns(j))

  }

  implicit def matrixIsVectorSpace[N <: Int with Singleton : Witness.Aux, R](implicit f: Field[R]): VectorSpace[Matrix[N, R], R] = new VectorSpace[Matrix[N, R], R] {

    override def scalar: Field[R] = f

    override def timesl(r: R, v: Matrix[N, R]): Matrix[N, R] = Matrix[N, R]((i: Int, j: Int) => r * v(i, j))

    override def divr(v: Matrix[N, R], r: R): Matrix[N, R] = Matrix[N, R]((i: Int, j: Int) => v(i, j) / r)

    override def negate(x: Matrix[N, R]): Matrix[N, R] = Ring[Matrix[N, R]].negate(x)

    override def zero: Matrix[N, R] = Ring[Matrix[N, R]].zero

    override def plus(x: Matrix[N, R], y: Matrix[N, R]): Matrix[N, R] = Ring[Matrix[N, R]].plus(x, y)

    override def minus(x: Matrix[N, R], y: Matrix[N, R]): Matrix[N, R] = Ring[Matrix[N, R]].minus(x, y)

  }

}
