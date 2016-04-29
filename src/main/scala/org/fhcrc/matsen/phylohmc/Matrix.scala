package org.fhcrc.matsen.phylohmc

import spire.algebra.Field
import spire.std.seq._
import spire.syntax.innerProductSpace._

class Matrix[@specialized(Double) R : Field](val size: Int, val values: Vector[R]) {

  require(values.length == size * size)

  @inline def index(i: Int, j: Int): Int = size * i + j

  def apply(i: Int, j: Int): R = values(index(i, j))

  def updated(i: Int, j: Int, r: R): Matrix[R] = new Matrix(size, values.updated(index(i, j), r))

  lazy val indices: Traversable[(Int, Int)] = (for (i <- 0 until size; j <- 0 until size) yield (i, j)).toVector

  lazy val rows: IndexedSeq[IndexedSeq[R]] = values.sliding(size, size).toVector

  lazy val columns: IndexedSeq[IndexedSeq[R]] = (for (j <- 0 until size) yield (for (i <- 0 until size) yield apply(i, j)).toVector).toVector

  def *(x: IndexedSeq[R]): IndexedSeq[R] = rows.map(_ dot x)

}

object Matrix {

  def apply[@specialized(Double) R : Field](values: R*): Matrix[R] = new Matrix(Math.sqrt(values.size).toInt, values.toVector)

  def apply[@specialized(Double) R : Field](n: Int)(f: (Int, Int) => R): Matrix[R] = new Matrix(n, Vector.tabulate(n * n)(k => f(k / n, k % n)))

}
