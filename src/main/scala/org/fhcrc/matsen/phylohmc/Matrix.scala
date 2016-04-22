package org.fhcrc.matsen.phylohmc

class Matrix[R](n: Int, values: Vector[R]) {

  require(values.length == n * n)

  @inline def index(i: Int, j: Int): Int = n * i + j

  def apply(i: Int, j: Int): R = values(index(i, j))

  def updated(i: Int, j: Int, r: R): Matrix[R] = new Matrix(n, values.updated(index(i, j), r))

}

object Matrix {

  def apply[R](values: R*): Matrix[R] = new Matrix(Math.sqrt(values.size).toInt, values.toVector)

  def apply[R](n: Int)(f: (Int, Int) => R): Matrix[R] = new Matrix(n, Vector.tabulate(n * n)(k => f(k / n, k % n)))

  def apply[R](n: Int)(f: Int => Int => R): Matrix[R] = new Matrix(n, Vector.tabulate(n * n)(k => f(k / n)(k % n)))

}
