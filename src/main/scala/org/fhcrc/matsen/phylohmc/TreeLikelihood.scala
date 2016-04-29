package org.fhcrc.matsen.phylohmc

import spire.algebra.{Field, Trig}
import spire.std.seq._
import spire.syntax.innerProductSpace._

import scala.collection.GenMap

class TreeLikelihood[@specialized(Double) R : Field : Trig, @specialized(Int) N](val patterns: Patterns, val model: SubstitutionModel[R], val mu: R) extends (Tree[R, N] => R) {

  override def apply(t: Tree[R, N]): R = {

    def recurse(parent: N, child: N): GenMap[Pattern, IndexedSeq[R]] = {

      @inline def internalInternal(x: IndexedSeq[R], y: IndexedSeq[R], A: Matrix[R], B: Matrix[R]): IndexedSeq[R] = (A.columns, B.columns).zipped.map((a, b) => (a dot x) * (b dot y))

      @inline def internalLeaf(x: IndexedSeq[R], j: Int, A: Matrix[R], B: Matrix[R]): IndexedSeq[R] = (A.columns, B.rows(j)).zipped.map((a, b) => (a dot x) * b)

      @inline def leafLeaf(i: Int, j: Int, A: Matrix[R], B: Matrix[R]): IndexedSeq[R] = (A.rows(i), B.rows(j)).zipped.map(Field[R].times)

      val children = t.neighbors(child) - parent
      val (left, right) = (children.head, children.tail.head)
      val leftMatrix = model(t(Branch(child, left)))
      val rightMatrix = model(t(Branch(child, right)))
      (t.isInternal(left), t.isInternal(right)) match {
        case (true, true) =>
          val leftPartials = recurse(child, left)
          val rightPartials = recurse(child, right)
          patterns.par.map(p => (p, leftPartials(p), rightPartials(p))).map(Function.tupled((p, x, y) => (p, internalInternal(x, y, leftMatrix, rightMatrix)))).toMap
        case (true, false) =>
          val leftPartials = recurse(child, left)
          patterns.par.map(p => (p, leftPartials(p), p(t.taxa(right)))).map(Function.tupled((p, x, j) => (p, internalLeaf(x, j, leftMatrix, rightMatrix)))).toMap
        case (false, true) =>
          val rightPartials = recurse(child, right)
          patterns.par.map(p => (p, rightPartials(p), p(t.taxa(left)))).map(Function.tupled((p, x, j) => (p, internalLeaf(x, j, rightMatrix, leftMatrix)))).toMap
        case (false, false) =>
          patterns.par.map(p => (p, p(t.taxa(left)), p(t.taxa(right)))).map(Function.tupled((p, i, j) => (p, leafLeaf(i, j, leftMatrix, rightMatrix)))).toMap
      }

    }

    val rho = t.nodes.filter(t.isLeaf).head
    val root = t.neighbors(rho).head
    val rrMatrix = model(t(Branch(rho, root)))
    Field[R].sum(
      recurse(rho, root).par.map(Function.tupled({ (pattern, partial) =>
        val i = pattern(t.taxa(rho))
        patterns.multiplicity(pattern) * Trig[R].log(model.stationaryDistribution(i) * Field[R].sum((rrMatrix.columns(i), partial).zipped.map(Field[R].times)))
      })).toIterator
    )

  }

}
