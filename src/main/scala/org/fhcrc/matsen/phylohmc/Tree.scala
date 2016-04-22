package org.fhcrc.matsen.phylohmc

import monocle.function.Index._
import monocle.std.map._
import spire.algebra.{AdditiveMonoid, Field, InnerProductSpace}

case class Tree[R : AdditiveMonoid, N](nodes: Set[N], branches: Set[Branch[N]], neighbors: Map[N, Set[N]], lengths: Map[Branch[N], R], taxa: PartialFunction[N, Taxon]) extends PartialFunction[Branch[N], R] {

  override def isDefinedAt(b: Branch[N]): Boolean = branches contains b

  override def apply(b: Branch[N]): R = lengths(b)

  def updated(b: Branch[N], r: R): Tree[R, N] = {
    require(branches contains b)
    copy(lengths = lengths + (b -> r))
  }

  def mapLengths[S : AdditiveMonoid](f: R => S): Tree[S, N] = copy(lengths = lengths.mapValues(f))

  def mapLengths[S : AdditiveMonoid](f: (Branch[N], R) => S) = copy(lengths = lengths.map(Function.tupled((b, l) => b -> f(b, l))))

  def nni(b: Branch[N], which: Boolean): Tree[R, N] = {
    require(isInternal(b))
    val u = b.head
    val v = b.tail
    val x = children(u, v).head
    val y = if (which) children(v, u).head else children(v, u).tail.head
    val (branchesp, neighborsp) = Tree.connect(u, y)(Tree.connect(v, x)(Tree.disconnect(u, x)(Tree.disconnect(v, y)(branches, neighbors))))
    copy(branches = branchesp, neighbors = neighborsp)
  }

  def children(n: N, p: N): Set[N] = neighbors(n) - p

  def isLeaf(n: N): Boolean = neighbors(n).size == 1

  def isInternal(n: N): Boolean = !isLeaf(n)

  def isPendant(b: Branch[N]): Boolean = !isInternal(b)

  def isInternal(b: Branch[N]): Boolean = isInternal(b.head) && isInternal(b.tail)

  lazy val treeLength: R = implicitly[AdditiveMonoid[R]].sum(lengths.values)

  lazy val newick: String = {

    def recurse(p: N, n: N): StringBuilder = (if (isInternal(n))
      new StringBuilder ++= "(" ++= children(n, p).map(recurse(n, _)).mkString(",") ++= ")"
    else
      new StringBuilder(taxa(n).name)
    ) ++= ":" ++= this(Branch(p, n)).toString

    val rho = nodes.view.filter(isLeaf).head
    val root = neighbors(rho).head
    val newick = recurse(rho, root)
    val i = newick.lastIndexOf(")")
    (newick.replace(i, i + 1, "," + taxa(rho)) ++= ");").toString()

  }

}

object Tree {

  def generateNeighbors[N](nodes: Set[N], branches: Set[Branch[N]]) = branches.foldLeft(Map[N, Set[N]]().withDefaultValue(Set())) { (neighbors, b) =>
    Tree.addNeighbor(b.head, b.tail)(Tree.addNeighbor(b.tail, b.head)(neighbors))
  }

  def connect[N](x: N, y: N)(bn: (Set[Branch[N]], Map[N, Set[N]])) =
    (bn._1 + Branch(x, y), addNeighbor(x, y)(addNeighbor(y, x)(bn._2)))

  def disconnect[N](x: N, y: N)(bn: (Set[Branch[N]], Map[N, Set[N]])) =
    (bn._1 - Branch(x, y), removeNeighbor(x, y)(removeNeighbor(y, x)(bn._2)))

  def addNeighbor[N](i: N, e: N) = index[Map[N, Set[N]], N, Set[N]](i).modify(_ + e)

  def removeNeighbor[N](i: N, e: N) = index[Map[N, Set[N]], N, Set[N]](i).modify(_ - e)

  def apply[R : AdditiveMonoid, N](nodes: Set[N], lengths: Map[Branch[N], R], taxa: PartialFunction[N, Taxon]): Tree[R, N] =
    Tree(nodes, lengths.keySet, generateNeighbors(nodes, lengths.keySet), lengths, taxa)

  def apply[R : AdditiveMonoid](taxa: TraversableOnce[Taxon], r: => R): Tree[R, Int] = {
    val leaves = taxa.toIndexedSeq.zipWithIndex.map(Function.tupled((t, i) => i -> t)).toMap
    val rho = leaves.keys.head
    val nodes = leaves.keys.tail.toSet
    val (root, branches) = (leaves.size until (2 * leaves.size - 2)).foldRight((nodes, Set[Branch[Int]]())) { (i, lb) =>
      val (nodes, branches) = lb
      (nodes.tail.tail + i, branches + Branch(nodes.head, i) + Branch(nodes.tail.head, i))
    }
    val branchesp = branches + Branch(rho, root.head)
    Tree((0 until (2 * leaves.size - 2)).toSet, branchesp.map(_ -> r).toMap, leaves)
  }

  implicit def TreeIsInnerProductSpace[R, N](implicit f: Field[R]) = new InnerProductSpace[Tree[R, N], R] {

    override def scalar: Field[R] = f

    override def timesl(r: R, v: Tree[R, N]): Tree[R, N] = v.copy(lengths = v.lengths.map(Function.tupled((b, l) => b -> f.times(l, r))))

    override def negate(x: Tree[R, N]): Tree[R, N] = x.copy(lengths = x.lengths.map(Function.tupled((b, l) => b -> f.negate(l))))

    override def zero: Tree[R, N] = Tree[R, N](Set[N](), Set[Branch[N]](), Map[N, Set[N]](), Map[Branch[N], R](), Map[N, Taxon]())

    override def plus(x: Tree[R, N], y: Tree[R, N]): Tree[R, N] = {
      import spire.std.map._
      import spire.syntax.vectorSpace._
      Tree[R, N](x.nodes ++ y.nodes, x.branches ++ y.branches, x.neighbors ++ y.neighbors, x.lengths + y.lengths, x.taxa.orElse(y.taxa))
    }

    override def minus(x: Tree[R, N], y: Tree[R, N]): Tree[R, N] = {
      import spire.std.map._
      import spire.syntax.vectorSpace._
      Tree[R, N](x.nodes ++ y.nodes, x.branches ++ y.branches, x.neighbors ++ y.neighbors, x.lengths - y.lengths, x.taxa.orElse(y.taxa))
    }

    override def dot(v: Tree[R, N], w: Tree[R, N]): R = {
      import spire.syntax.field._
      Field[R].sum((v.branches intersect w.branches).map(b => v(b) * w(b)))
    }

  }

}

case class Branch[N](head: N, tail: N) {

  def incident(n: N): Boolean = n match {
    case `head` => true
    case `tail` => true
    case _ => false
  }

  def getAdjacent(n: N): N = n match {
    case `head` => tail
    case `tail` => head
  }

  override val toString: String = s"($head, $tail)"

  override def equals(that: Any): Boolean = that match {
    case that: Branch[N] => (head == that.head && tail == that.tail) || (head == that.tail && tail == that.head)
    case _ => false
  }

  override val hashCode: Int = Set(head, tail).hashCode()

}
