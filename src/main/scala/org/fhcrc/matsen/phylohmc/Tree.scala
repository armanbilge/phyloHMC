package org.fhcrc.matsen.phylohmc

import monocle.function.Index._
import monocle.std.map._
import spire.algebra.AdditiveMonoid

case class Tree[@specialized(Double) R : AdditiveMonoid, @specialized(Int) N](nodes: Set[N], branchesToIndex: Map[Branch[N], Int], neighbors: Map[N, Set[N]], lengths: IndexedSeq[R], taxa: PartialFunction[N, Taxon]) extends PartialFunction[Branch[N], R] {

  val branches = branchesToIndex.keySet

  override def isDefinedAt(b: Branch[N]): Boolean = branchesToIndex contains b

  override def apply(b: Branch[N]): R = lengths(branchesToIndex(b))

  def modifyLengths(f: IndexedSeq[R] => IndexedSeq[R]): Tree[R, N] = copy(lengths = f(lengths))

  def indexToBranch(i: Int) = branchesToIndex.find(_._2 == i).get._1

  def nni(i: Int, which: Boolean): Tree[R, N] = nni(indexToBranch(i), which)

  def nni(b: Branch[N], which: Boolean): Tree[R, N] = {
    require(isInternal(b))
    val u = b.head
    val v = b.tail
    val x = children(u, v).head
    val y = if (which) children(v, u).head else children(v, u).tail.head
    val (branchesp, neighborsp) = Tree.connect(u, y, branchesToIndex(Branch(v, y)))(Tree.disconnect(v, y)(Tree.connect(v, x, branchesToIndex(Branch(u, x)))(Tree.disconnect(u, x)(branchesToIndex, neighbors))))
    copy(branchesToIndex = branchesp, neighbors = neighborsp)
  }

  def children(n: N, p: N): Set[N] = neighbors(n) - p

  def isLeaf(n: N): Boolean = neighbors(n).size == 1

  def isInternal(n: N): Boolean = !isLeaf(n)

  def isPendant(b: Branch[N]): Boolean = !isInternal(b)

  def isInternal(b: Branch[N]): Boolean = isInternal(b.head) && isInternal(b.tail)

  def isInternal(i: Int): Boolean = isInternal(indexToBranch(i))

  lazy val treeLength: R = implicitly[AdditiveMonoid[R]].sum(lengths)

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
    (newick.replace(i, i + 1, "," + taxa(rho).name) ++= ");").toString()

  }

  override lazy val toString = newick

}

object Tree {

  def generateNeighbors[N](nodes: Set[N], branches: Set[Branch[N]]) = branches.foldLeft(nodes.map(_ -> Set[N]()).toMap) { (neighbors, b) =>
    Tree.addNeighbor(b.head, b.tail)(Tree.addNeighbor(b.tail, b.head)(neighbors))
  }

  def connect[@specialized(Int) N](x: N, y: N, i: Int)(bn: (Map[Branch[N], Int], Map[N, Set[N]])) =
    (bn._1 + (Branch(x, y) -> i), addNeighbor(x, y)(addNeighbor(y, x)(bn._2)))

  def disconnect[@specialized(Int) N](x: N, y: N)(bn: (Map[Branch[N], Int], Map[N, Set[N]])) =
    (bn._1 - Branch(x, y), removeNeighbor(x, y)(removeNeighbor(y, x)(bn._2)))

  def addNeighbor[@specialized(Int) N](i: N, e: N) = index[Map[N, Set[N]], N, Set[N]](i).modify(_ + e)

  def removeNeighbor[@specialized(Int) N](i: N, e: N) = index[Map[N, Set[N]], N, Set[N]](i).modify(_ - e)

  def apply[R : AdditiveMonoid, N](nodes: Set[N], lengths: Map[Branch[N], R], taxa: PartialFunction[N, Taxon]): Tree[R, N] = {
    val branchesToIndex = lengths.keys.seq.zipWithIndex.toMap
    Tree(nodes, branchesToIndex, generateNeighbors(nodes, lengths.keySet), IndexedSeq.tabulate(lengths.size)(lengths.map(Function.tupled((k, v) => branchesToIndex(k) -> v))), taxa)
  }

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

}

case class Branch[@specialized(Int) N](head: N, tail: N) {

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
