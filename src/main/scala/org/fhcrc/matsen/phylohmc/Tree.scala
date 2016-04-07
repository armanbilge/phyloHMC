package org.fhcrc.matsen.phylohmc

import spire.algebra.Field

class Tree[R : Field, N](val nodes: Set[N], val branches: Set[Branch[N]])(
                        val neighbors: N => Set[N] = branches.foldLeft(Map[N, Set[N]]().withDefaultValue(Set())) { (neighbors, b) =>
                          neighbors + (b.head -> (neighbors(b.head) + b.tail)) + (b.tail -> (neighbors(b.tail) + b.head))
                        }
) {



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
    case _ => throw new IllegalArgumentException
  }

  override val toString: String = s"($head, $tail)"

  override def equals(that: Any): Boolean = that match {
    case that: Branch => (head == that.head && tail == that.tail) || (head == that.tail && tail == that.head)
    case _ => false
  }

  override val hashCode: Int = Set(head, tail).hashCode()

}
