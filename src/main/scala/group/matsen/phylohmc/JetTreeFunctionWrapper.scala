package group.matsen.phylohmc

import spire.algebra._
import spire.math.{Jet, JetDim}

import scala.reflect.ClassTag

class JetTreeFunctionWrapper[R : Eq : Field : NRoot : IsReal : Trig : ClassTag, N](f: Tree[Jet[R], N] => Jet[R]) extends (Tree[R, N] => (R, IndexedSeq[R])) {

  override def apply(t: Tree[R, N]): (R, IndexedSeq[R]) = {
    implicit val jetDim = JetDim(t.lengths.size)
    val y = f(Tree(t.nodes, t.branchesToIndex, t.neighbors, t.lengths.zipWithIndex.map(Function.tupled(Jet.apply(_, _))), t.taxa))
    (y.real, y.infinitesimal.toIndexedSeq)
  }

}
