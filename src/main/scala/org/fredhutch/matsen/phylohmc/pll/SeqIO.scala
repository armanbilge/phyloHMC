package org.fredhutch.matsen.phylohmc.pll

import org.fredhutch.matsen.phylohmc.Taxon

import scala.collection.Bag
import scala.io.Source

object SeqIO {

  def parseFasta(file: String): Map[Map[Taxon, Char], Int] = {
    val seqs = Source.fromFile(file).mkString.split(">").map(Source.fromString(_).getLines()).drop(1).map(s => (s.next, s.mkString)).toVector
    val taxa = seqs.map(Function.tupled((t, _) => new Taxon(t)))
    val bag = seqs(0)._2.indices.map(i => seqs.map(_._2(i)).zip(taxa).map(Function.tupled((i, t) => (t, i))).toMap).foldLeft(Bag(Bag.configuration.compact[Map[Taxon, Char]]))(_ + _)
    bag.toMap.toMap
  }

}
