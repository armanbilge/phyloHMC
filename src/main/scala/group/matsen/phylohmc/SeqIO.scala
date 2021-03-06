package group.matsen.phylohmc

import scala.io.Source

object SeqIO {

  def parseFasta(file: String): Patterns = {
    val seqs = Source.fromFile(file).mkString.split(">").map(Source.fromString(_).getLines()).drop(1).map(s => (s.next, s.mkString)).toVector
    val taxa = seqs.map(Function.tupled((t, _) => new Taxon(t)))
    seqs(0)._2.indices.map(i =>
      seqs.map(_._2(i)).map {
        case 'A' => 0
        case 'C' => 1
        case 'G' => 2
        case 'T' => 3
        case _ => 4
      }.zip(taxa).map(Function.tupled((i, t) => (t, i))).toMap
    ).foldLeft(Patterns(Patterns.configuration.compact[Pattern]))(_ + _)

  }

}
