package group.matsen

import scala.collection.Bag

package object phylohmc {

  type Pattern = Map[Taxon, Int]
  val Pattern = Map
  type Patterns = Bag[Pattern]
  val Patterns = Bag

}
