package group.matsen.phylohmc

import shapeless.Nat
import spire.algebra.{Field, Trig}
import spire.syntax.vectorSpace._

class HKY[R : Field : Trig](val stationaryDistribution: IndexedSeq[R], val kappa: R) extends SubstitutionModel[R] {

  val freqA = stationaryDistribution(0)
  val freqC = stationaryDistribution(1)
  val freqG = stationaryDistribution(2)
  val freqT = stationaryDistribution(3)

  val freqR = freqA + freqG
  val freqY = freqC + freqT

  val r1 = (1 / freqR) - 1
  val tab1A = freqA * r1

  val tab3A = freqA / freqR
  val tab2A = 1 - tab3A

  val r2 = 1 / r1
  val tab1C = freqC * r2

  val tab3C = freqC / freqY
  val tab2C = 1 - tab3C

  val tab1G = freqG * r1
  val tab3G = tab2A
  val tab2G = tab3A

  val tab1T = freqT * r2

  val tab3T = tab2C
  val tab2T = tab3C

  val beta = 1.0 / (2.0 * (freqR * freqY + kappa * (freqA * freqG + freqC * freqT)))
  val A_R = 1.0 + freqR * (kappa - 1)
  val A_Y = 1.0 + freqY * (kappa - 1)

  override val Q: Matrix[Nat._4, R] =
    beta *: Matrix[Nat._4, R](- freqC - kappa * freqG - freqT, freqA, kappa * freqA, freqA,
      freqC, - freqA - freqG - kappa * freqT, freqC, kappa * freqC,
      kappa * freqG, freqG, - kappa * freqA - freqC - freqT, freqG,
      freqT, kappa * freqT, freqT, - freqA - kappa * freqC - freqG
    )

  override def apply(t: R): Matrix[Nat._4, R] = {

    val xx = - beta * t
    val bbR = Trig[R].exp(xx * A_R)
    val bbY = Trig[R].exp(xx * A_Y)

    val aa = Trig[R].exp(xx)
    val oneminusa = 1 - aa

    val t1Aaa = tab1A * aa
    val t1Gaa = tab1G * aa
    val t1Caa = tab1C * aa
    val t1Taa = tab1T * aa

    Matrix[Nat._4, R](freqA + t1Aaa + (tab2A * bbR), freqA * oneminusa, freqA + t1Aaa - (tab3A * bbR), freqA * oneminusa,
      freqC * oneminusa, freqC + t1Caa + (tab2C * bbY), freqC * oneminusa, freqC + t1Caa - (tab3C * bbY),
      freqG + t1Gaa - (tab3G * bbR), freqG * oneminusa, freqG + t1Gaa + (tab2G * bbR), freqG * oneminusa,
      freqT * oneminusa, freqT + t1Taa - (tab3T * bbY), freqT * oneminusa, freqT + t1Taa + (tab2T * bbY))

  }

}

object HKY {

  def apply[R : Field : Trig](stationaryDistribution: IndexedSeq[R], kappa: R) = new HKY[R](stationaryDistribution, kappa)

}
