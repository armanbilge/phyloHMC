name := "phyloHMC"
enablePlugins(GitVersioning)
scalaVersion := "2.11.8"
scalacOptions += "-optimize"
libraryDependencies += "com.chuusai" %% "shapeless" % "2.3.2"
libraryDependencies += "com.github.julien-truffaut" %% "monocle-core" % "1.2.0"
libraryDependencies += "com.lihaoyi" % "ammonite" % "0.8.1" cross CrossVersion.full
libraryDependencies += "io.github.nicolasstucki" %% "multisets" % "0.3"
libraryDependencies += "net.java.dev.jna" % "jna" % "4.2.2"
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1"
libraryDependencies += "org.scala-lang" % "scala-compiler" % "2.11.8"
libraryDependencies += "org.spire-math" %% "spire" % "0.11.0"
