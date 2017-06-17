import Criterion.Main

import Fagin.Gff

main :: IO ()
main = defaultMain [
  bgroup "readInt" [ bench "123" $ whnf readInt'' "123" ]
  ]
