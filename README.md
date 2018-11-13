# max-flow

preflow push-relabel algorithmによって, max-flowを求める.

##  開発環境

+ maxOS High Sierra (version 10.13.6)
+ Python3.7



## ファイル



+ pure_*.py

  + optionなしの実装 (Genetic, FIFO, Highest preflow push-relabel algorithm)
  + 単純なアルゴリズムを見る用

  実行方法

  `python3 (ファイル名).py`でusageが表示される.

+ pureなし*.py

  + optionあり(Global labeling, Gap-relabeling, Freeze operation)



  実行方法

  `python3 (ファイル名).py`でusageが表示される.


+ graph_data*.txt
  + グラフデータ
    + nodeノード名 x座標 y座標
    + edge 始点ノード 終点ノード 容量
