# max-flow

preflow push-relabel algorithmによって, max-flowを求める.

##  開発環境

+ Python3.7



## ファイル



+ pure_*.py

  + optionなしの実装 (Genetic, FIFO, Highest preflow push-relabel algorithm)
  + 単純なアルゴリズムを見る用

  実行方法

  `python3 (ファイル名).py`でusageが表示される.

+ pureなし*.py

  + optionあり(Global labeling, Gap-relabeling, Freeze operation)
  + 各optionの[参考](https://qiita.com/nariaki3551/items/65baee3c6ef0a6ffa136)



  実行方法

  `python3 (ファイル名).py`でusageが表示される.


+ graph_data*.txt
  + グラフデータ
    + nodeノード名 x座標 y座標
    + edge 始点ノード 終点ノード 容量
