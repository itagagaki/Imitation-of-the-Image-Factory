3次元シーンにおける物体の幾何形状、光源、カメラのパラメーターを人間がテキストファイル（イメージ・スコア）に書き表し、
それをもとにコンピューターが映像を生成する。
（当時の）個人所有のマシンでは非力すぎるので、イメージ・スコアを映像にするシステム（イメージ・ファクトリー）を建造する。
……という大村皓一先生（日本におけるコンピューターグラフィックス草創期を代表するお一人）の構想に感化されて、ちょっと書いてみた小さなプログラムたちです。
いつ頃書いたかは覚えていないけど、1983年から1990年までの間のいつかです。
もう残っていないと思っていましたが、バックアップディスクの奥の奥から発掘されたので、とりあえずGitHubで保管しておくことにしました。

`Compiler`ディレクトリは、たぶん、イメージ・スコアをコンパイルして何らかのバイナリーを出力するもの。
構文解析にyaccを使っていたと記憶していましたが、どうやら使っていなかったですね。
当時、簡易的なCコンパイラも作っていて、それにyaccを使っていた記憶とごっちゃになっていたんだと思います。

そして`Renderer`ディレクトリは、たぶん、そのバイナリーからイメージを生成するもの。
アルゴリズムは単純なレイ・トレーシング法です。

なんかメモ書きもちらほらありますが、そのまま保管しておきます。
何か恥ずかしいこと書いて無ければいいのですが(笑)

---

今日（2024年11月13日）これを発掘したことで、大村先生のお名前をGoogleで検索したところ、今年亡くなっていたことを知りました。
先生の大きな功績とお人柄に深く敬意を表しますとともに、謹んで哀悼の意を表します。

---
